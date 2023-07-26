//
// Created by xiang on 22-12-20.
//
#include <yaml-cpp/yaml.h>
#include <execution>

#include <chrono>

#include "common/lidar_utils.h"
#include "fusion.h"

namespace sad {

/**
 * @description: 高精度融合定位（融合IMU、RTK、激光点云定位功能）
 * @param {string&} config_yaml
 * @return {*}
 */
Fusion::Fusion(const std::string& config_yaml){
    config_yaml_ = config_yaml;
    StaticIMUInit::Options imu_init_options;
    imu_init_options.use_speed_for_static_checking_ = false;  // 本节数据不需要轮速计
    imu_init_ = StaticIMUInit(imu_init_options);

    ndt_pcl_.setResolution(1.0); // PCL NDT

    // 【新增】第七章的NDT
    ndt_ = Ndt3d();
    ndt_.SetResolution(1.0); // 设置分辨率2.0m
}

/**
 * @description: 初始化，读取参数
 * @return {*}
 */
bool Fusion::Init() {
    // 从yaml文件中读取地图原点
    auto yaml = YAML::LoadFile(config_yaml_);
    auto origin_data = yaml["origin"].as<std::vector<double>>();
    map_origin_ = Vec3d(origin_data[0], origin_data[1], origin_data[2]);

    // 读取子地图存储文件夹
    pts_map_path_ = yaml["pts_map_data"].as<std::string>();

    ndt_map_path_10_ = yaml["ndt_map_data_10"].as<std::string>();
    ndt_map_path_5_ = yaml["ndt_map_data_5"].as<std::string>();
    ndt_map_path_4_ = yaml["ndt_map_data_4"].as<std::string>();
    ndt_map_path_3_ = yaml["ndt_map_data_3"].as<std::string>();
    ndt_map_path_1_ = yaml["ndt_map_data_1"].as<std::string>();

    use_pcl_ndt_ = yaml["use_pcl_ndt"].as<bool>();          // 是否使用PCL NDT库

    // 加载地图的索引文件
    LoadMapIndex();

    // lidar和IMU消息同步
    sync_ = std::make_shared<MessageSync>([this](const MeasureGroup& m) { ProcessMeasurements(m); });
    sync_->Init(config_yaml_);

    // lidar和IMU外参
    std::vector<double> ext_t = yaml["mapping"]["extrinsic_T"].as<std::vector<double>>();
    std::vector<double> ext_r = yaml["mapping"]["extrinsic_R"].as<std::vector<double>>();
    Vec3d lidar_T_wrt_IMU = math::VecFromArray(ext_t);
    Mat3d lidar_R_wrt_IMU = math::MatFromArray(ext_r);
    TIL_ = SE3(lidar_R_wrt_IMU, lidar_T_wrt_IMU);

    // ui
    ui_ = std::make_shared<ui::PangolinWindow>();
    ui_->Init();
    ui_->SetCurrentScanSize(50);
    return true;
}

void Fusion::ProcessRTK(GNSSPtr gnss) {
    gnss->utm_pose_.translation() -= map_origin_;  // 减掉地图原点
    last_gnss_ = gnss;
}

void Fusion::ProcessMeasurements(const MeasureGroup& meas) {
    measures_ = meas;

    if (imu_need_init_) {
        TryInitIMU();
        return;
    }

    /// 以下三步与LIO一致，只是align完成地图匹配工作
    if (status_ == Status::WORKING) {
        Predict();
        Undistort();
    } else {
        scan_undistort_ = measures_.lidar_;
    }

    Align();
}

void Fusion::TryInitIMU() {
    for (auto imu : measures_.imu_) {
        imu_init_.AddIMU(*imu);
    }

    if (imu_init_.InitSuccess()) {
        // 读取初始零偏，设置ESKF
        sad::ESKFD::Options options;
        // 噪声由初始化器估计
        // options.gyro_var_ = sqrt(imu_init_.GetCovGyro()[0]);
        // options.acce_var_ = sqrt(imu_init_.GetCovAcce()[0]);
        options.update_bias_acce_ = false;
        options.update_bias_gyro_ = false;
        eskf_.SetInitialConditions(options, imu_init_.GetInitBg(), imu_init_.GetInitBa(), imu_init_.GetGravity());
        imu_need_init_ = false;

        LOG(INFO) << "IMU初始化成功";
    }
}

void Fusion::Predict() {
    imu_states_.clear();
    imu_states_.emplace_back(eskf_.GetNominalState());

    /// 对IMU状态进行预测
    for (auto& imu : measures_.imu_) {
        eskf_.Predict(*imu);
        imu_states_.emplace_back(eskf_.GetNominalState());
    }
}

void Fusion::Undistort() {
    auto cloud = measures_.lidar_;
    auto imu_state = eskf_.GetNominalState();  // 最后时刻的状态
    SE3 T_end = SE3(imu_state.R_, imu_state.p_);

    /// 将所有点转到最后时刻状态上
    std::for_each(std::execution::par_unseq, cloud->points.begin(), cloud->points.end(), [&](auto& pt) {
        SE3 Ti = T_end;
        NavStated match;

        // 根据pt.time查找时间，pt.time是该点打到的时间与雷达开始时间之差，单位为毫秒
        math::PoseInterp<NavStated>(
            measures_.lidar_begin_time_ + pt.time * 1e-3, imu_states_, [](const NavStated& s) { return s.timestamp_; },
            [](const NavStated& s) { return s.GetSE3(); }, Ti, match);

        Vec3d pi = ToVec3d(pt);
        Vec3d p_compensate = TIL_.inverse() * T_end.inverse() * Ti * TIL_ * pi;

        pt.x = p_compensate(0);
        pt.y = p_compensate(1);
        pt.z = p_compensate(2);
    });
    scan_undistort_ = cloud;
}

void Fusion::Align() {
    FullCloudPtr scan_undistort_trans(new FullPointCloudType);
    pcl::transformPointCloud(*scan_undistort_, *scan_undistort_trans, TIL_.matrix());
    scan_undistort_ = scan_undistort_trans;
    current_scan_ = ConvertToCloud<FullPointType>(scan_undistort_);
    current_scan_ = VoxelCloud(current_scan_, 0.5);

    if (status_ == Status::WAITING_FOR_RTK) {
        // 若存在最近的RTK信号，则尝试初始化
        if (last_gnss_ != nullptr) {
            if (SearchRTK()) {
                // status_ == Status::WORKING;
                ui_->UpdateScan(current_scan_, eskf_.GetNominalSE3());
                ui_->UpdateNavState(eskf_.GetNominalState());
            }
        }
    } else {
        LidarLocalization();
        ui_->UpdateScan(current_scan_, eskf_.GetNominalSE3());
        ui_->UpdateNavState(eskf_.GetNominalState());
    }
}

/**
 * @description: RTK初始搜索
 * 
 *  当系统未收到第一个有效的RTK信号时，无法直到自车在地图中的位置，也就无法进行点云定位
 *  因为UCLT数据集中的RTK是单天线方案，只有位置，没有旋转信息，因此如果某个时刻收到首个RTK信号，就在它周边进行网格搜索（多分辨率NDT匹配），用于搜索车辆的初始航向角。
 * 
 * 
 * 
 * @return {*}
 */
bool Fusion::SearchRTK() {
    // 是否初始化失败过
    if (init_has_failed_) {
        // 若失败过，就跳过该RTK位置
        if ((last_gnss_->utm_pose_.translation() - last_searched_pos_.translation()).norm() < 20.0) {
            LOG(INFO) << "skip this position";
            return false;
        }
    }

    if (use_pcl_ndt_) 
        // 加载上一个RTK位置附近的地图
        LoadMap(last_gnss_->utm_pose_);
    else
        // TODO: 加载上一个RTK位置附近的NDT 体素地图
        LoadNdtMap(last_gnss_->utm_pose_);
    
    // 由于RTK不带姿态，我们必须先搜索一定的角度范围
    std::vector<GridSearchResult> search_poses;
    /// 由于RTK不带角度，这里按固定步长扫描RTK角度
    double grid_ang_range = 360.0, grid_ang_step = 10;  // 角度搜索范围与步长
    for (double ang = 0; ang < grid_ang_range; ang += grid_ang_step) {
        SE3 pose(SO3::rotZ(ang * math::kDEG2RAD), Vec3d(0, 0, 0) + last_gnss_->utm_pose_.translation());
        GridSearchResult gr;
        gr.pose_ = pose;
        search_poses.emplace_back(gr);
    }

    LOG(INFO) << "grid search poses: " << search_poses.size();
    std::for_each(std::execution::par_unseq, search_poses.begin(), search_poses.end(),
                  [this](GridSearchResult& gr) { 
                    AlignForGrid(gr); 
                });

    // 选择最优的匹配结果
    auto max_ele = std::max_element(search_poses.begin(), search_poses.end(),
                                    [](const auto& g1, const auto& g2) { return g1.score_ < g2.score_; });
    LOG(INFO) << "max score: " << max_ele->score_ << ", pose: \n" << max_ele->result_pose_.matrix();
    if (max_ele->score_ > rtk_search_min_score_) {
        LOG(INFO) << "初始化成功, score: " << max_ele->score_ << ">" << rtk_search_min_score_;
        status_ = Status::WORKING;

        /// 重置滤波器状态
        auto state = eskf_.GetNominalState();
        state.R_ = max_ele->result_pose_.so3();
        state.p_ = max_ele->result_pose_.translation();
        state.v_.setZero();
        eskf_.SetX(state, eskf_.GetGravity());

        ESKFD::Mat18T cov;
        cov = ESKFD::Mat18T::Identity() * 1e-4;
        cov.block<12, 12>(6, 6) = Eigen::Matrix<double, 12, 12>::Identity() * 1e-6;
        eskf_.SetCov(cov);

        return true;
    }

    init_has_failed_ = true;
    last_searched_pos_ = last_gnss_->utm_pose_;
    return false;
}

/**
 * @description: 多分辨率NDT匹配，用于初始RTK角度搜索
 * @param {GridSearchResult&} gr
 * @return {*}
 */
void Fusion::AlignForGrid(sad::Fusion::GridSearchResult& gr) {
    if (use_pcl_ndt_) {
        /// 多分辨率
        pcl::NormalDistributionsTransform<PointType, PointType> ndt_pcl;
        ndt_pcl.setTransformationEpsilon(0.05);
        ndt_pcl.setStepSize(0.7);
        ndt_pcl.setMaximumIterations(40);

        ndt_pcl.setInputSource(current_scan_);

        // Ndt3d ndt;
        // ndt.SetSource(current_scan_); // 待检索的当前帧点云

        auto map = ref_cloud_;

        CloudPtr output(new PointCloudType);
        // 分别在10m, 5m, 4m, 3m的分辨率下进行NDT配准
        std::vector<double> res{10.0, 5.0, 4.0, 3.0}; // 四种体素分辨率
        Mat4f T = gr.pose_.matrix().cast<float>();
        SE3 pose = gr.pose_;
        for (auto& r : res) {
            auto rough_map = VoxelCloud(map, r * 0.1);
            ndt_pcl.setInputTarget(rough_map);
            ndt_pcl.setResolution(r);
            ndt_pcl.align(*output, T);
            // 将上一个粗配准结果代入下一次配准中
            T = ndt_pcl.getFinalTransformation(); 
            // ndt.SetResolution(r);   // 设置分辨率 
            // ndt.SetTarget(rough_map);
            // ndt.AlignNdt(pose);

        }
        // 最后获取3米栅格分辨率下的的NDT匹配结果作为分值判定
        gr.score_ = ndt_pcl.getTransformationProbability();
        // gr.score_ = ndt.GetNdtMatchingScore();
        LOG(INFO) << "gr.score_ = " << gr.score_;
        // 获取最后一次配准得到的位姿
        gr.result_pose_ = Mat4ToSE3(ndt_pcl.getFinalTransformation());
        // gr.result_pose_ = pose;
    }
    else {
        Ndt3d::Options options;
        Ndt3d ndt(options);
        ndt.SetSource(current_scan_); // 待检索的当前帧点云

        // // 分别在10m, 5m, 4m, 3m的四种分辨率下进行NDT配准
        std::vector<double> res{10.0, 5.0, 4.0, 3.0}; // split_map划分的时候就按这三种分辨率来构建NDT地图
        SE3 pose = gr.pose_;
        for (auto& r : res) {
            LOG(INFO) << "r = " << r;
            if (r == 10.0) 
                ndt.SetNdtVoxel(ndt_map_grids_10_);
            else if (r == 5.0) 
                ndt.SetNdtVoxel(ndt_map_grids_5_);
            else if (r == 4.0) 
                ndt.SetNdtVoxel(ndt_map_grids_4_);
            else 
                ndt.SetNdtVoxel(ndt_map_grids_3_);
            
            ndt.SetResolution(r);   // 设置分辨率 
            ndt.AlignNdt(pose);     // 将上一个粗配准结果代入下一次配准中
        }
        // 最后获取2米栅格分辨率下的的NDT匹配结果作为分值判定
        gr.score_ = ndt.GetNdtMatchingScore();
        LOG(INFO) << "gr.score_ = " << gr.score_;
        // 获取最后一次配准得到的位姿
        gr.result_pose_ = pose;
    }
}

/**
 * @description: 雷达定位：新增第七章3D NDT的实现（匹配得分、加载地图区块内的NDT体素） 
 * @return {*}
 */
bool Fusion::LidarLocalization() {
    // 获取ESKF名义状态递推得到的预测位姿
    SE3 pred = eskf_.GetNominalSE3();

    SE3 pose;
    if (use_pcl_ndt_) {
        // 加载pred位置处周围相邻9个地图区块，拼接后形成局部地图，设置为NDT的目标点云
        auto t1 = std::chrono::steady_clock::now();
        LoadMap(pred);  
        auto t2 = std::chrono::steady_clock::now();
        auto time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1).count() * 1000;
        LOG(INFO) << "time used of ndt_pcl.LoadMap() = " << time_used;

        ndt_pcl_.setInputCloud(current_scan_);
        CloudPtr output(new PointCloudType);

        auto t3 = std::chrono::steady_clock::now();
        ndt_pcl_.align(*output, pred.matrix().cast<float>());
        auto t4 = std::chrono::steady_clock::now();
        time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t4-t3).count() * 1000;
        LOG(INFO) << "time used of ndt_pcl.align() = " << time_used;

        pose = Mat4ToSE3(ndt_pcl_.getFinalTransformation());

        LOG(INFO) << "ndt_pcl lidar loc score: " << ndt_pcl_.getTransformationProbability();
    }
    else {
        pose = pred;
        // 加载pred位置处周围相邻9个地图区块，拼接后形成局部地图，设置为NDT的目标点云
        auto t1 = std::chrono::steady_clock::now();
        LoadNdtMap(pred); // 【新增】加载指定位姿处周围地图区块内所有体素的均值和协方差信息
        auto t2 = std::chrono::steady_clock::now();
        auto time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1).count() * 1000;
        LOG(INFO) << "time used of ndt3d.LoadMap() = " << time_used;

        ndt_.SetSource(current_scan_);  //【新增】设置当前帧扫描点云为NDT的源点云
        auto t3 = std::chrono::steady_clock::now();
        ndt_.AlignNdt(pose);            //【新增】NDT配准，并计算匹配得分
        auto t4 = std::chrono::steady_clock::now();
        time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t4-t3).count() * 1000;
        LOG(INFO) << "time used of ndt3d.align() = " << time_used;
        
        LOG(INFO) << "ndt3d pts map lidar loc score: " << ndt_.GetNdtMatchingScore();
    }
    
    eskf_.ObserveSE3(pose, 1e-1, 1e-2); // 将NDT配准得到的位姿作为ESKF的观测值，更新ESKF状态
    return true;
}

/**
 * @description: 根据给定的位姿的平移信息，加载、卸载必要的地图区块
 * @param {SE3&} pose
 * @return {*}
 */
void Fusion::LoadMap(const SE3& pose) {
    // 将给定的上一次的RTK位置信息转换为地图索引
    int gx = floor((pose.translation().x() - 50.0) / 100); // 100米乘100米的地图区块
    int gy = floor((pose.translation().y() - 50.0) / 100);
    Vec2i key(gx, gy);

    // 一个区域的周边地图，我们认为9个就够了
    std::set<Vec2i, less_vec<2>> surrounding_index{
        key + Vec2i( 0, 0), key + Vec2i(-1,  0), key + Vec2i(-1, -1), 
        key + Vec2i(-1, 1), key + Vec2i( 0, -1), key + Vec2i( 0,  1), 
        key + Vec2i( 1, 0), key + Vec2i( 1, -1), key + Vec2i( 1,  1),
    };

    bool map_data_changed = false;
    int cnt_new_loaded = 0, cnt_unload = 0;
    // 加载包含给定位置在内的9个相邻地图区块
    for (auto& k : surrounding_index) {
        if (pts_map_data_index_.find(k) == pts_map_data_index_.end()) 
            // 没找到索引，该地图数据不存在
            continue;

        // 通过地图索引从map_data_中查找是否已经加载了这个区块
        if (map_data_.find(k) == map_data_.end()) {
            // 还未加载过，则加载该地图区块
            CloudPtr cloud(new PointCloudType);
            // 读取以地图区块索引命名的pcd点云数据
            pcl::io::loadPCDFile(pts_map_path_ + std::to_string(k[0]) + "_" + std::to_string(k[1]) + ".pcd", *cloud);
            
            map_data_.emplace(k, cloud);    // 将该地图区块的索引以及点云数据加入到map_data_容器中
            map_data_changed = true;        // 标记地图数据已经改变
            cnt_new_loaded++;               // 计数加一
        }
    }

    // 遍历map_data_容器内的地图区块
    // 卸载不需要的区域，这个稍微加大一点，不需要频繁卸载
    for (auto iter = map_data_.begin(); iter != map_data_.end();) {
        // 判断当前地图区块与给定RTK位置的距离是否大于3米
        if ((iter->first - key).cast<float>().norm() > 3.0) {
            // 若超过3米，太远了，则从map容器中卸载该地图区块
            iter = map_data_.erase(iter);
            cnt_unload++;                   // 计数加一
            map_data_changed = true;        // 标记地图数据已经改变
        } else 
            iter++;
    }

    LOG(INFO) << "new loaded: " << cnt_new_loaded << ", unload: " << cnt_unload;
    // 若地图数据发生改变
    if (map_data_changed) {
        // 重置目标点云 rebuild ndt target map
        ref_cloud_.reset(new PointCloudType);
        // 将map_data_容器内的所有地图区块拼接到ref_cloud_参考目标点云中
        for (auto& mp : map_data_) 
            *ref_cloud_ += *mp.second;

        LOG(INFO) << "rebuild global cloud, grids: " << map_data_.size();

        // 将拼接完成地图区块点云设置NDT的目标点云
        ndt_pcl_.setInputTarget(ref_cloud_); 
        // ndt_.SetTarget(ref_cloud_);
    }

    // 动态更新地图区块的点云显示
    ui_->UpdatePointCloudGlobal(map_data_);
}

// 加载NDT地图区块内所有体素的均值和协方差信息
void Fusion::LoadNdtMap(const SE3& pose) {
    // 将给定的上一次的RTK位置信息转换为地图索引
    int gx = floor((pose.translation().x() - 50.0) / 100); // 100米乘100米的地图区块
    int gy = floor((pose.translation().y() - 50.0) / 100);
    Vec2i key(gx, gy);

    // 一个区域的周边地图，我们认为9个就够了
    std::set<Vec2i, less_vec<2>> surrounding_index{
        key + Vec2i( 0, 0), key + Vec2i(-1,  0), key + Vec2i(-1, -1), 
        key + Vec2i(-1, 1), key + Vec2i( 0, -1), key + Vec2i( 0,  1), 
        key + Vec2i( 1, 0), key + Vec2i( 1, -1), key + Vec2i( 1,  1),
    };

    bool map_data_changed = false;
    bool ndt_map_data_changed_10 = false;
    bool ndt_map_data_changed_5 = false;
    bool ndt_map_data_changed_2 = false;
    bool ndt_map_data_changed_1 = false;
    int cnt_new_loaded = 0, cnt_unload = 0;

    // lambda表达式
    auto load_ndt_map_k = [](std::string ndt_map_path, Vec2i k, std::map<Vec2i, NdtmapVoxel, less_vec<2>> &ndt_map_data, bool &changed){
        std::ifstream fin(ndt_map_path + std::to_string(k[0]) + "_" + std::to_string(k[1]) + ".txt");
        if (ndt_map_data.find(k) == ndt_map_data.end()) {
            LOG(INFO) << "clark";
            NdtmapVoxel ndt_map_grids;      
            // 【新增】临时变量
            int x, y, z;
            double mu_x, mu_y, mu_z;
            double info_00, info_01, info_02, info_10, info_11, info_12, info_20, info_21, info_22;
            Eigen::Matrix<int, 3, 1> key_ndt;
            Vec3d mu;
            Mat3d info;
            while (!fin.eof()) {
                fin  >> x >> y >> z >> mu_x >> mu_y >> mu_z 
                    >> info_00 >> info_01 >> info_02 
                    >> info_10 >> info_11 >> info_12 
                    >> info_20 >> info_21 >> info_22;
                key_ndt << x, y, z; 
                mu = Vec3d(mu_x, mu_y, mu_z);
                info << info_00, info_01, info_02, info_10, info_11, info_12, info_20, info_21, info_22;
                if (ndt_map_grids.find(key_ndt) == ndt_map_grids.end()) {
                    Ndt3d::VoxelData voxelData;
                    voxelData.mu_ = mu;
                    voxelData.info_ = info;
                    // LOG(INFO) << "insert one ndt voxel";
                    ndt_map_grids.insert(std::make_pair(key_ndt, voxelData)); // 若不存在，则插入该栅格
                }
            }   
            ndt_map_data.emplace(k, ndt_map_grids);  
            changed = true;        // 标记地图数据已经改变
        }
        fin.close();
    };

    auto iterate_Erase_NdtVoxel = [](std::map<Vec2i, NdtmapVoxel, less_vec<2>> &ndt_map_data, Vec2i k, bool &changed){
        for (auto iter = ndt_map_data.begin(); iter != ndt_map_data.end();) {
        // 判断当前地图区块与给定RTK位置的距离是否大于3米
        if ((iter->first - k).cast<float>().norm() > 3.0) {
            // 若超过3米，太远了，则从ndt map容器中卸载该地图区块内所有的ndt体素信息
            iter = ndt_map_data.erase(iter);
            changed = true;
        } else 
            iter++; // 取下一个
    }
    };

    // 加载包含给定位置在内的9个相邻地图区块
    for (auto& k : surrounding_index) {
        if (pts_map_data_index_.find(k) != pts_map_data_index_.end()) {
            // 通过地图索引从map_data_中查找是否已经加载了这个区块
            if (map_data_.find(k) == map_data_.end()) {
                // 还未加载过，则加载该地图区块
                CloudPtr cloud(new PointCloudType);
                // 读取以地图区块索引命名的pcd点云数据
                pcl::io::loadPCDFile(pts_map_path_ + std::to_string(k[0]) + "_" + std::to_string(k[1]) + ".pcd", *cloud);
                
                map_data_.emplace(k, cloud);    // 将该地图区块的索引以及点云数据加入到map_data_容器中
                map_data_changed = true;        // 标记地图数据已经改变
            }
        }
            
        // TODO: 【新增】加载每个地图区块内所有NDT体素块的均值和协方差信息
        // 可能保存于以地图区块索引命名德txt文件或者二进制文件中
        // 加载以地图区块索引命名的ndt体素数据，每个地图区块内，有很多ndt体素，地图区块100x100 m^2，ndt体素分辨率1x1x1 m^3
        if (ndt_map_data_index_10_.find(k) != ndt_map_data_index_10_.end()) 
            load_ndt_map_k(ndt_map_path_10_, k, ndt_map_data_10_, ndt_map_data_changed_10);

        if (ndt_map_data_index_5_.find(k) != ndt_map_data_index_5_.end()) 
            load_ndt_map_k(ndt_map_path_5_, k, ndt_map_data_5_, ndt_map_data_changed_5);

        if (ndt_map_data_index_4_.find(k) != ndt_map_data_index_4_.end())
            load_ndt_map_k(ndt_map_path_4_, k, ndt_map_data_4_, ndt_map_data_changed_4);
        
        if (ndt_map_data_index_3_.find(k) != ndt_map_data_index_3_.end())
            load_ndt_map_k(ndt_map_path_3_, k, ndt_map_data_3_, ndt_map_data_changed_3);

        if (ndt_map_data_index_1_.find(k) != ndt_map_data_index_1_.end())
            load_ndt_map_k(ndt_map_path_1_, k, ndt_map_data_1_, ndt_map_data_changed_1);
    }

    // 遍历map_data_容器内的地图区块
    // 卸载不需要的区域，这个稍微加大一点，不需要频繁卸载
    for (auto iter = map_data_.begin(); iter != map_data_.end();) {
        // 判断当前地图区块与给定RTK位置的距离是否大于3米
        if ((iter->first - key).cast<float>().norm() > 3.0) {
            // 若超过3米，太远了，则从map容器中卸载该地图区块
            iter = map_data_.erase(iter);
            map_data_changed = true;        // 标记地图数据已经改变
        } else 
            iter++;
    }

    // 遍历ndt_map_grids_容器内的地图区块
    // 卸载不需要的区域，这个稍微加大一点，不需要频繁卸载
    iterate_Erase_NdtVoxel(ndt_map_data_10_, key, ndt_map_data_changed_10);
    iterate_Erase_NdtVoxel(ndt_map_data_5_, key, ndt_map_data_changed_5);
    iterate_Erase_NdtVoxel(ndt_map_data_4_, key, ndt_map_data_changed_4);
    iterate_Erase_NdtVoxel(ndt_map_data_3_, key, ndt_map_data_changed_3);
    iterate_Erase_NdtVoxel(ndt_map_data_1_, key, ndt_map_data_changed_1);

    // 若地图数据发生改变
    if (ndt_map_data_changed_10) {
        for (auto& nmp : ndt_map_data_10_) 
            ndt_map_grids_10_.insert(nmp.second.begin(), nmp.second.end());
    }

    if (ndt_map_data_changed_5) {
        for (auto& nmp : ndt_map_data_5_)
            ndt_map_grids_5_.insert(nmp.second.begin(), nmp.second.end());
    }

    if (ndt_map_data_changed_4) {
        for (auto& nmp : ndt_map_data_4_)
            ndt_map_grids_4_.insert(nmp.second.begin(), nmp.second.end());
    }

    if (ndt_map_data_changed_3) {
        for (auto& nmp : ndt_map_data_3_)
            ndt_map_grids_3_.insert(nmp.second.begin(), nmp.second.end());
    }

    if (ndt_map_data_changed_1) {
        LOG(INFO) << " ndt_map_data_changed_1: " << ndt_map_data_changed_1 << " ndt_map_grids_1_.size(): "<<ndt_map_grids_1_.size()<< " ndt_map_data_1_.size(): "<<ndt_map_data_1_.size();
        for (auto& nmp : ndt_map_data_1_)
            ndt_map_grids_1_.insert(nmp.second.begin(), nmp.second.end());
        LOG(INFO) << "ndt_map_grids_1_.size(): " << ndt_map_grids_1_.size();
        ndt_.SetNdtVoxel(ndt_map_grids_1_);
    }

    // 若地图数据发生改变
    if (map_data_changed) {
        // 重置目标点云 rebuild ndt target map
        ref_cloud_.reset(new PointCloudType);
        // 将map_data_容器内的所有地图区块拼接到ref_cloud_参考目标点云中
        for (auto& mp : map_data_) 
            *ref_cloud_ += *mp.second;

        LOG(INFO) << "rebuild global cloud, grids: " << map_data_.size();

        // 将拼接完成地图区块点云设置NDT的目标点云
        // ndt_pcl_.setInputTarget(ref_cloud_); 
        // ndt_.SetTarget(ref_cloud_);
    }

    // 动态更新地图区块的点云显示
    ui_->UpdatePointCloudGlobal(map_data_);
}

/**
 * @description: 
 * 加载地图索引文件
    -3 0
    -2 -3
    -2 -2
    -2 -1
    -2 0
    -1 -2
    -1 -1
    -1 0
    0 -1
 * @return {*}
 */
void Fusion::LoadMapIndex() {

    std::ifstream fin_pts(pts_map_path_ + "/pts_map_index.txt");
    // 【新增】读取NDT地图区块索引，与上面点云地图的区块类似，但可能不完全一样
    std::ifstream fin_ndt_10(ndt_map_path_10_ + "/ndt_map_index.txt");
    std::ifstream fin_ndt_5(ndt_map_path_5_ + "/ndt_map_index.txt");
    std::ifstream fin_ndt_4(ndt_map_path_4_ + "/ndt_map_index.txt");
    std::ifstream fin_ndt_3(ndt_map_path_3_ + "/ndt_map_index.txt");
    std::ifstream fin_ndt_1(ndt_map_path_1_ + "/ndt_map_index.txt");
    while (!fin_pts.eof()) {
        int x, y;   // 块地图的索引
        fin_pts >> x >> y;
        pts_map_data_index_.emplace(Vec2i(x, y));
    }
    while (!fin_ndt_10.eof()) {
        int x, y;   // 块地图的索引
        fin_ndt_10 >> x >> y;
        ndt_map_data_index_10_.emplace(Vec2i(x, y));
    }
    while (!fin_ndt_5.eof()) {
        int x, y;   // 块地图的索引
        fin_ndt_5 >> x >> y;
        ndt_map_data_index_5_.emplace(Vec2i(x, y));
    }
    while (!fin_ndt_4.eof()) {
        int x, y;   // 块地图的索引
        fin_ndt_4 >> x >> y;
        ndt_map_data_index_4_.emplace(Vec2i(x, y));
    }
    while (!fin_ndt_3.eof()) {
        int x, y;   // 块地图的索引
        fin_ndt_3 >> x >> y;
        ndt_map_data_index_3_.emplace(Vec2i(x, y));
    }
    while (!fin_ndt_1.eof()) {
        int x, y;   // 块地图的索引
        fin_ndt_1 >> x >> y;
        ndt_map_data_index_1_.emplace(Vec2i(x, y));
    }
    fin_pts.close();
    fin_ndt_10.close();
    fin_ndt_5.close();
    fin_ndt_4.close();
    fin_ndt_3.close();
    fin_ndt_1.close();
}

/**
 * @description: 处理IMU数据
 * @param {IMUPtr} imu
 * @return {*}
 */
void Fusion::ProcessIMU(IMUPtr imu) { 
    sync_->ProcessIMU(imu); 
}

/**
 * @description: 处理点云数据
 * @param {Ptr} cloud
 * @return {*}
 */
void Fusion::ProcessPointCloud(sensor_msgs::PointCloud2::Ptr cloud) {
    sync_->ProcessCloud(cloud); 
}

}  // namespace sad