#include <pcl/common/transforms.h>
#include <yaml-cpp/yaml.h>
#include <execution>
#include <fstream>

#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/robust_kernel.h>
#include <g2o/core/sparse_block_matrix.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>

// #include "ch4/g2o_types_preinteg.h"
#include "common/g2o_types.h"

#include "common/lidar_utils.h"
#include "common/point_cloud_utils.h"
#include "common/timer/timer.h"

#include "lio_iekf.h"

namespace sad {

/**
 * @description: 构造函数
 * @param {Options} options
 * @return {*}
 */
LioIEKF::LioIEKF(Options options) : options_(options) {
    StaticIMUInit::Options imu_init_options;
    imu_init_options.use_speed_for_static_checking_ = false;  // 本节数据不需要轮速计
    imu_init_ = StaticIMUInit(imu_init_options);
}

/**
 * @description: 初始化
 * @param {string} &config_yaml 配置文件路径
 * @return {*}
 */
bool LioIEKF::Init(const std::string &config_yaml) {
    // 判断配置文件是否读取成功
    if (!LoadFromYAML(config_yaml)) {
        LOG(INFO) << "init failed.";
        return false;
    }

    if (options_.with_ui_) {
        ui_ = std::make_shared<ui::PangolinWindow>();
        ui_->Init();
    }

    return true;
}

/**
 * @description: 处理同步之后的IMU和雷达数据
 * @param {MeasureGroup} &meas
 * @return {*}
 */
void LioIEKF::ProcessMeasurements(const MeasureGroup &meas) {
    LOG(INFO) << "call meas, imu: " << meas.imu_.size() << ", lidar pts: " << meas.lidar_->size();
    measures_ = meas;

    if (imu_need_init_) {
        // 初始化IMU系统
        TryInitIMU();
        return;
    }

    // 利用IMU数据进行状态预测
    Predict();

    // 利用同步的IMU数据对点云进行去畸变，
    // ESKF预测得到的每个IMU数据的名义状态，找到对应的IMU数据，
    // 利用SLERP球面线性插值得到旋转，利用线性插值得到位移
    Undistort();

    // 配准
    Align();
}

/**
 * @description: 从yaml文件中读取参数
 * @param {string} &yaml_file
 * @return {*}
 */
bool LioIEKF::LoadFromYAML(const std::string &yaml_file) {
    // 给消息同步器添加回调函数
    sync_ = std::make_shared<MessageSync>([this](const MeasureGroup &m) {  // lamgda函数赋值给std::function进行回调函数包装
                                                    ProcessMeasurements(m);  // 处理同步好的IMU和雷达数据
                                                });
    // 用yaml文件初始化消息同步器
    sync_->Init(yaml_file);

    /// 自身参数主要是雷达与IMU外参
    auto yaml = YAML::LoadFile(yaml_file);
    std::vector<double> ext_t = yaml["mapping"]["extrinsic_T"].as<std::vector<double>>();
    std::vector<double> ext_r = yaml["mapping"]["extrinsic_R"].as<std::vector<double>>();

    Vec3d lidar_T_wrt_IMU = math::VecFromArray(ext_t);
    Mat3d lidar_R_wrt_IMU = math::MatFromArray(ext_r);
    // 由旋转矩阵R和平移向量t构造SE3
    TIL_ = SE3(lidar_R_wrt_IMU, lidar_T_wrt_IMU);
    return true;
}

/**
 * @description: 局部地图处理
 * @return {*}
 */
void LioIEKF::MapIncremental() {
    PointVec points_to_add;
    PointVec point_no_need_downsample;

    int cur_pts = scan_down_body_->size();
    points_to_add.reserve(cur_pts);
    point_no_need_downsample.reserve(cur_pts);

    std::vector<size_t> index(cur_pts);
    for (size_t i = 0; i < cur_pts; ++i) 
        index[i] = i;

    // 并发处理
    std::for_each(  std::execution::unseq, 
                    index.begin(), index.end(), 
                    [&](const size_t &i) {
                        /* transform to world frame */
                        // 雷达系转换到世界系
                        // PointBodyToWorld(&(scan_down_body_->points[i]), &(scan_down_world_->points[i]));

                        /* decide if need add to map */
                        // 判断是否需要加入到局部地图中
                        PointType &point_world = scan_down_world_->points[i];

                        // 判断第i个点的近邻点集是否为空
                        if (!nearest_points_[i].empty() && flg_ESKF_inited_) {
                            // 取出第i个点的近邻点集
                            const PointVec &points_near = nearest_points_[i];

                            // 计算中心坐标
                            Eigen::Vector3f center = ((point_world.getVector3fMap() / filter_size_map_min_).array().floor() + 0.5) * filter_size_map_min_;

                            // 计算第i个点到中心点的L1距离
                            Eigen::Vector3f dis_2_center = points_near[0].getVector3fMap() - center;

                            // 判断距离是否大于阈值
                            if (fabs(dis_2_center.x()) > 0.5 * filter_size_map_min_ &&
                                fabs(dis_2_center.y()) > 0.5 * filter_size_map_min_ &&
                                fabs(dis_2_center.z()) > 0.5 * filter_size_map_min_) {
                                // 若是，则加入到无需降采样点集中
                                point_no_need_downsample.emplace_back(point_world);
                                return; // 程序返回？因为这里是lambda函数内部，所以返回的是lambda函数，而不是MapIncremental函数
                            }

                            // 此时，标记改为需要增加
                            bool need_add = true;
                            // 计算第i个点到中心点的L2距离
                            float dist = math::calc_dist(point_world.getVector3fMap(), center); // 【在math_utils.h】中添加了两个函数实现
                            // 判断近邻点数是否多于5个
                            if (points_near.size() >= options_.NUM_MATCH_POINTS) {
                                // 遍历所有近邻点
                                for (int readd_i = 0; readd_i < options_.NUM_MATCH_POINTS; readd_i++) {
                                    // 判断这些近邻点距离中心点的距离是否小于阈值
                                    if (math::calc_dist(points_near[readd_i].getVector3fMap(), center) < dist + 1e-6) {
                                        need_add = false; // 只要有一个距离很小的，就不需要增加了，直接跳出循环
                                        break;
                                    }
                                }
                            }
                            // 判断是否需要增加
                            if (need_add) 
                                // 加入到需要增加的点集中
                                points_to_add.emplace_back(point_world);
                        } else 
                            points_to_add.emplace_back(point_world);
                    });

    ivox_->AddPoints(points_to_add);
    ivox_->AddPoints(point_no_need_downsample);
}

/**
 * @description: 计算残差和雅克比矩阵【新增】
 * @param {SE3&} input_pose
 * @param {Mat18d&} HTVH
 * @param {Vec18d&} HTVr
 * @return {*}
 */
void LioIEKF::ComputeResidualAndJacobians_P2Plane(const SE3& input_pose, Mat18d& HTVH, Vec18d& HTVr) {
    LOG(INFO) << "aligning with point to plane";
    int cnt_pts = scan_down_body_->size();
    // assert(target_ != nullptr && source_ != nullptr);

    // 大部分流程和前面的AlignP2Plane()是一样的，只是会把z, H, R三者抛出去，而非自己处理

    // 输入位姿，来自ESKF的Predict()函数预测得到的名义旋转R_、名义位移T_
    SE3 pose = input_pose;
    // if (!options_.use_initial_translation_) 
    //     pose.translation() = target_center_ - source_center_;  // 设置平移初始值

    // 初始化索引，0，1，2，3，4。。。
    std::vector<int> index(cnt_pts);
    for (int i = 0; i < index.size(); ++i) 
        index[i] = i;

    std::vector<bool> effect_pts(index.size(), false);                    // 用于标记有效点
    std::vector<Eigen::Matrix<double, 1, 18>> jacobians(index.size());    // 用于存储雅可比矩阵
    std::vector<double> errors(index.size());                             // 用于存储残差
    // std::vector<Mat3d> infos(index.size());                            // 用于存储信息矩阵

    // gauss-newton 迭代
    // 最近邻，可以并发
    std::for_each(  std::execution::par_unseq, 
                    index.begin(), index.end(), 
                    [&](int idx) {
                        // 并发遍历到点云中的某个点，不是按顺序遍历的
                        auto q = ToVec3d(scan_down_body_->points[idx]);

                        PointType &point_world = scan_down_world_->points[idx];
                        
                        // 利用ESKF预测的名义位姿对该点进行转换
                        // 雷达系转换到IMU系：P_I = R_IL * P_L + T_IL
                        Vec3d qs = pose * q;  // R * q + t

                        point_world = ToPointType(qs);

                        auto &points_near = nearest_points_[idx];
                        // kd树中查找转换后点的5个最近邻
                        // std::vector<int> nn;
                        // kdtree_->GetClosestPoint(ToPointType(qs), nn, 5);  
                        
                        // [新增]kdtree search五个近邻点，存在nearest_points_[i]中
                        ivox_->GetClosestPoint(point_world, points_near, options_.NUM_MATCH_POINTS);

                        // 判断查找到近邻点数是否多于3个，平面方程拟合，a*x+b*y+c*z+d=0，最少需要4个点才能拟合出平面系数
                        if (points_near.size() > options_.MIN_NUM_MATCH_POINTS) {
                            std::vector<Vec3d> nn_eigen;
                            // 遍历近邻点集
                            for (int i = 0; i < points_near.size(); ++i) 
                                // 将近邻点转换为Vec3d类型存储
                                nn_eigen.emplace_back(ToVec3d(points_near[i]));
                            
                            Vec4d n;
                            // 对这几个近邻点执行平面拟合，平面系数a,b,c,d存储在四维向量n中
                            if (!math::FitPlane(nn_eigen, n)) {
                                effect_pts[idx] = false; // 平面拟合失败，标记为无效点
                                return;
                            }

                            // 计算点到平面的距离
                            double dis = n.head<3>().dot(qs) + n[3]; 
                            // 添加阈值检查判断拟合出的平面是否合理
                            if (fabs(dis) > options_.max_plane_distance_) {
                                // 点离的太远了不要
                                effect_pts[idx] = false;
                                return;
                            }

                            // 构建雅可比矩阵，对应公式（7.7）
                            Eigen::Matrix<double, 1, 18> J;
                            J.setZero(); // 其它四项1x3的块矩阵均为零矩阵
                            J.block<1, 3>(0, 0) = n.head<3>().transpose();
                            J.block<1, 3>(0, 3) = -n.head<3>().transpose() * pose.so3().matrix() * SO3::hat(q);

                            jacobians[idx] = J;
                            errors[idx] = dis;
                            effect_pts[idx] = true; // 标记为有效点
                        } else 
                            effect_pts[idx] = false;
                    });

    // 累加Hessian和error,计算dx
    double total_res = 0;
    int effective_num = 0;

    HTVH.setZero();
    HTVr.setZero();

    // 每个点反馈的info信息矩阵因子
    // 由于NDT点数明显多于预测方程，可能导致估计结果向NDT倾斜，
    // 给信息矩阵添加一个乘积因子0.01，让更新部分更加平滑一些。
    const double info_ratio = 0.01;  

    for (int idx = 0; idx < effect_pts.size(); ++idx) {
        if (!effect_pts[idx]) 
            continue;

        total_res += errors[idx] * errors[idx];
        effective_num++;

        HTVH += jacobians[idx].transpose() * jacobians[idx] * info_ratio;    // 18x18
        HTVr += -jacobians[idx].transpose() * errors[idx] * info_ratio;      // 18x1
    }

    LOG(INFO) << "effective: " << effective_num;
}
/**
 * @description: 对去畸变后的scan点云，执行一次配准和观测
 * @return {*}
 */
void LioIEKF::Align() {
    FullCloudPtr scan_undistort_trans(new FullPointCloudType);
    // 将去畸变后的雷达点云转到IMU坐标系下
    pcl::transformPointCloud(*scan_undistort_fullcloud_, *scan_undistort_trans, TIL_.matrix().cast<float>());
    scan_undistort_fullcloud_ = scan_undistort_trans;

    // FullPoint类型的点云转换为PointCloud类型的点云
    scan_undistort_ = ConvertToCloud<FullPointType>(scan_undistort_fullcloud_);

    /// the first scan
    if (flg_first_scan_) {

        ndt_.AddCloud(scan_undistort_);

        // 新增：将第一帧点云放入ivox中
        // ivox_->AddPoints(scan_undistort_->points);
        first_lidar_time_ = measures_.lidar_begin_time_;

        flg_first_scan_ = false;

        // nav_state_prior_ = std::make_shared<PriorNavState>(last_nav_state_, Mat15d::Identity() * 1e4);
        return;
    }

    flg_ESKF_inited_ = (measures_.lidar_begin_time_ - first_lidar_time_) >= options_.INIT_TIME;

    // 后续的scan，使用NDT配合pose进行更新
    LOG(INFO) << "=== frame " << frame_num_;
    // pred_nav_state_ = ieskf_.GetNominalState();

    // 点云降采样
    pcl::VoxelGrid<PointType> voxel;
    voxel.setLeafSize(0.5, 0.5, 0.5);
    voxel.setInputCloud(scan_undistort_);
    voxel.filter(*scan_down_body_); // 体素滤波，降采样

    int cur_pts = scan_down_body_->size(); // 降采样后的去畸变点云数量

    scan_down_world_->resize(cur_pts);
    nearest_points_.resize(cur_pts);

    ndt_.SetSource(scan_down_body_);
    ieskf_.UpdateUsingCustomObserve([this](const SE3 &input_pose, Mat18d &HTVH, Vec18d &HTVr) {
                                        // 使用输入位姿R，t（ESKF的名义状态），得到HTVH
                                        ndt_.ComputeResidualAndJacobians(input_pose, HTVH, HTVr);
                                        // ComputeResidualAndJacobians_P2Plane(input_pose, HTVH, HTVr); // 【新增】
                                    });

    auto current_nav_state = ieskf_.GetNominalState();
    // Optimize();
    // nominal state反馈给ieskf
    // ieskf_.SetX(current_nav_state);

    // 若运动了一定范围，则把点云放入地图中
    SE3 current_pose = ieskf_.GetNominalSE3();
    SE3 delta_pose = last_pose_.inverse() * current_pose;

    // 判断车辆运动是否超过1米或者角度超过10度
    if (delta_pose.translation().norm() > 1.0 || delta_pose.so3().log().norm() > math::deg2rad(10.0)) {
        // 将地图合入NDT中
        // CloudPtr scan_down_world_(new PointCloudType);
        pcl::transformPointCloud(*scan_down_body_, *scan_down_world_, current_pose.matrix());
        ndt_.AddCloud(scan_down_world_); // 更新NDT地图

        /// 更新ivox中的局部地图
        // MapIncremental(); 
        // LOG(INFO) << "[ mapping ]: In num: " << scan_undistort_->points.size() 
        //           << " downsamp " << cur_pts
        //           << " Map grid num: " << ivox_->NumValidGrids();

        // 计算耗时
        // sad::evaluate_and_call([&]() { 
        //                         MapIncremental(); 
        //                     }, 
        //                     "Incremental Mapping", 1);
        last_pose_ = current_pose;
    }

    // 放入UI
    if (ui_) {
        ui_->UpdateScan(scan_undistort_, current_nav_state.GetSE3());  // 转成Lidar Pose传给UI
        ui_->UpdateNavState(current_nav_state);
    }

    frame_num_++;
    return;
}

/**
 * @description: 静止初始化估计重力矢量，陀螺、加计零偏
 * @return {*}
 */
void LioIEKF::TryInitIMU() {
    for (auto imu : measures_.imu_) 
        imu_init_.AddIMU(*imu);

    if (imu_init_.InitSuccess()) {
        // 读取初始零偏，设置ESKF
        sad::IESKFD::Options options;
        // 噪声由初始化器估计
        options.gyro_var_ = sqrt(imu_init_.GetCovGyro()[0]);
        options.acce_var_ = sqrt(imu_init_.GetCovAcce()[0]);
        ieskf_.SetInitialConditions(options, imu_init_.GetInitBg(), imu_init_.GetInitBa(), imu_init_.GetGravity());
        imu_need_init_ = false;

        LOG(INFO) << "IMU初始化成功";
    }
}

/**
 * @description: 点云去畸变，将所有点转换到最后时刻状态上
 * @return {*}
 */
void LioIEKF::Undistort() {
    auto cloud = measures_.lidar_;                  // 原始雷达点云
    auto imu_state = ieskf_.GetNominalState();      // IESKF最后时刻的名义状态
    SE3 T_end = SE3(imu_state.R_, imu_state.p_);    // 位姿R, t

    if (options_.save_motion_undistortion_pcd_) {   // 保存去畸变前的点云
        sad::SaveCloudToFile("./data/ch7/before_undist.pcd", *cloud);
    }

    /// 将所有点转到最后时刻状态上
    std::for_each(std::execution::par_unseq, 
                cloud->points.begin(), cloud->points.end(), 
                [&](auto &pt) {
                    SE3 Ti = T_end;     // 第i时刻的位姿 
                    NavStated match;    // 用于存储匹配到的状态

                    // 根据pt.time查找时间，pt.time是该点打到的时间与雷达开始时间之差，单位为毫秒
                    // 位姿插值：旋转用SLERP球面线性插值，位移用线性插值
                    math::PoseInterp<NavStated>(
                        measures_.lidar_begin_time_ + pt.time * 1e-3, 
                        imu_states_, 
                        [](const NavStated &s) { return s.timestamp_; },
                        [](const NavStated &s) { return s.GetSE3(); }, 
                        Ti,     // 插值得到该点pt的SE3位姿
                        match);

                    Vec3d pi = ToVec3d(pt); // 去畸变前坐标
                    // T_end.inverse() * Ti：最后时刻位姿 到 第i个点的位姿之间的位姿变换，
                    // T_end.inverse() * Ti * TIL_ * pi：相当于将第i个点转换到了最后时刻IESKF状态上
                    // 去畸变后的坐标：从右往左看，先将点云坐标转到IMU坐标系下，然后转换到最后时刻的IESKF状态上，最后转换到雷达坐标系下
                    Vec3d p_compensate = TIL_.inverse() * T_end.inverse() * Ti * TIL_ * pi;
                    pt.x = p_compensate(0);
                    pt.y = p_compensate(1);
                    pt.z = p_compensate(2);
                });
    scan_undistort_ = cloud;    // 去畸变后的点云

    if (options_.save_motion_undistortion_pcd_) {   // 保存去畸变后的点云
        sad::SaveCloudToFile("./data/ch7/after_undist.pcd", *cloud);
    }
}

/**
 * @description: IESKF预测得到的名义状态，并保存每个IMU数据预测的名义状态
 * @return {*}
 */
void LioIEKF::Predict() {
    imu_states_.clear();
    imu_states_.emplace_back(ieskf_.GetNominalState()); // 加入上一次预测的名义状态

    /// 对IMU状态进行预测
     // 遍历同步好的所有IMU数据
    for (auto &imu : measures_.imu_) {
        ieskf_.Predict(*imu);   // 滤波器预测递推得到名义状态 
        imu_states_.emplace_back(ieskf_.GetNominalState()); // 加入名义状态
    }
}

/**
 * @description: 将ROS的点云类型转换为FullCloudPtr类型
 * @param {ConstPtr} &msg
 * @return {*}
 */
void LioIEKF::PCLCallBack(const sensor_msgs::PointCloud2::ConstPtr &msg) { sync_->ProcessCloud(msg); }

void LioIEKF::LivoxPCLCallBack(const livox_ros_driver::CustomMsg::ConstPtr &msg) { sync_->ProcessCloud(msg); }

void LioIEKF::IMUCallBack(IMUPtr msg_in) { sync_->ProcessIMU(msg_in); }

void LioIEKF::Finish() {
    if (ui_) {
        ui_->Quit();
    }
    LOG(INFO) << "finish done";
}

}  // namespace sad