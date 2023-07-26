#include <yaml-cpp/yaml.h>
#include <execution>

#include "common/lidar_utils.h"
#include "common/timer/timer.h"
#include "loosely_lio.h"

namespace sad {

/**
 * @description: 构造函数
 * @param {Options} options
 * @return {*}
 */
LooselyLIO::LooselyLIO(Options options) : options_(options) {
    StaticIMUInit::Options imu_init_options;
    imu_init_options.use_speed_for_static_checking_ = false;  // 本节数据不需要轮速计
    imu_init_ = StaticIMUInit(imu_init_options);
}

/**
 * @description: 初始化
 * @param {string} &config_yaml
 * @return {*}
 */
bool LooselyLIO::Init(const std::string &config_yaml) {
    /// 初始化自身的参数
    if (!LoadFromYAML(config_yaml)) 
        return false;

    /// 初始化NDT LO的参数
    sad::IncrementalNDTLO::Options indt_options;
    indt_options.display_realtime_cloud_ = false;  // 这个程序自己有UI，不用PCL中的
    inc_ndt_lo_ = std::make_shared<sad::IncrementalNDTLO>(indt_options);

    /// 初始化UI
    if (options_.with_ui_) {
        ui_ = std::make_shared<ui::PangolinWindow>();
        ui_->Init();
    }

    return true;
}

/**
 * @description: 从yaml配置文件中提取参数
 * @param {string} &yaml_file
 * @return {*}
 */
bool LooselyLIO::LoadFromYAML(const std::string &yaml_file) {
    // get params from yaml
    sync_ = std::make_shared<MessageSync>([this](const MeasureGroup &m) { ProcessMeasurements(m); });
    sync_->Init(yaml_file);

    /// 自身参数主要是雷达与IMU外参
    auto yaml = YAML::LoadFile(yaml_file);
    std::vector<double> ext_t = yaml["mapping"]["extrinsic_T"].as<std::vector<double>>();
    std::vector<double> ext_r = yaml["mapping"]["extrinsic_R"].as<std::vector<double>>();

    Vec3d lidar_T_wrt_IMU = math::VecFromArray(ext_t);
    Mat3d lidar_R_wrt_IMU = math::MatFromArray(ext_r);
    TIL_ = SE3(lidar_R_wrt_IMU, lidar_T_wrt_IMU);   // TIL_是雷达相对于IMU的外参
    return true;
}

void LooselyLIO::ProcessMeasurements(const MeasureGroup &meas) {
    LOG(INFO) << "call meas, imu: " << meas.imu_.size() << ", lidar pts: " << meas.lidar_->size();
    measures_ = meas;

    if (imu_need_init_) {
        // 初始化IMU系统
        TryInitIMU();
        return;
    }

    // 利用IMU数据进行状态预测
    Predict();

    // 对点云去畸变
    Undistort();

    // 配准
    Align();
}

/**
 * @description: ESKF名义状态预测
 * @return {*}
 */
void LooselyLIO::Predict() {
    imu_states_.clear(); // 清空IMU状态
    imu_states_.emplace_back(eskf_.GetNominalState());

    /// 对IMU状态进行预测
    for (auto &imu : measures_.imu_) {
        eskf_.Predict(*imu);
        imu_states_.emplace_back(eskf_.GetNominalState());
    }
}

void LooselyLIO::TryInitIMU() {
    for (auto imu : measures_.imu_) {
        imu_init_.AddIMU(*imu);
    }

    if (imu_init_.InitSuccess()) {
        // 读取初始零偏，设置ESKF
        sad::ESKFD::Options options;
        // 噪声由初始化器估计
        options.gyro_var_ = sqrt(imu_init_.GetCovGyro()[0]);
        options.acce_var_ = sqrt(imu_init_.GetCovAcce()[0]);
        eskf_.SetInitialConditions(options, imu_init_.GetInitBg(), imu_init_.GetInitBa(), imu_init_.GetGravity());
        imu_need_init_ = false;

        LOG(INFO) << "IMU初始化成功";
    }
}

/**
 * @description: 去畸变
 *  
 *  单次扫描需要ts时间，雷达在不同时刻接收到的所有激光点，全部转换到最后时刻ts上
 * 
 * @return {*}
 */
void LooselyLIO::Undistort() {
    auto cloud = measures_.lidar_;
    auto imu_state = eskf_.GetNominalState();  // 最后时刻的状态
    SE3 T_end = SE3(imu_state.R_, imu_state.p_);

    if (options_.save_motion_undistortion_pcd_) 
        pcl::io::savePCDFileBinary("./data/ch7/before_undist.pcd", *cloud);

    /// 单次扫描需要ts时间，雷达在不同时刻接收到的所有激光点，全部转换到最后时刻ts上
    std::for_each(  std::execution::par_unseq, 
                    cloud->points.begin(), cloud->points.end(), 
                    [&](auto &pt) { 
                        SE3 Ti = T_end;
                        NavStated match;

                        // 某个点的pt，雷达返回的时间是pt.time，
                        // 激光雷达单次扫描时间为ts，那么pt.time位于（0, ts）之间
                        // 根据pt.time查找时间，目的是通过ESKF名义状态递推预测得到的每个IMU时刻的位姿来插值得到pt.time时刻的位姿
                        // 旋转采用SLERP球面线性插值，平移采用线性插值
                        math::PoseInterp<NavStated>(measures_.lidar_begin_time_ + pt.time * 1e-3, 
                                                    imu_states_, 
                                                    [](const NavStated &s) { return s.timestamp_; },    // 时间戳
                                                    [](const NavStated &s) { return s.GetSE3(); },      // SE3位姿
                                                    Ti, // 插值得到的pt.time时刻的位姿 
                                                    match);

                        Vec3d pi = ToVec3d(pt);
                        // TIL_ * pi： 雷达系转换到IMU系
                        // T_end.inverse() * Ti： pt.time时刻与最后时刻之间的位姿差
                        // 连起来就是，pt.time时刻雷达系下的点转换到IMU系下，再转换到最后时刻的IMU系下，最后转换到雷达系下
                        Vec3d p_compensate = TIL_.inverse() * T_end.inverse() * Ti * TIL_ * pi;

                        pt.x = p_compensate(0);
                        pt.y = p_compensate(1);
                        pt.z = p_compensate(2);
                    });
    scan_undistort_ = cloud; // 运动补偿后的激光点，已经都在最后时刻的雷达系下

    if (options_.save_motion_undistortion_pcd_) 
        pcl::io::savePCDFileBinary("./data/ch7/after_undist.pcd", *cloud);
}

/**
 * @description: 点云配准
 *          
 *      将运动补偿去畸变后的点云进行配准，激光里程计LO得到的位姿再给ESKF。
 *          
 * @return {*}
 */
void LooselyLIO::Align() {
    FullCloudPtr scan_undistort_trans(new FullPointCloudType);
    pcl::transformPointCloud(*scan_undistort_, *scan_undistort_trans, TIL_.matrix());
    scan_undistort_ = scan_undistort_trans;

    auto current_scan = ConvertToCloud<FullPointType>(scan_undistort_);

    // voxel 之
    pcl::VoxelGrid<PointType> voxel;
    voxel.setLeafSize(0.5, 0.5, 0.5);
    voxel.setInputCloud(current_scan);

    CloudPtr current_scan_filter(new PointCloudType);
    voxel.filter(*current_scan_filter);

    /// 处理首帧雷达数据
    if (flg_first_scan_) {
        SE3 pose;
        inc_ndt_lo_->AddCloud(current_scan_filter, pose);
        flg_first_scan_ = false;
        return;
    }

    /// 从ESKF中获取预测pose，放入LO，获取LO位姿，最后合入EKF
    SE3 pose_predict = eskf_.GetNominalSE3();

    // 增量式NDT激光里程计估计位姿
    inc_ndt_lo_->AddCloud(current_scan_filter, pose_predict, true);
    pose_of_lo_ = pose_predict;

    // LO激光里程计得到的SE3位姿，放入ESKF进行观测更新
    eskf_.ObserveSE3(pose_of_lo_, 1e-2, 1e-2); // 平移和旋转噪声都设置为1e-2

    if (options_.with_ui_) {
        // 放入UI
        ui_->UpdateScan(current_scan, eskf_.GetNominalSE3());  // 转成Lidar Pose传给UI
        ui_->UpdateNavState(eskf_.GetNominalState());
    }
    frame_num_++;
}

void LooselyLIO::PCLCallBack(const sensor_msgs::PointCloud2::ConstPtr &msg) { sync_->ProcessCloud(msg); }

void LooselyLIO::LivoxPCLCallBack(const livox_ros_driver::CustomMsg::ConstPtr &msg) { sync_->ProcessCloud(msg); }

void LooselyLIO::IMUCallBack(IMUPtr msg_in) { sync_->ProcessIMU(msg_in); }

void LooselyLIO::Finish() {
    if (options_.with_ui_) {
        while (ui_->ShouldQuit() == false) 
            usleep(1e5);
        ui_->Quit();
    }
    LOG(INFO) << "finish done";
}

}  // namespace sad