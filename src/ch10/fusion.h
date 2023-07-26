//
// Created by xiang on 22-12-20.
//

#ifndef SLAM_IN_AUTO_DRIVING_FUSION_H
#define SLAM_IN_AUTO_DRIVING_FUSION_H

#include "ch3/eskf.hpp"
#include "ch3/static_imu_init.h"
#include "common/eigen_types.h"
#include "common/gnss.h"
#include "common/imu.h"
#include "common/message_def.h"
#include "common/point_types.h"

#include "ch7/loosely_coupled_lio/cloud_convert.h"
#include "ch7/loosely_coupled_lio/measure_sync.h"

#include "tools/ui/pangolin_window.h"

#include <pcl/registration/ndt.h>
#include <sensor_msgs/PointCloud2.h>

#include "ch7/ndt_3d.h" // 【新增】多分辨率NDT，匹配得分计算

namespace sad {

using KeyType = Eigen::Matrix<int, 3, 1>;  // 体素的索引
/**
 * 第10章显示的高精度融合定位，融合IMU、RTK、激光点云定位功能
 *
 * - NOTE 一些IMU的异常处理没有加在这里，有可能会被IMU带歪。
 */
class Fusion {
public:
    explicit Fusion(const std::string& config_yaml);

    enum class Status {
        WAITING_FOR_RTK,  // 需要等待第一个有效的RTK信号来控制点云定位的搜索范围
        WORKING,          // 正常工作
    };

    /// 初始化，读取参数
    bool Init();

    /// 处理输入的RTK
    void ProcessRTK(GNSSPtr gnss);
    void ProcessIMU(IMUPtr imu);
    void ProcessPointCloud(sensor_msgs::PointCloud2::Ptr cloud);

private:
    /// 加载某个点附近的地图
    void LoadMap(const SE3& pose);

    /// 【新增】加载某个点附近的NDT地图
    void LoadNdtMap(const SE3& pose);

    /// 读取地图的索引文件
    void LoadMapIndex();

    /// 处理同步之后的IMU和雷达数据
    void ProcessMeasurements(const MeasureGroup& meas);

    /// 网格搜索时的结构
    struct GridSearchResult {
        SE3 pose_;              // 搜索的位姿
        SE3 result_pose_;       // 配准后的位姿
        double score_ = 0.0;    // NDT配准得分
    };

    /// 在初始RTK附近搜索车辆位置
    bool SearchRTK();

    /// 对网格搜索的某个点进行配准，得到配准后位姿与配准分值
    void AlignForGrid(GridSearchResult& gr);

    /// 激光定位
    bool LidarLocalization();

    /// 使用IMU初始化
    void TryInitIMU();

    /// 利用IMU预测状态信息
    /// 这段时间的预测数据会放入imu_states_里
    void Predict();

    /// 对measures_中的点云去畸变
    void Undistort();

    /// 执行一次配准和观测
    void Align();

    /// 标志位
    Status status_ = Status::WAITING_FOR_RTK; // 默认为等待初始RTK信号

    /// 数据
    std::string config_yaml_;                           // 配置文件路径
    Vec3d map_origin_ = Vec3d::Zero();                  // 地图原点
    std::string pts_map_path_;                          // 地图数据目录
    std::set<Vec2i, less_vec<2>> pts_map_data_index_;   // 哪些格子存在点云地图数据
    std::map<Vec2i, CloudPtr, less_vec<2>> map_data_;   // 第9章建立的地图数据

    //【新增】NDT地图数据目录
    std::string ndt_map_path_10_, ndt_map_path_5_, ndt_map_path_4_, 
                ndt_map_path_3_, ndt_map_path_1_;                        

    //【新增】哪些格子存在NDT地图数据
    std::set<Vec2i, less_vec<2>> ndt_map_data_index_10_, ndt_map_data_index_5_, ndt_map_data_index_4_, 
                                 ndt_map_data_index_3_, ndt_map_data_index_1_; 

    using NdtmapVoxel = std::unordered_map<KeyType, Ndt3d::VoxelData, hash_vec<3>>; //【新增】

    NdtmapVoxel ndt_map_grids_10_, ndt_map_grids_5_, ndt_map_grids_4_
                ndt_map_grids_3_, ndt_map_grids_1_;    // 【新增】

    // 以网格ID为索引的地图数据
    std::map<Vec2i, NdtmapVoxel, less_vec<2>>   ndt_map_data_10_, ndt_map_data_5_, ndt_map_data_2_,
                                                ndt_map_data_3_, ndt_map_data_1_; 

    std::shared_ptr<MessageSync> sync_ = nullptr;  // 消息同步器
    StaticIMUInit imu_init_;                       // IMU静止初始化

    /// 用到第三章中的ESKF滤波器
    ESKFD eskf_;
    std::vector<NavStated> imu_states_;  // ESKF预测期间的状态

    // 去畸变后的点云 scan after undistortion 
    FullCloudPtr scan_undistort_{new FullPointCloudType()};  
    CloudPtr current_scan_ = nullptr; // 当前帧扫描点云

    SE3 TIL_;   // 雷达到IMU的外参
    MeasureGroup measures_;  // sync IMU and lidar scan
    GNSSPtr last_gnss_ = nullptr;   // 上一个RTK数据

    bool init_has_failed_ = false;  // 初始化是否失败过
    SE3 last_searched_pos_;         // 上次搜索的GNSS位置

    /// 激光定位
    bool imu_need_init_ = true;     // 是否需要估计IMU初始零偏
    CloudPtr ref_cloud_ = nullptr;  // NDT用于参考的点云

    // PCL NDT
    pcl::NormalDistributionsTransform<PointType, PointType> ndt_pcl_;

    // TODO: 换成自己实现的3D NDT，需要实现多分辨率NDT，以及匹配得分的计算
    Ndt3d ndt_;

    /// 参数
    double rtk_search_min_score_ = 0.16; //3.5; // 最小得分需修改为合适的，适配自己实现的NDT匹配得分

    // ui
    std::shared_ptr<ui::PangolinWindow> ui_ = nullptr;

    bool use_pcl_ndt_ = false;   // 是否使用PCL NDT库
};

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_FUSION_H
