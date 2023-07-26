//
// Created by xiang on 22-12-6.
//

#include "frontend.h"
#include "ch8/lio-iekf/lio_iekf.h"
#include "common/io_utils.h"

#include <yaml-cpp/yaml.h>

namespace sad {

Frontend::Frontend(const std::string& config_yaml) { config_yaml_ = config_yaml; }

bool Frontend::Init() {
    LOG(INFO) << "load yaml from " << config_yaml_;
    auto yaml = YAML::LoadFile(config_yaml_);
    try {
        auto n = yaml["bag_path"];
        LOG(INFO) << Dump(n);
        bag_path_ = yaml["bag_path"].as<std::string>();
        lio_yaml_ = yaml["lio_yaml"].as<std::string>();
    } catch (...) {
        LOG(ERROR) << "failed to parse yaml";
        return false;
    }

    system("rm -rf ./data/ch9/*.pcd");
    system("rm -rf ./data/ch9/keyframes.txt");

    LioIEKF::Options options;
    options.with_ui_ = false;  // 跑建图不需要打开前端UI
    lio_ = std::make_shared<LioIEKF>(options);
    lio_->Init(lio_yaml_);
    return true;
}

/**
 * @description: 运行前端，保存关键帧点云到本地
 * @return {*}
 */
void Frontend::Run() {
    sad::RosbagIO rosbag_io(bag_path_, DatasetType::NCLT);

    // 先提取RTK pose，注意NCLT只有平移部分（单天线RTK方案，仅有平移信息不含位移信息）
    rosbag_io
        .AddAutoRTKHandle([this](GNSSPtr gnss) {
            // 将ROS包中的RTK数据提取出来，放在RTK消息队列中，按照采集时间进行排序
            gnss_.emplace(gnss->unix_time_, gnss);
            return true;
        })
        .Go();

    // 此时队列中只有一个RTK回调处理函数，将其清空，不再需要处理RTK
    rosbag_io.CleanProcessFunc();  

    // 从RTK数据中找到第一个有效的RTK数据，将其作为地图原点，
    // 将其它RTK减去原点后，作为RTK的位置观测值
    RemoveMapOrigin();

    // 再用IMU和激光点云数据运行LIO
    rosbag_io
        .AddAutoPointCloudHandle([&](sensor_msgs::PointCloud2::Ptr cloud) -> bool {
            // 将ROS的点云类型转换为FullCloudPtr类型；
            // 添加新的点云缓存和时间戳缓存；
            // 调用Sync()函数执行lidar与IMU数据的同步
            lio_->PCLCallBack(cloud);
            
            // 按照距离和角度阈值抽取关键帧，并保存关键帧点云pcd文件到本地。
            ExtractKeyFrame(lio_->GetCurrentState());
            return true;
        })
        .AddImuHandle([&](IMUPtr imu) {
            lio_->IMUCallBack(imu);
            return true;
        })
        .Go();
    lio_->Finish();

    // 保存运行结果
    SaveKeyframes();

    LOG(INFO) << "done.";
}

/**
 * @description: 按照距离和角度阈值抽取关键帧
 * @param {NavStated&} state
 * @return {*}
 */
void Frontend::ExtractKeyFrame(const sad::NavStated& state) {
    // 判断上一个关键帧是否为空
    if (last_kf_ == nullptr) {
        // 若为空，判断当前scan扫描得到的点云是否为空
        if (!lio_->GetCurrentScan()) 
            return; // 若为空，说明LIO还没完成初始化，直接返回
        // 构建第一个关键帧
        auto kf = std::make_shared<Keyframe>(state.timestamp_,          // 时间戳
                                             kf_id_++,                  // 关键帧ID
                                             state.GetSE3(),            // SE3位姿
                                             lio_->GetCurrentScan());   // 当前帧激光点云
        // 按照时间戳查询RTK位姿列表，调用math中的插值函数得到关键帧点云的位姿
        FindGPSPose(kf);
        kf->SaveAndUnloadScan("./data/ch9/");   // 保存关键帧点云到指定路径，以关键帧ID命名pcd文件
        keyframes_.emplace(kf->id_, kf);        // 保存关键帧点云到指定路径，以关键帧ID命名pcd文件  
        last_kf_ = kf;                          // 更新上一个关键帧
    } else {    // 若上一个关键帧非空
        // 计算当前IESKF状态state与上一次kf关键帧之间的相对运动SE3位姿
        SE3 delta = last_kf_->lidar_pose_.inverse() * state.GetSE3();
        // 判断相对运动是否超过距离或者角度阈值
        if (delta.translation().norm() > kf_dis_th_ || 
            delta.so3().log().norm() > kf_ang_th_deg_ * math::kDEG2RAD) {
            // 距离或者角度其中一个超过阈值，构建新的关键帧
            auto kf = std::make_shared<Keyframe>(state.timestamp_,          // 时间戳
                                                 kf_id_++,                  // 关键帧ID
                                                 state.GetSE3(),            // SE3位姿
                                                 lio_->GetCurrentScan());   // 当前帧激光点云
            // 插值RTK位姿队列获取关键帧时间戳处的SE3位姿（旋转用SLERP，平移用线性插值）
            FindGPSPose(kf); // 这样就为每个关键帧点云找到了对应的LIO、RTK位姿以及扫描到的点云
            keyframes_.emplace(kf->id_, kf);        // 将关键帧加入关键帧列表，按照关键帧ID排序
            kf->SaveAndUnloadScan("./data/ch9/");   // 保存关键帧点云到指定路径，以关键帧ID命名pcd文件  
            LOG(INFO) << "生成关键帧" << kf->id_;
            last_kf_ = kf;                          // 更新上一个关键帧
        }
    }
}

void Frontend::FindGPSPose(std::shared_ptr<Keyframe> kf) {
    SE3 pose;
    GNSSPtr match;
    // 插值RTK位姿队列获取关键帧时间戳处的SE3位姿，旋转用SLERP，平移用线性插值
    if (math::PoseInterp<GNSSPtr>(  kf->timestamp_, 
                                    gnss_, 
                                    [](const GNSSPtr& gnss) -> SE3 {    // 用户提供获取位姿的lambda方法
                                        return gnss->utm_pose_;         // 对于RTK位姿读数来说，返回utm_pose即可
                                    }, 
                                    pose, 
                                    match)) {
        kf->rtk_pose_ = pose;
        kf->rtk_valid_ = true;
    } else {
        // 插值失败，说明RTK数据不足，或者，关键帧查询时间戳晚于最后一个RTK数据，RTK数据无效
        kf->rtk_valid_ = false;
    }
}

void Frontend::SaveKeyframes() {
    std::ofstream fout("./data/ch9/keyframes.txt");
    for (auto& kfp : keyframes_) {
        kfp.second->Save(fout);
    }
    fout.close();
}

/**
 * @description: 从RTK数据中找到第一个有效的RTK数据，将其作为地图原点，将其它RTK减去原点后，作为RTK的位置观测值
 * @return {*}
 */
void Frontend::RemoveMapOrigin() {
    // 判断是否有RTK数据
    if (gnss_.empty()) 
        return;

    bool origin_set = false;
    // 遍历所有RTK数据
    for (auto& p : gnss_) {
        // 判断RTK数据是否有效
        if (p.second->status_ == GpsStatusType::GNSS_FIXED_SOLUTION) {
            map_origin_ = p.second->utm_pose_.translation();    // 获取第一个有效RTK数据的位置
            origin_set = true;  // 标记已经找到第一个有效RTK数据，可用于设置地图原点

            // 地图原点设置在第一个有效RTK数据的位置
            LOG(INFO) << "map origin is set to " << map_origin_.transpose();

            auto yaml = YAML::LoadFile(config_yaml_);
            std::vector<double> ori{map_origin_[0], map_origin_[1], map_origin_[2]};
            yaml["origin"] = ori;
            std::ofstream fout(config_yaml_);
            fout << yaml;   // 将地图原点写入yaml文件

            break;  // 找到第一个有效RTK数据后，就退出循环
        }
    }

    if (origin_set) {
        LOG(INFO) << "removing origin from rtk";
        for (auto& p : gnss_) {
            p.second->utm_pose_.translation() -= map_origin_; // 将所有RTK数据减去地图原点
        }
    }
}

}  // namespace sad