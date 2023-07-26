//
// Created by xiang on 2022/7/18.
//

#include "ch7/direct_ndt_lo.h"
#include "common/math_utils.h"
#include "tools/pcl_map_viewer.h"

#include <pcl/common/transforms.h>

namespace sad {

/**
 * @description: 添加点云，估计位姿
 * @param {CloudPtr} scan
 * @param {SE3&} pose
 * @return {*}
 */
void DirectNDTLO::AddCloud(CloudPtr scan, SE3& pose) {
    // 判断局部地图是否为空，还没有初始化
    if (local_map_ == nullptr) {
        // 第一个帧，直接加入local map
        local_map_.reset(new PointCloudType);
        *local_map_ += *scan;   // 往局部地图中添加点云，+是pcl中的运算符重载，表示将两个点云拼接起来
        pose = SE3();           // 位姿初始化
        last_kf_pose_ = pose;   // 上一次关键帧位姿

        if (options_.use_pcl_ndt_) 
            ndt_pcl_.setInputTarget(local_map_);    // PCL的NDT
        else 
            ndt_.SetTarget(local_map_);            // SAD中Ndt3d

        return;
    }

    // 计算当前帧scan点云相对于local map的位姿
    pose = AlignWithLocalMap(scan);

    // 将当前帧点云转换到世界坐标系下
    CloudPtr scan_world(new PointCloudType);
    pcl::transformPointCloud(*scan, *scan_world, pose.matrix().cast<float>());

    // 每隔一定距离或者一定角度取一个关键帧，再将最近的若干个关键帧拼成一个局部地图，作为NDT算法的配准目标。
    if (IsKeyframe(pose)) {
        // 仅对关键帧进行处理
        last_kf_pose_ = pose; // 记录当前关键帧的位姿

        // 重建local map
        // 将当前帧世界系下的点云加入到局部地图中
        scans_in_local_map_.emplace_back(scan_world);
        // 判断局部地图中关键帧的数量是否超过阈值，超过则删除最早的关键帧
        if (scans_in_local_map_.size() > options_.num_kfs_in_local_map_) 
            scans_in_local_map_.pop_front();

        local_map_.reset(new PointCloudType); // 重置局部地图
        for (auto& scan : scans_in_local_map_) 
            *local_map_ += *scan; // 将局部地图中的所有关键帧点云拼接起来

        // 将最近若干个关键帧拼成的局部地图作为NDT算法的配准目标
        if (options_.use_pcl_ndt_) 
            ndt_pcl_.setInputTarget(local_map_);
        else 
            ndt_.SetTarget(local_map_);
    }

    if (viewer_ != nullptr)  
        viewer_->SetPoseAndCloud(pose, scan_world); // 若viewer非空，则显示当前帧点云和位姿
}

/**
 * @description: 判定是否为关键帧
 * @param {SE3&} current_pose
 * @return {*}
 */
bool DirectNDTLO::IsKeyframe(const SE3& current_pose) {
    SE3 delta = last_kf_pose_.inverse() * current_pose; // 计算上一个关键帧与当前帧之间的相对位姿
    // 只要与上一帧相对运动超过一定距离或角度，就记为关键帧
    return delta.translation().norm() > options_.kf_distance_ ||                // 距离是否超过阈值
           delta.so3().log().norm() > options_.kf_angle_deg_ * math::kDEG2RAD;  // 角度是否超过阈值
}

/**
 * @description: 与local map进行配准
 * @param {CloudPtr} scan
 * @return {*}
 */
SE3 DirectNDTLO::AlignWithLocalMap(CloudPtr scan) {
    // 将当前帧点云设置为NDT算法的输入，待配准
    if (options_.use_pcl_ndt_) 
        ndt_pcl_.setInputSource(scan);
    else 
        ndt_.SetSource(scan);

    CloudPtr output(new PointCloudType());  // 配准后的点云

    SE3 guess;
    bool align_success = true;
    // 判断估计的位姿数目是否小于2个
    if (estimated_poses_.size() < 2) {
        if (options_.use_pcl_ndt_) {
            ndt_pcl_.align(*output, guess.matrix().cast<float>());
            guess = Mat4ToSE3(ndt_pcl_.getFinalTransformation().cast<double>().eval()); 
        } else 
            align_success = ndt_.AlignNdt(guess);   // SAD中Ndt3d
    } else { // 多于2个时，
        // 从最近两个pose来推断
        SE3 T1 = estimated_poses_[estimated_poses_.size() - 1]; // 倒数第一个最近的
        SE3 T2 = estimated_poses_[estimated_poses_.size() - 2]; // 倒数第二个次近的
        guess = T1 * (T2.inverse() * T1);   // 计算

        if (options_.use_pcl_ndt_) {
            ndt_pcl_.align(*output, guess.matrix().cast<float>());
            guess = Mat4ToSE3(ndt_pcl_.getFinalTransformation().cast<double>().eval());
        } else 
            align_success = ndt_.AlignNdt(guess);
    }

    LOG(INFO) << "pose: " << guess.translation().transpose() << ", "
              << guess.so3().unit_quaternion().coeffs().transpose();

    if (options_.use_pcl_ndt_) 
        LOG(INFO) << "trans prob: " << ndt_pcl_.getTransformationProbability();

    estimated_poses_.emplace_back(guess);
    return guess;
}

void DirectNDTLO::SaveMap(const std::string& map_path) {
    if (viewer_) 
        viewer_->SaveMap(map_path);
}

}  // namespace sad