//
// Created by xiang on 2022/7/25.
//

#include "ch7/loam-like/loam_like_odom.h"
#include "ch7/loam-like/feature_extraction.h"
#include "common/lidar_utils.h"
#include "common/math_utils.h"

#include <execution>



namespace sad {

/**
 * @description: 构造函数
 * @param {Options} options
 * @return {*}
 */
LoamLikeOdom::LoamLikeOdom(LoamLikeOdom::Options options) : options_(options), feature_extraction_(new FeatureExtraction), global_map_(new PointCloudType()) {

    if (options_.display_realtime_cloud_) 
        viewer_ = std::make_shared<PCLMapViewer>(0.1);

    kdtree_edge_.SetEnableANN();
    kdtree_surf_.SetEnableANN();
    kdtree_ground_.SetEnableANN(); // 【新增】
}

/**
 * @description: 处理当前帧点云
 * @param {FullCloudPtr} cloud
 * @return {*}
 */
void LoamLikeOdom::ProcessPointCloud(FullCloudPtr cloud) {
    LOG(INFO) << "processing frame " << cnt_frame_++;
    /// 【step 1】提特征
    CloudPtr current_edge(new PointCloudType),      // 当前帧的边缘点云
             current_surf(new PointCloudType),      // 当前帧的平面点云
             current_ground(new PointCloudType);  // 地面点云【新增】 LegoLOAM方法
             
    // 提取角点和平面
    feature_extraction_->Extract(cloud, current_edge, current_surf, current_ground);

    if (current_edge->size() < options_.min_edge_pts_ || 
        current_surf->size() < options_.min_surf_pts_ ||
        current_ground->size() < options_.min_ground_pts_) {
        LOG(ERROR) << "not enough edge/surf/ground pts: " << current_edge->size() << "," << current_surf->size()<< "," << current_ground->size();
        return;
    }

    LOG(INFO) << "edge: " << current_edge->size() << ", surf: " << current_surf->size() << ", ground: " << current_ground->size();

    // 首帧特殊处理
    if (local_map_edge_ == nullptr || local_map_surf_ == nullptr || local_map_ground_ == nullptr) {
        
        local_map_edge_ = current_edge;
        local_map_surf_ = current_surf;
        local_map_ground_ = current_ground; // 【新增】

        // 分别对角点和平面点云构建kd树
        kdtree_edge_.BuildTree(local_map_edge_);    
        kdtree_surf_.BuildTree(local_map_surf_);
        kdtree_ground_.BuildTree(local_map_ground_); // 【新增】

        edges_.emplace_back(current_edge);
        surfs_.emplace_back(current_surf);
        ground_.emplace_back(current_ground); // 【新增】
        return;
    }

    /// 与局部地图配准
    SE3 pose = AlignWithLocalMap(current_edge, current_surf, current_ground);


    CloudPtr scan_world(new PointCloudType);
    pcl::transformPointCloud(*ConvertToCloud<FullPointType>(cloud), *scan_world, pose.matrix());

    CloudPtr edge_world(new PointCloudType), 
             surf_world(new PointCloudType),
             ground_world(new PointCloudType); // 【新增】

    pcl::transformPointCloud(*current_edge, *edge_world, pose.matrix());
    pcl::transformPointCloud(*current_surf, *surf_world, pose.matrix());
    pcl::transformPointCloud(*current_ground, *ground_world, pose.matrix()); // 【新增】

    if (IsKeyframe(pose)) {
        LOG(INFO) << "inserting keyframe";
        last_kf_pose_ = pose;
        last_kf_id_ = cnt_frame_;

        // 重建local map
        edges_.emplace_back(edge_world);
        surfs_.emplace_back(surf_world);
        ground_.emplace_back(ground_world); // 【新增】

        if (edges_.size() > options_.num_kfs_in_local_map_) {
            edges_.pop_front();
        }
        if (surfs_.size() > options_.num_kfs_in_local_map_) {
            surfs_.pop_front();
        }
        if (ground_.size() > options_.num_kfs_in_local_map_) { // 【新增】
            ground_.pop_front();
        }

        local_map_surf_.reset(new PointCloudType);
        local_map_edge_.reset(new PointCloudType);
        local_map_ground_.reset(new PointCloudType); // 【新增】

        for (auto& s : edges_) 
            *local_map_edge_ += *s;
        
        for (auto& s : surfs_) 
            *local_map_surf_ += *s;

        for (auto& s : ground_) // 【新增】
            *local_map_ground_ += *s;

        local_map_surf_ = VoxelCloud(local_map_surf_, 1.0);
        local_map_edge_ = VoxelCloud(local_map_edge_, 1.0);
        local_map_ground_ = VoxelCloud(local_map_ground_, 1.0); // 【新增】

        LOG(INFO) << "insert keyframe, surf pts: " << local_map_surf_->size() 
                  << ", edge pts: " << local_map_edge_->size()
                  << ", ground pts: " << local_map_ground_->size();

        kdtree_surf_.BuildTree(local_map_surf_);
        kdtree_edge_.BuildTree(local_map_edge_);
        kdtree_ground_.BuildTree(local_map_ground_); // 【新增】

        *global_map_ += *scan_world;
    }

    LOG(INFO) << "current pose: " << pose.translation().transpose() << ", "
              << pose.so3().unit_quaternion().coeffs().transpose();

    if (viewer_ != nullptr) 
        viewer_->SetPoseAndCloud(pose, scan_world);
}

bool LoamLikeOdom::IsKeyframe(const SE3& current_pose) {
    if ((cnt_frame_ - last_kf_id_) > 30) 
        return true;
    
    // 只要与上一帧相对运动超过一定距离或角度，就记关键帧
    SE3 delta = last_kf_pose_.inverse() * current_pose;
    return delta.translation().norm() > options_.kf_distance_ ||
           delta.so3().log().norm() > options_.kf_angle_deg_ * math::kDEG2RAD;
}

/**
 * @description: 与局部地图进行配准
 * @param {CloudPtr} edge
 * @param {CloudPtr} surf
 * @param {CloudPtr} ground 【新增】
 * @return {*}
 */
SE3 LoamLikeOdom::AlignWithLocalMap(CloudPtr edge, CloudPtr surf, CloudPtr ground) {
    // 这部分的ICP需要自己写
    SE3 pose;
    if (estimated_poses_.size() >= 2) {
        // 从最近两个pose来推断
        SE3 T1 = estimated_poses_[estimated_poses_.size() - 1];
        SE3 T2 = estimated_poses_[estimated_poses_.size() - 2];
        pose = T1 * (T2.inverse() * T1);
    }

    int edge_size = edge->size();
    int surf_size = surf->size();
    int ground_size = ground->size(); // 【新增】

    // 我们来写一些并发代码
    for (int iter = 0; iter < options_.max_iteration_; ++iter) {
        std::vector<bool> effect_surf(surf_size, false);
        std::vector<Eigen::Matrix<double, 1, 6>> jacob_surf(surf_size);  // 点面的残差是1维的
        std::vector<double> errors_surf(surf_size);

        std::vector<bool> effect_ground(ground_size, false);                  // 【新增】
        std::vector<Eigen::Matrix<double, 1, 6>> jacob_ground(ground_size);  // 点面的残差是1维的
        std::vector<double> errors_ground(ground_size);

        std::vector<bool> effect_edge(edge_size, false);
        std::vector<Eigen::Matrix<double, 3, 6>> jacob_edge(edge_size);  // 点线的残差是3维的
        std::vector<Vec3d> errors_edge(edge_size);

        std::vector<int> index_surf(surf_size);
        std::iota(index_surf.begin(), index_surf.end(), 0);  // 填入 0,1,2,3,4...
        std::vector<int> index_edge(edge_size);
        std::iota(index_edge.begin(), index_edge.end(), 0);  // 填入 0,1,2,3,4...
        std::vector<int> index_ground(ground_size);
        std::iota(index_ground.begin(), index_ground.end(), 0);  // 填入 0,1,2,3,4...

        // gauss-newton 迭代
        // 最近邻，角点部分
        if (options_.use_edge_points_) {
            std::for_each(std::execution::par_unseq, index_edge.begin(), index_edge.end(), [&](int idx) {
                Vec3d q = ToVec3d(edge->points[idx]);
                Vec3d qs = pose * q;

                // 检查最近邻
                std::vector<int> nn_indices;
                // 搜索5个最近邻
                kdtree_edge_.GetClosestPoint(ToPointType(qs), nn_indices, 5);
                effect_edge[idx] = false;

                if (nn_indices.size() >= 3) {
                    std::vector<Vec3d> nn_eigen;
                    for (auto& n : nn_indices) 
                        nn_eigen.emplace_back(ToVec3d(local_map_edge_->points[n]));

                    // point to point residual 点到点残差
                    Vec3d d, p0;
                    if (!math::FitLine(nn_eigen, p0, d, options_.max_line_distance_)) 
                        return;

                    Vec3d err = SO3::hat(d) * (qs - p0);
                    if (err.norm() > options_.max_line_distance_) 
                        return;

                    effect_edge[idx] = true;

                    // 雅可比矩阵，对应公式（7.11）
                    Eigen::Matrix<double, 3, 6> J;
                    J.block<3, 3>(0, 0) = -SO3::hat(d) * pose.so3().matrix() * SO3::hat(q);
                    J.block<3, 3>(0, 3) = SO3::hat(d);

                    jacob_edge[idx] = J;
                    errors_edge[idx] = err;
                }
            });
        }

        /// 最近邻，平面点部分，点面ICP
        if (options_.use_surf_points_) {
            std::for_each(std::execution::par_unseq, index_surf.begin(), index_surf.end(), [&](int idx) {
                Vec3d q = ToVec3d(surf->points[idx]);
                Vec3d qs = pose * q;

                // 检查最近邻
                std::vector<int> nn_indices;

                kdtree_surf_.GetClosestPoint(ToPointType(qs), nn_indices, 5);
                effect_surf[idx] = false;

                if (nn_indices.size() == 5) {
                    std::vector<Vec3d> nn_eigen;
                    for (auto& n : nn_indices) 
                        nn_eigen.emplace_back(ToVec3d(local_map_surf_->points[n]));

                    // 点面残差
                    Vec4d n;
                    if (!math::FitPlane(nn_eigen, n)) 
                        return;

                    double dis = n.head<3>().dot(qs) + n[3];
                    if (fabs(dis) > options_.max_plane_distance_) 
                        return;

                    effect_surf[idx] = true;

                    // build residual
                    Eigen::Matrix<double, 1, 6> J;
                    J.block<1, 3>(0, 0) = -n.head<3>().transpose() * pose.so3().matrix() * SO3::hat(q);
                    J.block<1, 3>(0, 3) = n.head<3>().transpose();

                    jacob_surf[idx] = J;
                    errors_surf[idx] = dis;
                }
            });
        }

        /// 最近邻，地面点部分，点面ICP
        if (options_.use_ground_points_) {
            std::for_each(std::execution::par_unseq, index_ground.begin(), index_ground.end(), [&](int idx) {
                Vec3d q = ToVec3d(ground->points[idx]);
                Vec3d qs = pose * q;

                // 检查最近邻
                std::vector<int> nn_indices;

                kdtree_ground_.GetClosestPoint(ToPointType(qs), nn_indices, 5);
                effect_ground[idx] = false;

                if (nn_indices.size() == 5) {
                    std::vector<Vec3d> nn_eigen;
                    for (auto& n : nn_indices) 
                        nn_eigen.emplace_back(ToVec3d(local_map_ground_->points[n]));

                    // 点面残差
                    Vec4d n;
                    if (!math::FitPlane(nn_eigen, n)) 
                        return;

                    double dis = n.head<3>().dot(qs) + n[3];
                    if (fabs(dis) > options_.max_plane_distance_) 
                        return;

                    effect_ground[idx] = true;

                    // build residual
                    Eigen::Matrix<double, 1, 6> J;
                    J.block<1, 3>(0, 0) = -n.head<3>().transpose() * pose.so3().matrix() * SO3::hat(q);
                    J.block<1, 3>(0, 3) = n.head<3>().transpose();

                    jacob_ground[idx] = J;
                    errors_ground[idx] = dis;
                }
            });
        }

        // 累加Hessian和error,计算dx
        // 原则上可以用reduce并发，写起来比较麻烦，这里写成accumulate
        double total_res = 0;
        int effective_num = 0;

        Mat6d H = Mat6d::Zero();
        Vec6d err = Vec6d::Zero();

        for (const auto& idx : index_surf) {
            if (effect_surf[idx]) {
                H += jacob_surf[idx].transpose() * jacob_surf[idx];
                err += -jacob_surf[idx].transpose() * errors_surf[idx];
                effective_num++;
                total_res += errors_surf[idx] * errors_surf[idx];
            }
        }

        for (const auto& idx : index_edge) {
            if (effect_edge[idx]) {
                H += jacob_edge[idx].transpose() * jacob_edge[idx];
                err += -jacob_edge[idx].transpose() * errors_edge[idx];
                effective_num++;
                total_res += errors_edge[idx].norm();
            }
        }

        for (const auto& idx : index_ground) {
            if (effect_ground[idx]) {
                H += jacob_ground[idx].transpose() * jacob_ground[idx];
                err += -jacob_ground[idx].transpose() * errors_ground[idx];
                effective_num++;
                total_res += errors_ground[idx] * errors_ground[idx];
            }
        }

        if (effective_num < options_.min_effective_pts_) {
            LOG(WARNING) << "effective num too small: " << effective_num;
            return pose;
        }

        Vec6d dx = H.inverse() * err;
        pose.so3() = pose.so3() * SO3::exp(dx.head<3>());
        pose.translation() += dx.tail<3>();

        // 更新
        LOG(INFO) << "iter " << iter << " total res: " << total_res << ", eff: " << effective_num
                  << ", mean res: " << total_res / effective_num << ", dxn: " << dx.norm();

        if (dx.norm() < options_.eps_) {
            LOG(INFO) << "converged, dx = " << dx.transpose();
            break;
        }
    }

    estimated_poses_.emplace_back(pose);
    return pose;
}

void LoamLikeOdom::SaveMap(const std::string& path) {
    if (global_map_ && global_map_->empty() == false) 
        pcl::io::savePCDFileBinaryCompressed(path, *global_map_);
}

}  // namespace sad