//
// Created by xiang on 2022/7/25.
//

#ifndef SLAM_IN_AUTO_DRIVING_LOAM_LIKE_ODOM_H
#define SLAM_IN_AUTO_DRIVING_LOAM_LIKE_ODOM_H

#include "ch5/kdtree.h"
#include "ch7/icp_3d.h"
#include "common/eigen_types.h"
#include "common/point_types.h"
#include "tools/pcl_map_viewer.h"

#include <deque>

namespace sad {
class FeatureExtraction;

/**
 * 类loam的里程计方法
 * 首先使用feature extraction提出一个点云中的边缘点和平面点
 * 然后分别对于edge点和surface点，使用不同的ICP方法
 */
class LoamLikeOdom {
public:
    struct Options {
        Options() {}

        int min_edge_pts_ = 20;               // 最小边缘点数
        int min_ground_pts_ = 20;             // 最小地面点数 【新增】
        int min_surf_pts_ = 20;               // 最小平面点数
        double kf_distance_ = 1.0;            // 关键帧距离
        double kf_angle_deg_ = 15;            // 旋转角度
        int num_kfs_in_local_map_ = 30;       // 局部地图含有多少个关键帧
        bool display_realtime_cloud_ = true;  // 是否显示实时点云

        // ICP 参数
        int max_iteration_ = 5;             // 最大迭代次数
        double max_plane_distance_ = 0.05;  // 平面最近邻查找时阈值
        double max_ground_distance_ = 0.05;  // 地面最近邻查找时阈值【新增】
        double max_line_distance_ = 0.5;    // 点线最近邻查找时阈值
        int min_effective_pts_ = 10;        // 最近邻点数阈值
        double eps_ = 1e-3;                 // 收敛判定条件

        bool use_edge_points_ = false;  // 是否使用边缘点
        bool use_surf_points_ = false;  // 是否使用平面点
        bool use_ground_points_ = true; // 是否使用地面点 【新增】
    };

    explicit LoamLikeOdom(Options options = Options());

    /**
     * 往里程计中添加一个点云，内部会分为角点和平面点
     * @param pcd_edge
     * @param pcd_surf
     */
    void ProcessPointCloud(FullCloudPtr full_cloud);

    void SaveMap(const std::string& path);

private:
    /// 与局部地图进行配准
    SE3 AlignWithLocalMap(CloudPtr edge, CloudPtr surf, CloudPtr ground);

    /// 判定是否为关键帧
    bool IsKeyframe(const SE3& current_pose);

    Options options_;

    int cnt_frame_ = 0;
    int last_kf_id_ = 0;

    CloudPtr local_map_edge_ = nullptr, 
             local_map_surf_ = nullptr;  // 局部地图的local map
             
    CloudPtr local_map_ground_ = nullptr; // 局部地面地图的 local map  【新增】


    std::vector<SE3> estimated_poses_;    // 所有估计出来的pose，用于记录轨迹和预测下一个帧
    SE3 last_kf_pose_;                    // 上一关键帧的位姿
    std::deque<CloudPtr> edges_,    // 缓存的角点和平面点
                         surfs_,
                         ground_;  // 缓存的地面点 【新增】

    CloudPtr global_map_ = nullptr;  // 用于保存的全局地图

    std::shared_ptr<FeatureExtraction> feature_extraction_ = nullptr;

    std::shared_ptr<PCLMapViewer> viewer_ = nullptr;
    KdTree kdtree_edge_, 
           kdtree_surf_,
           kdtree_ground_; // 地面点云的kd树 【新增】
};

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_LOAM_LIKE_ODOM_H
