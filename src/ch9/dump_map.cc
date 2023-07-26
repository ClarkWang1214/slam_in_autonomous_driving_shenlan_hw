//
// Created by xiang on 22-12-7.
//

#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_double(voxel_size, 0.1, "导出地图分辨率");
DEFINE_string(pose_source, "lidar", "使用的pose来源:lidar / rtk / opti1 / opti2");
DEFINE_string(dump_to, "./data/ch9/", "导出的目标路径");

#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>

#include "keyframe.h"
#include "common/point_cloud_utils.h"

/**
 * 将keyframes.txt中的地图和点云合并为一个pcd
 */
int main(int argc, char** argv) {
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;

    using namespace sad;
    std::map<IdType, KFPtr> keyframes;  // 关键帧集合
    if (!LoadKeyFrames("./data/ch9/keyframes.txt", keyframes)) { // 从文件中读取关键帧
        LOG(ERROR) << "failed to load keyframes.txt";
        return -1;
    }

    // 检查关键帧集合是否为空
    if (keyframes.empty()) {
        LOG(INFO) << "keyframes are empty";
        return 0;
    }

    // dump kf cloud and merge
    // 关键帧点云合并成整个地图文件
    LOG(INFO) << "merging";
    CloudPtr global_cloud(new PointCloudType);      // 新建一个全局地图点云对象
    pcl::VoxelGrid<PointType> voxel_grid_filter;    // 体素滤波器
    float resolution = FLAGS_voxel_size;            // 体素滤波分辨率
    voxel_grid_filter.setLeafSize(resolution, resolution, resolution);

    int cnt = 0;
    // 按照关键帧id从小到大遍历关键帧集合
    for (auto& kfp : keyframes) {
        auto kf = kfp.second;   // 当前关键帧
        SE3 pose;
        // 按照指定的pose来源选择位姿
        if (FLAGS_pose_source == "rtk") 
            pose = kf->rtk_pose_;
        else if (FLAGS_pose_source == "lidar") 
            pose = kf->lidar_pose_;
        else if (FLAGS_pose_source == "opti1") 
            pose = kf->opti_pose_1_;
        else if (FLAGS_pose_source == "opti2") 
            pose = kf->opti_pose_2_;

        // 从指定路径中加载以当前关键帧id命名的点云pcd文件
        kf->LoadScan("./data/ch9/");

        CloudPtr cloud_trans(new PointCloudType);
        // 【将当前帧雷达点云转换到全局坐标系下】，变换位姿采用上面指定的位姿源
        pcl::transformPointCloud(*kf->cloud_, *cloud_trans, pose.matrix());

        // voxel size
        CloudPtr kf_cloud_voxeled(new PointCloudType);
        voxel_grid_filter.setInputCloud(cloud_trans);   // 设置输入点云
        voxel_grid_filter.filter(*kf_cloud_voxeled);    // 体素滤波后的点云

        // 将当前帧【点云拼接】到全局点云中，+=操作符重载
        *global_cloud += *kf_cloud_voxeled; 
        kf->cloud_ = nullptr; // 清除内存，指针置空

        LOG(INFO) << "merging " << cnt << " in " << keyframes.size() << ", pts: " << kf_cloud_voxeled->size()
                  << " global pts: " << global_cloud->size();
        cnt++;  // 计数加一
    }

    if (!global_cloud->empty()) {
        // 保存合并后的全局地图点云到本地
        sad::SaveCloudToFile(FLAGS_dump_to + "/map.pcd", *global_cloud);
    }

    LOG(INFO) << "done.";
    return 0;
}
