//
// Created by xiang on 2022/3/15.
//

#ifndef SLAM_IN_AUTO_DRIVING_LIDAR_UTILS_H
#define SLAM_IN_AUTO_DRIVING_LIDAR_UTILS_H

#include <sensor_msgs/LaserScan.h>
#include <sensor_msgs/MultiEchoLaserScan.h>
#include <sensor_msgs/PointCloud.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/point_cloud_conversion.h>

#include <pcl/filters/voxel_grid.h>
#include <pcl_conversions/pcl_conversions.h>

#include "common/point_types.h"
#include "velodyne_msgs/VelodyneScan.h"

/// 雷达扫描的一些消息定义和工具函数
using Scan2d = sensor_msgs::LaserScan;
using MultiScan2d = sensor_msgs::MultiEchoLaserScan;
using PacketsMsg = velodyne_msgs::VelodyneScan;
using PacketsMsgPtr = boost::shared_ptr<PacketsMsg>;

namespace sad {

inline Scan2d::Ptr MultiToScan2d(MultiScan2d::Ptr mscan) {
    Scan2d::Ptr scan(new Scan2d);
    scan->header = mscan->header;
    scan->range_max = mscan->range_max;
    scan->range_min = mscan->range_min;
    scan->angle_increment = mscan->angle_increment;
    scan->angle_max = mscan->angle_max;
    scan->angle_min = mscan->angle_min;
    for (auto r : mscan->ranges) {
        if (r.echoes.empty()) {
            scan->ranges.emplace_back(scan->range_max + 1.0);
        } else {
            scan->ranges.emplace_back(r.echoes[0]);
        }
    }
    for (auto i : mscan->intensities) {
        if (i.echoes.empty()) {
            scan->intensities.emplace_back(0);
        } else {
            scan->intensities.emplace_back(i.echoes[0]);
        }
    }
    scan->scan_time = mscan->scan_time;
    scan->time_increment = mscan->time_increment;

    // limit range max
    scan->range_max = 20.0;
    return scan;
}

/// ROS PointCloud2 转通常的pcl PointCloud
inline CloudPtr PointCloud2ToCloudPtr(sensor_msgs::PointCloud2::Ptr msg) {
    CloudPtr cloud(new PointCloudType);
    pcl::fromROSMsg(*msg, *cloud);
    return cloud;
}

/**
 * 其他类型点云转到PointType点云
 * 用的最多的是全量点云转到XYZI点云
 * @tparam PointT
 * @param input
 * @return
 */
template <typename PointT = FullPointType>
CloudPtr ConvertToCloud(typename pcl::PointCloud<PointT>::Ptr input) {
    CloudPtr cloud(new PointCloudType);
    for (auto& pt : input->points) {
        PointType p;
        p.x = pt.x;
        p.y = pt.y;
        p.z = pt.z;
        p.intensity = pt.intensity;
        cloud->points.template emplace_back(p);
    }
    cloud->width = input->width;
    return cloud;
}

/**
 * @description: 对点云进行体素滤波voxel filter，指定分辨率
 * @param {CloudPtr} cloud
 * @param {float} voxel_size
 * @return {*}
 */
inline CloudPtr VoxelCloud(CloudPtr cloud, float voxel_size = 0.1) {
    CloudPtr output(new PointCloudType); 
    pcl::VoxelGrid<PointType> voxel;                        // 体素滤波器
    voxel.setLeafSize(voxel_size, voxel_size, voxel_size);  // 设置分辨率
    voxel.setInputCloud(cloud);     // 设置待滤波的输入点云
    voxel.filter(*output);          // 执行滤波，将滤波结果存储到output中
    return output;
}

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_LIDAR_UTILS_H
