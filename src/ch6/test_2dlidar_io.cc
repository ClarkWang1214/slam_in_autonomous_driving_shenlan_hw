//
// Created by xiang on 2022/3/15.
//
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <opencv2/highgui.hpp>

#include "ch6/lidar_2d_utils.h"
#include "common/io_utils.h"

DEFINE_string(bag_path, "./dataset/sad/2dmapping/test_2d_lidar.bag", "数据包路径");

/// 测试从rosbag中读取2d scan并plot的结果

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);

    sad::RosbagIO rosbag_io(fLS::FLAGS_bag_path);
    rosbag_io.AddScan2DHandle("/pavo_scan_bottom", // 2D激光话题
                         [](Scan2d::Ptr scan) { // lambda匿名函数
                             cv::Mat image; // 用于显式图像
                             // 调用2D激光可视化函数，绘制在image图像上
                             sad::Visualize2DScan(scan, SE2(), image, Vec3b(255, 0, 0));
                             cv::imshow("scan", image);
                             cv::waitKey(20);
                             return true;
                         }).Go();

    return 0;
}