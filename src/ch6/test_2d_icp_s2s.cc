//
// Created by xiang on 2022/3/15.
//
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <opencv2/highgui.hpp>

#include "ch6/icp_2d.h"
#include "ch6/lidar_2d_utils.h"
#include "common/io_utils.h"

DEFINE_string(bag_path, "./dataset/sad/2dmapping/floor1.bag", "数据包路径");
DEFINE_string(method, "point2point", "2d icp方法：point2point/point2plane");

/// 测试从rosbag中读取2d scan并plot的结果
/// 通过选择method来确定使用点到点或点到面的ICP
int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);

    sad::RosbagIO rosbag_io(fLS::FLAGS_bag_path);
    Scan2d::Ptr last_scan = nullptr, current_scan = nullptr;

    /// 我们将上一个scan与当前scan进行配准
    rosbag_io
        .AddScan2DHandle("/pavo_scan_bottom",
                         [&](Scan2d::Ptr scan) {
                             current_scan = scan;

                             if (last_scan == nullptr) {
                                 last_scan = current_scan;
                                 return true;
                             }

                             sad::Icp2d icp;
                             icp.SetTarget(last_scan);
                             icp.SetSource(current_scan);

                             SE2 pose;
                             if (fLS::FLAGS_method == "point2point") {
                                 icp.AlignGaussNewton(pose);
                             } else if (fLS::FLAGS_method == "point2line") {
                                 icp.AlignGaussNewtonPoint2Plane(pose);
                             } else if (fLS::FLAGS_method == "point2point_g2o") {
                                 LOG(INFO) << "icp.AlignG2OP2P(pose): ";
                                 icp.AlignG2OP2P(pose);     // 【新增】
                             } else if (fLS::FLAGS_method == "point2point_g2o_2") {
                                 LOG(INFO) << "icp.AlignG2OP2P_2(pose): ";
                                 icp.AlignG2OP2P_2(pose);   // 【新增】
                             } else if (fLS::FLAGS_method == "point2point_g2o_3") {
                                 LOG(INFO) << "icp.AlignG2OP2P_3(pose): ";
                                 icp.AlignG2OP2P_3(pose);   // 【新增】
                             } else if (fLS::FLAGS_method == "point2line_g2o") {
                                 LOG(INFO) << "icp.AlignG2OP2L(pose): ";
                                 icp.AlignG2OP2L(pose);     // 【新增】
                             } else if (fLS::FLAGS_method == "point2line_g2o_2") {
                                 LOG(INFO) << "icp.AlignG2OP2L_2(pose): ";
                                 icp.AlignG2OP2L_2(pose);   // 【新增】
                             }

                             cv::Mat image;
                             sad::Visualize2DScan(last_scan, SE2(), image, Vec3b(255, 0, 0));    // target是蓝的
                             sad::Visualize2DScan(current_scan, pose, image, Vec3b(0, 0, 255));  // source是红的
                             cv::imshow("scan", image);
                             cv::waitKey(20);

                             last_scan = current_scan;
                             return true;
                         })
        .Go();

    return 0;
}