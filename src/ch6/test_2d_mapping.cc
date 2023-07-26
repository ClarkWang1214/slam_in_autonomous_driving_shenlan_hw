//
// Created by xiang on 2022/3/15.
//
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <opencv2/highgui.hpp>

#include "ch6/lidar_2d_utils.h"
#include "ch6/mapping_2d.h"
#include "common/io_utils.h"

DEFINE_string(bag_path, "./dataset/sad/2dmapping/floor1.bag", "数据包路径");
DEFINE_bool(with_loop_closing, false, "是否使用回环检测");

/// 测试2D lidar SLAM

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);

    sad::RosbagIO rosbag_io(fLS::FLAGS_bag_path);
    sad::Mapping2D mapping;

    // std::system("rm -rf ./data/ch6/*"); // 删除ch6文件夹下的所有文件

    // 初始化
    if (mapping.Init(FLAGS_with_loop_closing) == false) 
        return -1;

    rosbag_io.AddScan2DHandle(  FLAGS_scan_name,// 只有floor4.bag是"pavo_scan_bottom"，其它三个"/pavo_scan_bottom"
                                [&](Scan2d::Ptr scan) {     // Lambda函数
                                    return mapping.ProcessScan(scan);   // 2D激光处理
                                }
                            )
                            .Go(); // 遍历rosbag数据包
    // cv::imwrite("./data/ch6/global_map.png", mapping.ShowGlobalMap(2000));  // 保存全局地图
    return 0;
}