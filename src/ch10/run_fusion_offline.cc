//
// Created by xiang on 22-12-20.
//

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <yaml-cpp/yaml.h>

#include <chrono>

#include "common/io_utils.h"
#include "fusion.h"



DEFINE_string(config_yaml, "./config/mapping.yaml", "配置文件");

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);

    sad::Fusion fusion(FLAGS_config_yaml);
    if (!fusion.Init()) {
        return -1;
    }

    auto yaml = YAML::LoadFile(FLAGS_config_yaml);
    auto bag_path = yaml["bag_path"].as<std::string>();
    sad::RosbagIO rosbag_io(bag_path, sad::DatasetType::NCLT);

    auto t1 = std::chrono::steady_clock::now();

    /// 把各种消息交给fusion
    rosbag_io
        .AddAutoRTKHandle([&fusion](GNSSPtr gnss) {
            fusion.ProcessRTK(gnss);
            return true;
        })
        .AddAutoPointCloudHandle([&](sensor_msgs::PointCloud2::Ptr cloud) -> bool {
            fusion.ProcessPointCloud(cloud);
            return true;
        })
        .AddImuHandle([&](IMUPtr imu) {
            fusion.ProcessIMU(imu);
            return true;
        })
        .Go();

    auto t2 = std::chrono::steady_clock::now();
    auto time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1).count() * 1000;
    LOG(INFO) << "time used of this dataset = " << time_used;

    LOG(INFO) << "done.";
    sleep(10);
    return 0;
}