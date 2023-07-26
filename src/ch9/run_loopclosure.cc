//
// Created by xiang on 22-12-19.
//

#include <gflags/gflags.h>
#include <glog/logging.h>
#include "loopclosure.h"

DEFINE_string(config_yaml, "./config/mapping.yaml", "配置文件");
// DEFINE_bool(use_pcl_ndt, false, "use pcl ndt to align?");

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);

    sad::LoopClosure lc(FLAGS_config_yaml);
    // 从keyframe.txt文件中获取关键帧信息
    // 从yaml文件中加载回环检测相关的参数值
    lc.Init();

    // 开始回环检测
    lc.Run();

    return 0;
}