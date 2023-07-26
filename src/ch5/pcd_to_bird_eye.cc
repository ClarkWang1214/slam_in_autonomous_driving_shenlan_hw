//
// Created by xiang on 2021/8/9.
//

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using PointType = pcl::PointXYZI;
using PointCloudType = pcl::PointCloud<PointType>;

DEFINE_string(pcd_path, "./data/ch5/map_example.pcd", "点云文件路径");
DEFINE_double(image_resolution, 0.1, "俯视图分辨率");
DEFINE_double(min_z, 0.2, "俯视图最低高度");
DEFINE_double(max_z, 2.5, "俯视图最高高度");

/// 本节演示如何将一个pcd格式点云转换为俯视图像
void GenerateBEVImage(PointCloudType::Ptr cloud) {
    // 【计算点云边界】
    // std::minmax_element返回一个pair，包含最小值和最大值的迭代器
    auto minmax_x = std::minmax_element(cloud->points.begin(), cloud->points.end(), 
                                        [](const PointType& p1, const PointType& p2) { return p1.x < p2.x; }); // lambda表达式
    auto minmax_y = std::minmax_element(cloud->points.begin(), cloud->points.end(),
                                        [](const PointType& p1, const PointType& p2) { return p1.y < p2.y; }); // lambda表达式
    // x和y的最小值和最大值
    double min_x = minmax_x.first->x;   
    double max_x = minmax_x.second->x;
    double min_y = minmax_y.first->y;   
    double max_y = minmax_y.second->y;

    const double inv_r = 1.0 / FLAGS_image_resolution;      // 图像分辨率的倒数

    const int image_rows = int((max_y - min_y) * inv_r);    // 将最大值和最小值的差值乘以分辨率的倒数，得到图像的行数、列数
    const int image_cols = int((max_x - min_x) * inv_r);    

    float x_center = 0.5 * (max_x + min_x);                 // 点云中心坐标
    float y_center = 0.5 * (max_y + min_y);                 
    float x_center_image = image_cols / 2;                  // 图像中心坐标
    float y_center_image = image_rows / 2;

    // 生成图像
    cv::Mat image(image_rows, image_cols, CV_8UC3, cv::Scalar(255, 255, 255));

    // 遍历点云，按照设置的分辨率计算点云中每个点(x,y,z)在图像中的位置(u,v)
    for (const auto& pt : cloud->points) {
        int x = int((pt.x - x_center) * inv_r + x_center_image);
        int y = int((pt.y - y_center) * inv_r + y_center_image);
        // 认为车辆周围一定高度范围内的障碍物（默认为0.2到2.5米）是有效的，再把这些障碍物信息放置在俯视图中，就得到了一张俯视图像
        if (x < 0 || x >= image_cols || y < 0 || y >= image_rows || pt.z < FLAGS_min_z || pt.z > FLAGS_max_z) {
            // 若点云在图像中的位置越界，或者点云的高度不在有效范围内，则跳过该点
            continue;
        }
        // 给符合高度条件的点赋予颜色
        image.at<cv::Vec3b>(y, x) = cv::Vec3b(227, 143, 79);
    }
    // 保存图像到本地
    cv::imwrite("./bev.png", image);
}

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);

    if (FLAGS_pcd_path.empty()) {
        LOG(ERROR) << "pcd path is empty";
        return -1;
    }

    // 读取点云
    PointCloudType::Ptr cloud(new PointCloudType);
    pcl::io::loadPCDFile(FLAGS_pcd_path, *cloud);

    if (cloud->empty()) {
        LOG(ERROR) << "cannot load cloud file";
        return -1;
    }

    LOG(INFO) << "cloud points: " << cloud->size();
    GenerateBEVImage(cloud);

    return 0;
}
