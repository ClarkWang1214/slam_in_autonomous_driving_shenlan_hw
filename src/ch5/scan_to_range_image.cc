// // Created by xiang on 2021/8/9.
//

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <opencv2/opencv.hpp>

using PointType = pcl::PointXYZI;
using PointCloudType = pcl::PointCloud<PointType>;

DEFINE_string(pcd_path, "./data/ch5/scan_example.pcd", "点云文件路径");
DEFINE_double(azimuth_resolution_deg, 0.3, "方位角分辨率（度）");
DEFINE_int32(elevation_rows, 16, "俯仰角对应的行数");
DEFINE_double(elevation_range, 15.0, "俯仰角范围");  // VLP-16 上下各15度范围
DEFINE_double(lidar_height, 1.128, "雷达安装高度");

/**
 * @description: PCL RANSAC 分割地面点云
 * @param {Ptr} cloud 点云
 * @return {*}
 */
void ExtractGroundPoints_RANSAC(PointCloudType::Ptr cloud) {
    //创建分割时所需要的模型系数对象coefficients及存储内点的点索引集合对象inliers。
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    pcl::PointIndices::Ptr inliers(new pcl::PointIndices());
    // 创建分割对象
    pcl::SACSegmentation<pcl::PointXYZ> seg;    // 可选择配置，设置模型系数需要优化
    seg.setOptimizeCoefficients(true);          // 设置对估计的模型参数进行优化处理
    seg.setModelType(pcl::SACMODEL_PLANE);      // 设置分割的模型类型为平面
    seg.setMethodType(pcl::SAC_RANSAC);         // 设置用ransac随机参数估计方法
    seg.setMaxIterations(1000);
    seg.setDistanceThreshold(0.15);             // 平面厚度距离阈值，表示距离平面多少米之内的点都算平面内点inlier。
    seg.setInputCloud(cloud);                   // 输入点云
    seg.segment(*inliers, *coefficients);       // 实现分割，并存储分割结果到点集合inliers及存储平面模型系数coefficients

    PointCloudType::Ptr cloud_ground(new PointCloudType),       // 地面点云
                        cloud_no_ground(new PointCloudType);    // 非地面点云
                        
    pcl::ExtractIndices<pcl::PointXYZ> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(inliers);
    extract.setNegative(true);  // true表示剔除地面点
    extract.filter(*cloud_no_ground);
    extract.setNegative(false); // false表示保留平面内点
    extract.filter(*cloud_ground);

    pcl::io::savePCDFileBinaryCompressed("./data/ch5/scan_example_ground.pcd", *cloud_ground);
    pcl::io::savePCDFileBinaryCompressed("./data/ch5/scan_example_no_ground.pcd", *cloud_no_ground);
}

/**
 * @description: 生成距离图
 * @param {Ptr} cloud 点云
 * @return {*}
 */
void GenerateRangeImage(PointCloudType::Ptr cloud) {
    int image_cols = int(360 / FLAGS_azimuth_resolution_deg);  // 水平为360度，按分辨率切分即可，360/0.3=1200
    int image_rows = FLAGS_elevation_rows;                     // 图像行数（激光雷达线数）固定 ，16
    LOG(INFO) << "range image: " << image_rows << "x" << image_cols; // 生成的距离图大小为：16x1200

    // 我们生成一个HSV图像以更好地显示图像
    cv::Mat image(image_rows, image_cols, CV_8UC3, cv::Scalar(0, 0, 0));

    double ele_resolution = FLAGS_elevation_range * 2 / FLAGS_elevation_rows;  // elevation俯仰角分辨率

    for (const auto& pt : cloud->points) {
        double azimuth = atan2(pt.y, pt.x) * 180 / M_PI;    // 根据方位角A，arctan(y/x)  horizontal angle
        // double range = sqrt(pt.x * pt.x + pt.y * pt.y);     
        double range = sqrt(pt.x * pt.x + pt.y * pt.y + pt.z * pt.z);     // 计算距离r，对应公式（5.2）
        double elevation = asin((pt.z - FLAGS_lidar_height) / range) * 180 / M_PI;  // 计算俯仰角E，arcsin(z/r)

        // keep in 0~360
        if (azimuth < 0) {
            azimuth += 360; // 保证方位角在0~360度之间
        }

        // 计算当前点在距离图中的行索引和列索引
        int x = int(azimuth / FLAGS_azimuth_resolution_deg);                      // 方位角除以分辨率，得到列索引，
        int y = int((elevation + FLAGS_elevation_range) / ele_resolution + 0.5);  // 列

        if (x >= 0 && x < image.cols && y >= 0 && y < image.rows) {
            // 将满足条件的点赋予颜色，得到距离图
            image.at<cv::Vec3b>(y, x) = cv::Vec3b(uchar(range / 100 * 255.0), 255, 127);
        }
    }

    // 沿Y轴翻转，因为我们希望Z轴朝上时Y朝上
    cv::Mat image_flipped;
    cv::flip(image, image_flipped, 0);

    // hsv to rgb
    cv::Mat image_rgb;
    cv::cvtColor(image_flipped, image_rgb, cv::COLOR_HSV2BGR); // 将HSV图像转换为RGB图像
    cv::imwrite("./range_image.png", image_rgb);
    
    // cv::Mat image_rgb_resized;
    // cv::resize(image_rgb, image_rgb_resized, cv::Size(300, 16), 0, 0, cv::INTER_NEAREST); // 将图像缩小一半

    // // cv::Rect m_select = cv::Rect(0, 0, 256, 16); // 将图像裁剪为原始大小
    // // cv::Mat image_rgb_resized2 = image_rgb_resized(m_select);
    // cv::imwrite("./range_image_resized.png", image_rgb_resized);
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

    // 地面点云提取 PCL RANSAC 拟合
    ExtractGroundPoints_RANSAC(cloud);

    // 生成距离图
    // GenerateRangeImage(cloud);

    return 0;
}
