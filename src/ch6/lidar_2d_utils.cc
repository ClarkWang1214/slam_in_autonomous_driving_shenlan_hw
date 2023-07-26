//
// Created by xiang on 2022/3/15.
//

#include "ch6/lidar_2d_utils.h"
#include <opencv2/imgproc.hpp>

namespace sad {

/**
 * @description: 2D激光数据可视化函数
 * @param {Ptr} scan            当前扫描一圈的2D激光消息数据
 * @param {SE2&} pose           当前scan的位姿（本体系下的扫描数据，输入位姿给定为原点）R(0), t(0,0)
 * @param {Mat&} image          将2D激光数据绘制到图像上
 * @param {Vec3b&} color        颜色
 * @param {int} image_size      图像大小
 * @param {float} resolution    分辨率
 * @param {SE2&} pose_submap    子地图的位姿
 * @return {*}
 */
void Visualize2DScan(Scan2d::Ptr scan, const SE2& pose, cv::Mat& image, const Vec3b& color, int image_size,
                     float resolution, const SE2& pose_submap) {
    if (image.data == nullptr) 
        image = cv::Mat(image_size, image_size, CV_8UC3, cv::Vec3b(255, 255, 255));

    // 遍历扫描一圈的所有点，点的个数取决于激光的角分辨率
    for (size_t i = 0; i < scan->ranges.size(); ++i) {
        // 判断每个点的距离是否越界
        if (scan->ranges[i] < scan->range_min || scan->ranges[i] > scan->range_max) 
            continue;

        // 获取第i个点与车辆自身的角度
        double real_angle = scan->angle_min + i * scan->angle_increment; 

        // 将雷达的距离和角度转换为笛卡尔坐标
        double x = scan->ranges[i] * std::cos(real_angle);  
        double y = scan->ranges[i] * std::sin(real_angle);

        // 判断角度是否越界
        if (real_angle < scan->angle_min + 30 * M_PI / 180.0 || real_angle > scan->angle_max - 30 * M_PI / 180.0) 
            continue;
        
        // 
        Vec2d psubmap = pose_submap.inverse() * (pose * Vec2d(x, y));

        int image_x = int(psubmap[0] * resolution + image_size / 2);
        int image_y = int(psubmap[1] * resolution + image_size / 2);
        if (image_x >= 0 && image_x < image.cols && image_y >= 0 && image_y < image.rows) 
            image.at<cv::Vec3b>(image_y, image_x) = cv::Vec3b(color[0], color[1], color[2]);
    }

    // 同时画出pose自身所在位置
    Vec2d pose_in_image = pose_submap.inverse() * (pose.translation()) * double(resolution) + Vec2d(image_size / 2, image_size / 2);
    cv::circle(image, cv::Point2f(pose_in_image[0], pose_in_image[1]), 5, cv::Scalar(color[0], color[1], color[2]), 2);
}

}  // namespace sad