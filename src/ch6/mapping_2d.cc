//
// Created by xiang on 2022/3/23.
//

#include "ch6/mapping_2d.h"
#include "ch6/lidar_2d_utils.h"
#include "ch6/loop_closing.h"
#include "ch6/submap.h"

#include <glog/logging.h>
#include <execution>
#include <opencv2/opencv.hpp>

namespace sad {

/**
 * @description: 初始化
 * @param {bool} with_loop_closing  是否使用回环检测
 * @return {*}
 */
bool Mapping2D::Init(bool with_loop_closing) {
    keyframe_id_ = 0;
    current_submap_ = std::make_shared<Submap>(SE2());  // 第一个submap的pose为_R(0),_t(0,0)
    all_submaps_.emplace_back(current_submap_);         // 将第一个submap加入到all_submaps_中

    // 判断是否使用回环检测
    if (with_loop_closing) {
        loop_closing_ = std::make_shared<LoopClosing>();    // 创建回环检测对象，为什么要用make_shared? 为了在多个类中共享这个对象
        loop_closing_->AddNewSubmap(current_submap_);       // 将第一个submap加入到回环检测中
    }

    return true;
}

bool Mapping2D::ProcessScan(MultiScan2d::Ptr scan) { return ProcessScan(MultiToScan2d(scan)); }

/**
 * @description: 处理一帧激光数据scan（单回波）
 * @param {Ptr} scan
 * @return {*}
 */
bool Mapping2D::ProcessScan(Scan2d::Ptr scan) {
    current_frame_ = std::make_shared<Frame>(scan); // 用接收到的2D激光消息构造当前帧
    current_frame_->id_ = frame_id_++;              // 当前帧id加1

    // 判断是否有上一帧
    if (last_frame_) {
        // 若有，则从上一帧中获取位姿估计
        // set pose from last frame
        // current_frame_->pose_ = last_frame_->pose_;
        current_frame_->pose_ = last_frame_->pose_ * motion_guess_; // 当前帧在世界系下的位姿   T_w_c
        current_frame_->pose_submap_ = last_frame_->pose_submap_;   // 当前帧在子地图中的位姿   T_s_c
    }

    // 利用scan matching来匹配地图
    if (!first_scan_) 
        // 第一帧无法匹配，直接加入到occupancy map
        current_submap_->MatchScan(current_frame_);
    

    // current_submap_->AddScanInOccupancyMap(current_frame_);
    first_scan_ = false;

    // 如果机器人发生移动或者转动，按一定距离和角度阈值来取关键帧
    bool is_kf = IsKeyFrame();

    if (is_kf) {
        // 若当前帧是关键帧，则将当前关键帧加入到子地图中
        AddKeyFrame();

        // 将当前关键帧加入到占据栅格地图中，并采用模板法或者BresenhamFilling算法更新占据栅格地图
        current_submap_->AddScanInOccupancyMap(current_frame_);

        // 判断是否使用回环检测
        if (loop_closing_) 
            // 处理回环检测，将当前关键帧加入到回环检测中
            loop_closing_->AddNewFrame(current_frame_);
        
        if (current_submap_->HasOutsidePoints()     // 判断机器人移动范围是否超出了当前子地图
            || (current_submap_->NumFrames()) > 50) // 或者，当前子地图包含的关键帧数量是否超过了一定数量（50）
            // 走出了submap或者单个submap中的关键帧较多，则建立扩展新的子地图。
            ExpandSubmap(); 
    }

    /// 可视化输出
    auto occu_image = current_submap_->GetOccuMap().GetOccupancyGridBlackWhite();
    Visualize2DScan(current_frame_->scan_, current_frame_->pose_, occu_image, Vec3b(0, 0, 255), 1000, 20.0, current_submap_->GetPose());
    cv::putText(occu_image, "submap " + std::to_string(current_submap_->GetId()), cv::Point2f(20, 20), cv::FONT_HERSHEY_COMPLEX, 0.5, cv::Scalar(0, 255, 0));
    cv::putText(occu_image, "keyframes " + std::to_string(current_submap_->NumFrames()), cv::Point2f(20, 50), cv::FONT_HERSHEY_COMPLEX, 0.5, cv::Scalar(0, 255, 0));
    cv::imshow("occupancy map", occu_image);

    auto field_image = current_submap_->GetLikelihood().GetFieldImage();
    Visualize2DScan(current_frame_->scan_, current_frame_->pose_, field_image, Vec3b(0, 0, 255), 1000, 20.0, current_submap_->GetPose());
    cv::imshow("likelihood", field_image);

    /// global map
    if (is_kf) 
        cv::imshow("global map", ShowGlobalMap());

    cv::waitKey(10);

    if (last_frame_)
        motion_guess_ = last_frame_->pose_.inverse() * current_frame_->pose_;

    last_frame_ = current_frame_;

    return true;
}

/**
 * @description: 判定当前帧是否为关键帧
 * @return {*}
 */
bool Mapping2D::IsKeyFrame() {
    // 判断上一帧是否为空，也就是当前帧是否为第一帧
    if (last_keyframe_ == nullptr) 
        return true;    // 第一帧肯定是关键帧

    // 计算上一帧的位姿与当前帧位姿之差： T_c1_c2 = (T_w_c1)^-1 * T_w_c2 = T_c1_w * T_w_c2
    SE2 delta_pose = last_keyframe_->pose_.inverse() * current_frame_->pose_;
    if (delta_pose.translation().norm() > keyframe_pos_th_  // 两帧之间的位移是否大于阈值
        || fabs(delta_pose.so2().log()) > keyframe_ang_th_) // 或 两帧之间的旋转角度的绝对值是否大于阈值
        return true;    // 如果机器人发生移动或者转动，按一定距离和角度阈值来取关键帧

    return false;
}

/**
 * @description: 增加一个关键帧
 * @return {*}
 */
void Mapping2D::AddKeyFrame() {
    LOG(INFO) << "add keyframe " << keyframe_id_;
    current_frame_->keyframe_id_ = keyframe_id_++;  // 当前帧的关键帧id加1
    current_submap_->AddKeyFrame(current_frame_);   // 将当前帧加入到当前子地图中
    last_keyframe_ = current_frame_;                // 更新上一关键帧
}

/**
 * @description: 建立新的子地图submap
 * @return {*}
 */
void Mapping2D::ExpandSubmap() {
    
    if (loop_closing_) 
        // 当前submap作为历史地图放入loop closing
        loop_closing_->AddFinishedSubmap(current_submap_); 

    // 将当前submap替换成新的
    auto last_submap = current_submap_; // 保存当前submap为上一个

    // debug调试时使用，保存用于查看所有子地图
    cv::imwrite("./data/ch6/submap_" + std::to_string(last_submap->GetId()) + ".png", 
                last_submap->GetOccuMap()       // 获取上一个子地图的占据栅格地图对象
                .GetOccupancyGridBlackWhite()); // 根据占据栅格值，绘制黑白灰三色的占据栅格地图

    // 生成新的子地图对象，以当前帧为中心，位姿取为T_ws = T_wc
    current_submap_ = std::make_shared<Submap>(current_frame_->pose_);
    current_submap_->SetId(++submap_id_); // 子地图id加1

    // 当前帧在子地图中的位姿置零
    current_frame_->pose_submap_ = SE2();
    current_submap_->AddKeyFrame(current_frame_);   // 将当前帧加入到子地图中

    // 此时，新的子地图中没有任何数据，可以将旧的子地图中最近（10个）关键帧拷贝至新的子地图中，不让新的submap显得太空
    current_submap_->SetOccuFromOtherSubmap(last_submap);

    current_submap_->AddScanInOccupancyMap(current_frame_); // 将当前帧加入到占据栅格地图中

    // 添加当前子地图到所有子地图中
    all_submaps_.emplace_back(current_submap_);

    if (loop_closing_) 
        // 加入到回环检测中
        loop_closing_->AddNewSubmap(current_submap_);

    LOG(INFO) << "create submap " << current_submap_->GetId()
              << " with pose: " << current_submap_->GetPose().translation().transpose() << ", "
              << current_submap_->GetPose().so2().log();
}

cv::Mat Mapping2D::ShowGlobalMap(int max_size) {
    //// TODO 全局地图固定大小，使用动态分辨率
    Vec2f top_left = Vec2f(999999, 999999);
    Vec2f bottom_right = Vec2f(-999999, -999999);

    const float submap_resolution = 20.0;  // 子地图分辨率（1米多少个像素）
    const float submap_size = 50.0;        // 单个submap大小

    /// 计算全局地图物理边界
    for (auto m : all_submaps_) {
        Vec2d c = m->GetPose().translation();
        if (top_left[0] > c[0] - submap_size / 2) {
            top_left[0] = c[0] - submap_size / 2;
        }
        if (top_left[1] > c[1] - submap_size / 2) {
            top_left[1] = c[1] - submap_size / 2;
        }

        if (bottom_right[0] < c[0] + submap_size / 2) {
            bottom_right[0] = c[0] + submap_size / 2;
        }
        if (bottom_right[1] < c[1] + submap_size / 2) {
            bottom_right[1] = c[1] + submap_size / 2;
        }
    }

    if (top_left[0] > bottom_right[0] || top_left[1] > bottom_right[1]) {
        return cv::Mat();
    }

    /// 全局地图物理中心
    Vec2f global_center = Vec2f((top_left[0] + bottom_right[0]) / 2.0, (top_left[1] + bottom_right[1]) / 2.0);
    float phy_width = bottom_right[0] - top_left[0];   // 物理尺寸
    float phy_height = bottom_right[1] - top_left[1];  // 物理尺寸
    float global_map_resolution = 0;

    if (phy_width > phy_height) {
        global_map_resolution = max_size / phy_width;
    } else {
        global_map_resolution = max_size / phy_height;
    }

    Vec2f c = global_center;
    int c_x = global_center[0] * global_map_resolution;
    int c_y = global_center[1] * global_map_resolution;
    global_center = Vec2f(c_x / global_map_resolution, c_y / global_map_resolution);  // 全局地图图像中心

    int width = int((bottom_right[0] - top_left[0]) * global_map_resolution + 0.5);
    int height = int((bottom_right[1] - top_left[1]) * global_map_resolution + 0.5);

    Vec2f center_image = Vec2f(width / 2, height / 2);
    cv::Mat output_image(height, width, CV_8UC3, cv::Scalar(127, 127, 127));

    std::vector<Vec2i> render_data;
    render_data.reserve(width * height);
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            render_data.emplace_back(Vec2i(x, y));
        }
    }

    std::for_each(std::execution::par_unseq, render_data.begin(), render_data.end(), [&](const Vec2i& xy) {
        int x = xy[0], y = xy[1];
        Vec2f pw = (Vec2f(x, y) - center_image) / global_map_resolution + c;  // 世界坐标

        for (auto& m : all_submaps_) {
            Vec2f ps = m->GetPose().inverse().cast<float>() * pw;  // in submap
            Vec2i pt = (ps * submap_resolution + Vec2f(500, 500)).cast<int>();

            if (pt[0] < 0 || pt[0] >= 1000 || pt[1] < 0 || pt[1] >= 1000) {
                continue;
            }

            uchar value = m->GetOccuMap().GetOccupancyGrid().at<uchar>(pt[1], pt[0]);
            if (value > 127) {
                if (m == current_submap_) {
                    output_image.at<cv::Vec3b>(y, x) = cv::Vec3b(235, 250, 230);
                } else {
                    output_image.at<cv::Vec3b>(y, x) = cv::Vec3b(255, 255, 255);
                }
                break;
            } else if (value < 127) {
                if (m == current_submap_) {
                    output_image.at<cv::Vec3b>(y, x) = cv::Vec3b(230, 20, 30);
                } else {
                    output_image.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 0, 0);
                }
                break;
            }
        }
    });

    for (auto& m : all_submaps_) {
        /// submap pose 在全局地图中的投影
        SE2f submap_pose = m->GetPose().cast<float>();
        Vec2f submap_center = submap_pose.translation();
        Vec2f submap_xw = submap_pose * Vec2f(1.0, 0);
        Vec2f submap_yw = submap_pose * Vec2f(0, 1.0);

        Vec2f center_map = (submap_center - global_center) * global_map_resolution + center_image;
        Vec2f x_map = (submap_xw - global_center) * global_map_resolution + center_image;
        Vec2f y_map = (submap_yw - global_center) * global_map_resolution + center_image;

        // x轴和y轴
        cv::line(output_image, cv::Point2f(center_map.x(), center_map.y()), cv::Point2f(x_map.x(), x_map.y()),
                 cv::Scalar(0, 0, 255), 2);
        cv::line(output_image, cv::Point2f(center_map.x(), center_map.y()), cv::Point2f(y_map.x(), y_map.y()),
                 cv::Scalar(0, 255, 0), 2);
        cv::putText(output_image, std::to_string(m->GetId()), cv::Point2f(center_map.x() + 10, center_map.y() - 10),
                    cv::FONT_HERSHEY_COMPLEX, 0.5, cv::Scalar(255, 0, 0));

        // 轨迹
        for (const auto& frame : m->GetFrames()) {
            Vec2f p_map =
                (frame->pose_.translation().cast<float>() - global_center) * global_map_resolution + center_image;
            cv::circle(output_image, cv::Point2f(p_map.x(), p_map.y()), 1, cv::Scalar(0, 0, 255), 1);
        }
    }

    if (loop_closing_) {
        /// 回环检测的pose graph
        auto loops = loop_closing_->GetLoops();
        for (auto lc : loops) {
            auto first_id = lc.first.first;
            auto second_id = lc.first.second;

            Vec2f c1 = all_submaps_[first_id]->GetPose().translation().cast<float>();
            Vec2f c2 = all_submaps_[second_id]->GetPose().translation().cast<float>();

            Vec2f c1_map = (c1 - global_center) * global_map_resolution + center_image;
            Vec2f c2_map = (c2 - global_center) * global_map_resolution + center_image;

            cv::line(output_image, cv::Point2f(c1_map.x(), c1_map.y()), cv::Point2f(c2_map.x(), c2_map.y()),
                     cv::Scalar(255, 0, 0), 2);
        }
    }

    return output_image;
}

}  // namespace sad