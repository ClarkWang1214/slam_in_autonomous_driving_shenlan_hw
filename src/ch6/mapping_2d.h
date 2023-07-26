//
// Created by xiang on 2022/3/23.
//

#ifndef SLAM_IN_AUTO_DRIVING_MAPPING_2D_H
#define SLAM_IN_AUTO_DRIVING_MAPPING_2D_H

#include "ch6/frame.h"
#include "common/eigen_types.h"
#include "common/lidar_utils.h"

#include <memory>
#include <opencv2/core.hpp>

namespace sad {

class Submap;       // 子地图
class LoopClosing;  // 回环检测

/**
 * 2D 激光建图的主要类
 */
class Mapping2D {
public:
    /// @brief  初始化
    /// @param with_loop_closing 
    /// @return 
    bool Init(bool with_loop_closing = true);

    /// @brief  处理一帧激光数据scan（单回波）
    /// @param scan 
    /// @return 
    bool ProcessScan(Scan2d::Ptr scan); 

    /// @brief 处理一帧激光数据scan（多回波） 暂时没用到
    /// @param scan 
    /// @return 
    bool ProcessScan(MultiScan2d::Ptr scan);

    /**
     * 显示全局地图
     * @param max_size 全局地图最大长宽
     * @return 全局地图图像
     */
    cv::Mat ShowGlobalMap(int max_size = 500);

   private:
    /// 判定当前帧是否为关键帧
    bool IsKeyFrame();

    /// 增加一个关键帧
    void AddKeyFrame();

    /// 扩展新的submap
    void ExpandSubmap();

    /// 数据成员
    size_t frame_id_ = 0;       // 当前帧id
    size_t keyframe_id_ = 0;    // 当前关键帧id
    size_t submap_id_ = 0;      // 当前子地图id

    bool first_scan_ = true;                                // 是否是第一帧激光数据
    std::shared_ptr<Frame> current_frame_ = nullptr;        // 当前帧
    std::shared_ptr<Frame> last_frame_ = nullptr;           // 上一帧
    SE2 motion_guess_;                                      // 位姿估计

    std::shared_ptr<Frame> last_keyframe_ = nullptr;        // 上一关键帧
    std::shared_ptr<Submap> current_submap_ = nullptr;      // 当前子地图
    std::vector<std::shared_ptr<Submap>> all_submaps_;      // 所有子地图
    
    std::shared_ptr<LoopClosing> loop_closing_ = nullptr;   // 回环检测

    // 参数
    inline static constexpr double keyframe_pos_th_ = 0.3;              // 关键帧位移量
    inline static constexpr double keyframe_ang_th_ = 15 * M_PI / 180;  // 关键帧角度量
};

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_MAPPING_2D_H
