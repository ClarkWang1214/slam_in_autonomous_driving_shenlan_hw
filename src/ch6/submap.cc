//
// Created by xiang on 2022/3/23.
//

#include "ch6/submap.h"
#include <glog/logging.h>

namespace sad {

/**
 * @description: 把另一个submap中的占据栅格复制到本地图中
 * @param {shared_ptr<Submap>} other
 * @return {*}
 */
void Submap::SetOccuFromOtherSubmap(std::shared_ptr<Submap> other) {
    auto frames_in_other = other->GetFrames();  // 取出另一个submap中的所有帧点云
    // 取出vector容器中的最近10个帧（倒数10个）
    for (size_t i = frames_in_other.size() - 10; i < frames_in_other.size(); ++i) {
        if (i > 0) 
            occu_map_.AddLidarFrame(frames_in_other[i]);    // 将帧点云加入到栅格地图中
    }
    field_.SetFieldImageFromOccuMap(occu_map_.GetOccupancyGrid());
}

/**
 * @description: 将当前帧scan扫描数据与当前子地图进行匹配，得到scan在当前子地图中的位姿，T_s_c
 * @param {shared_ptr<Frame>} frame
 * @return {*}
 */
bool Submap::MatchScan(std::shared_ptr<Frame> frame) {
    field_.SetSourceScan(frame->scan_);             // 设置似然场的源点云    
    field_.AlignG2O(frame->pose_submap_);           // 调用基于g2o的似然场方法实现当前帧与子地图似然场的匹配，得到当前帧在子地图中的位姿 T_s_c 
    frame->pose_ = pose_ * frame->pose_submap_;     // 当前帧在世界坐标系中的位姿 T_w_c = T_w_s * T_s_c    

    return true;
}

/**
 * @description: 在栅格地图中增加一个关键帧
 * @param {shared_ptr<Frame>} frame
 * @return {*}
 */
void Submap::AddScanInOccupancyMap(std::shared_ptr<Frame> frame) {
    occu_map_.AddLidarFrame(frame, OccupancyMap::GridMethod::MODEL_POINTS);  // 更新栅格地图中的格子，模板法
    field_.SetFieldImageFromOccuMap(occu_map_.GetOccupancyGrid());           // 更新似然场函数图像
}

/**
 * @description: 判断是否有落在占据栅格外部
 * @param {return} occu_map_
 * @return {*}
 */
bool Submap::HasOutsidePoints() const { return occu_map_.HasOutsidePoints(); }

void Submap::SetPose(const SE2& pose) {
    pose_ = pose;
    occu_map_.SetPose(pose);
    field_.SetPose(pose);
}

void Submap::UpdateFramePoseWorld() {
    for (auto& frame : frames_) {
        frame->pose_ = pose_ * frame->pose_submap_;
    }
}

}  // namespace sad