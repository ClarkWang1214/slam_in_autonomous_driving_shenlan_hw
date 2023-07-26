//
// Created by xiang on 2022/3/23.
//

#ifndef SLAM_IN_AUTO_DRIVING_SUBMAP_H
#define SLAM_IN_AUTO_DRIVING_SUBMAP_H

#include "ch6/frame.h"
#include "ch6/likelihood_filed.h"
#include "ch6/occupancy_map.h"

namespace sad {

/**
 * 子地图类
 * 子地图关联到若干个关键帧，也会维护自己的栅格地图与似然场
 * 往子地图中添加关键帧时，会更新它的栅格地图与似然场
 * 子地图有自己的pose，记为 Tws，每个frame的世界位姿可以由子地图的pose乘frame在子地图中的相对pose得到
 */
class Submap {
public:
    Submap(const SE2& pose) : pose_(pose) {
        occu_map_.SetPose(pose_);
        field_.SetPose(pose_);
    }

    /// 把另一个submap中的占据栅格复制到本地图中
    void SetOccuFromOtherSubmap(std::shared_ptr<Submap> other);

    /// 将frame与本submap进行匹配，计算frame->pose
    bool MatchScan(std::shared_ptr<Frame> frame);

    /// 判定当前的scan是否有位于submap外部的点
    bool HasOutsidePoints() const;

    /// 在栅格地图中增加一个帧
    void AddScanInOccupancyMap(std::shared_ptr<Frame> frame);

    /// 在子地图中增加一个关键帧
    void AddKeyFrame(std::shared_ptr<Frame> frame) { frames_.emplace_back(frame); }

    /// 当子地图的位姿更新时，重设每个frame的世界位姿
    void UpdateFramePoseWorld();

    /// accessors
    OccupancyMap& GetOccuMap() { return occu_map_; }
    LikelihoodField& GetLikelihood() { return field_; }

    std::vector<std::shared_ptr<Frame>>& GetFrames() { return frames_; }
    size_t NumFrames() const { return frames_.size(); }

    void SetId(size_t id) { id_ = id; }
    size_t GetId() const { return id_; }

    void SetPose(const SE2& pose);
    SE2 GetPose() const { return pose_; }

private:
    SE2 pose_;      // submap的pose, Tws，当前子地图在世界系下的位姿    
    size_t id_ = 0; // 子地图的id

    std::vector<std::shared_ptr<Frame>> frames_;  // 子地图submap中所有的关键帧扫描点云scan

    // 一个子地图对应一个似然场图像与一个栅格地图
    LikelihoodField field_;                       // 似然场，用于匹配
    OccupancyMap occu_map_;                       // 栅格地图
};

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_SUBMAP_H
