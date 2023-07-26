//
// Created by xiang on 2022/3/18.
//

#ifndef SLAM_IN_AUTO_DRIVING_LIKELIHOOD_FILED_H
#define SLAM_IN_AUTO_DRIVING_LIKELIHOOD_FILED_H

#include <opencv2/core.hpp>

#include "common/eigen_types.h"
#include "common/lidar_utils.h"

namespace sad {

/**
 * @description: 高斯似然场 Gaussian Likihood Field：将扫描数据与栅格数据进行匹配
 * 
 *          在ICP方法中，每迭代一次，必须将这些弹簧重新安装一遍，比较费时
 * 
 *          可以换一种思路，如果部位这些点云安装弹簧，而是认为这些点云是在空间中形成了一个磁场。
 *          磁场会吸引附近的点云，吸力随着距离平方衰减。实际上是一种似然场法的思想。
 * 
 *          可以在地图上每个点附近定义一个不断向外衰减的场。这里的场field与物理中的场不同，
 *          存在着一定的有效范围和分辨率，场可以随距离呈平方衰减，也可以是高斯衰减。
 *          当一个被测量点落在场附近时，可以用场的读数来作为该点的误差函数。
 * 
 * @return {*}
 */
class LikelihoodField {
public:
    /// 2D 场的模板，在设置target scan或map的时候生成
    struct ModelPoint {
        ModelPoint(int dx, int dy, float res) : dx_(dx), dy_(dy), residual_(res) {}
        int dx_ = 0;            // 偏移量
        int dy_ = 0;            // 偏移量
        float residual_ = 0;    // 距离平方
    };

    LikelihoodField() { BuildModel(); }

    /// 增加一个2D的目标scan
    void SetTargetScan(Scan2d::Ptr scan);

    /// 设置被配准的那个scan
    void SetSourceScan(Scan2d::Ptr scan);

    /// 从占据栅格地图生成一个似然场地图
    void SetFieldImageFromOccuMap(const cv::Mat& occu_map);

    /// 使用高斯牛顿法配准
    bool AlignGaussNewton(SE2& init_pose);

    /**
     * 使用g2o配准
     * @param init_pose 初始位姿 NOTE 使用submap时，给定相对于该submap的位姿，估计结果也是针对于这个submap的位姿
     * @return
     */
    bool AlignG2O(SE2& init_pose);

    /// 获取场函数，转换为RGB图像
    cv::Mat GetFieldImage();

    bool HasOutsidePoints() const { return has_outside_pts_; }

    void SetPose(const SE2& pose) { pose_ = pose; }

private:
    void BuildModel();

    SE2 pose_;  // T_W_S
    Scan2d::Ptr target_ = nullptr;
    Scan2d::Ptr source_ = nullptr;

    std::vector<ModelPoint> model_;  // 2D 模板
    cv::Mat field_;                  // 场函数
    bool has_outside_pts_ = false;   // 是否含有出了这个场的点

    // 参数配置
    inline static const float resolution_ = 20;  // 每米多少个像素
};

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_LIKELIHOOD_FILED_H
