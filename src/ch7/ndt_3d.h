//
// Created by xiang on 2022/7/14.
//

#ifndef SLAM_IN_AUTO_DRIVING_NDT_3D_H
#define SLAM_IN_AUTO_DRIVING_NDT_3D_H

#include "common/eigen_types.h"
#include "common/point_types.h"
#include <glog/logging.h>

namespace sad {

/**
 * 3D 形式的NDT
 * 
 * 大概思路：
 *  1. 将目标点云按一定分辨率分成若干体素；
 *  2. 计算每个体素中的点云高斯分布。第k个体素的均值和方差；
 *  3. 配准时，先计算每个点落在哪个体素中，然后建立该点与体素的均值和方差之间构成的残差；
 *  4. 利用高斯牛顿或者列文伯格马夸尔特LM算法进行位姿估计优化与更新。
 */
class Ndt3d {
public:
    enum class NearbyType {
        CENTER,   // 只考虑中心
        NEARBY6,  // 上下左右前后
    };

    struct Options {
        int max_iteration_ = 20;        // 最大迭代次数
        double voxel_size_ = 1.0;       // 体素大小，默认1m
        double inv_voxel_size_ = 1.0;   // 体素大小倒数
        int min_effective_pts_ = 10;    // 最近邻点数阈值
        int min_pts_in_voxel_ = 3;      // 每个栅格中最小点数
        double eps_ = 1e-2;             // 收敛判定条件
        double res_outlier_th_ = 20.0;  // 异常值拒绝阈值
        bool remove_centroid_ = false;  // 是否计算两个点云中心并移除中心？

        NearbyType nearby_type_ = NearbyType::NEARBY6;
    };

    using KeyType = Eigen::Matrix<int, 3, 1>;  // 体素的索引

    // 体素结构体数据
    struct VoxelData {
        VoxelData() {}
        VoxelData(size_t id) { idx_.emplace_back(id); }

        std::vector<size_t> idx_;      // 点云中点的索引
        Vec3d mu_ = Vec3d::Zero();     // 均值
        Mat3d sigma_ = Mat3d::Zero();  // 协方差
        Mat3d info_ = Mat3d::Zero();   // 协方差之逆，信息矩阵
    };

    // 构造函数
    Ndt3d() {
        options_.inv_voxel_size_ = 1.0 / options_.voxel_size_;
        GenerateNearbyGrids();
        ComputeNDTScorePara(); // 【新增】
    }

    // 构造函数
    Ndt3d(Options options) : options_(options) {
        options_.inv_voxel_size_ = 1.0 / options_.voxel_size_;
        GenerateNearbyGrids();
        ComputeNDTScorePara(); // 【新增】
    }

    /// 设置目标的Scan
    void SetTarget(CloudPtr target) {
        target_ = target;

        // 利用target构建NDT体素
        BuildVoxels();

        // 计算点云中心
        target_center_ = std::accumulate(target->points.begin(), target_->points.end(), Vec3d::Zero().eval(),
                                         [](const Vec3d& c, const PointType& pt) -> Vec3d { return c + ToVec3d(pt); }) / target_->size();
    }

    /// 设置被配准的Scan
    void SetSource(CloudPtr source) {
        source_ = source;

        source_center_ = std::accumulate(source_->points.begin(), source_->points.end(), Vec3d::Zero().eval(),
                                         [](const Vec3d& c, const PointType& pt) -> Vec3d { return c + ToVec3d(pt); }) / source_->size();
    }

    void SetGtPose(const SE3& gt_pose) {
        gt_pose_ = gt_pose;
        gt_set_ = true;
    }

    /// 使用gauss-newton方法进行ndt配准
    bool AlignNdt(SE3& init_pose);

    //【新增】设置内部NDT网格结构的体素分辨率
    void SetResolution(double resolution) {
        options_.voxel_size_ = resolution;
        options_.inv_voxel_size_ = 1.0 / options_.voxel_size_;
        ComputeNDTScorePara(); // 【新增】
    }

    // 【新增】预先计算NDT匹配得分计算相关参数
    void ComputeNDTScorePara() {
        gauss_c1_ = 10.0 * (1 - outlier_ratio_);
        gauss_c2_ = outlier_ratio_ / pow(options_.voxel_size_, 3);
        gauss_d3_ = -std::log(gauss_c2_);
        gauss_d1_ = -std::log(gauss_c1_ + gauss_c2_) - gauss_d3_;
        gauss_d2_ = -2 * std::log((-std::log(gauss_c1_ * std::exp(-0.5) + gauss_c2_) - gauss_d3_) / gauss_d1_);
    }

    // 【新增】获取NDT匹配得分
    double GetNdtMatchingScore() {
        return ndt_matching_score_;
    }

    // 【新增】为NDT设置新的体素栅格 
    void SetNdtVoxel(std::unordered_map<KeyType, VoxelData, hash_vec<3>> grids) {
        grids_.clear();
        grids_.swap(grids);
    }

private:
    void BuildVoxels();

    /// 根据最近邻的类型，生成附近网格
    void GenerateNearbyGrids();

    CloudPtr target_ = nullptr;
    CloudPtr source_ = nullptr;

    Vec3d target_center_ = Vec3d::Zero();
    Vec3d source_center_ = Vec3d::Zero();

    SE3 gt_pose_;
    bool gt_set_ = false;

    Options options_;

    std::unordered_map<KeyType, VoxelData, hash_vec<3>> grids_;  // 栅格数据
    std::vector<KeyType> nearby_grids_;                          // 附近的栅格

    // 【新增】
    // NDT匹配得分计算相关
    /** \brief The ratio of outliers of points w.r.t. a normal distribution,
     * Equation 6.7 [Magnusson 2009]. */
    double outlier_ratio_ = 0.55; // PCL中默认设为了0.55
    /** \brief The normalization constants used fit the point distribution to a
     * normal distribution, Equation 6.8 [Magnusson 2009]. */
    double gauss_d1_, gauss_d2_, gauss_d3_;
    double gauss_c1_, gauss_c2_;
    /** \brief The likelihood score of the transform applied to the input cloud,
     * Equation 6.9 and 6.10 [Magnusson 2009]. */
    
    double ndt_matching_score_;
};

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_NDT_3D_H
