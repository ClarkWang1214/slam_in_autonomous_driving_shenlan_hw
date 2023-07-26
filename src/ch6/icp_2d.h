//
// Created by xiang on 2022/3/15.
//

#ifndef SLAM_IN_AUTO_DRIVING_ICP_2D_H
#define SLAM_IN_AUTO_DRIVING_ICP_2D_H

#include "common/eigen_types.h"
#include "common/lidar_utils.h"

#include <pcl/search/kdtree.h>
#include <pcl/search/impl/kdtree.hpp>
#include <pcl/kdtree/kdtree_flann.h>

#include <pcl/segmentation/extract_clusters.h> // 【新增】

#include <g2o/core/base_unary_edge.h> // 【新增】
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/robust_kernel.h>
#include <g2o/core/robust_kernel_impl.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>

#include "ch6/g2o_types.h" // 【新增】

namespace sad {

using Point2d = pcl::PointXY;
using Point3d = pcl::PointXYZ;
using Cloud2d = pcl::PointCloud<Point2d>;
using Cloud3d = pcl::PointCloud<Point3d>;

/**
 * 第六章谈到的各种类型的ICP代码实现
 * 用法：先SetTarget, 此时构建target点云的KD树；再SetSource，然后调用Align*方法
 */
class Icp2d {
public:
    Icp2d() {}

    /// 设置目标的Scan，上一帧的Scan
    void SetTarget(Scan2d::Ptr target) {
        target_scan_ = target;
        BuildTargetKdTree();
    }

    /// 设置被配准的Scan
    void SetSource(Scan2d::Ptr source) { source_scan_ = source; }

    /// 使用高斯牛顿法进行配准
    bool AlignGaussNewton(SE2& init_pose);

    /// 使用高斯牛顿法进行配准, Point-to-Plane
    bool AlignGaussNewtonPoint2Plane(SE2& init_pose);

    /// 【新增】使用g2o进行点到点的ICP配准
    bool AlignG2OP2P(SE2& init_pose);   // 外部for循环迭代10次，内部for循环遍历scan点云，kdtree搜索current_pose转换的pw点得到qw，optimize(1)，每次迭代后更新current_pose
    bool AlignG2OP2P_2(SE2& init_pose); // for循环遍历scan点云，kdtree传入位姿边构造函数，g2o内部每次迭代更新的位姿转换得到pw点进行近邻搜索得到qw，optimize(10)
    bool AlignG2OP2P_3(SE2& init_pose); // for循环遍历scan点云，每个点都使用initial_pose转换的pw点，去进行近邻搜索得到qw，optimize(10)

    /// 【新增】使用g2o进行点到线的ICP配准（2D激光场景中，并不存在面，可以将点到线的形式看成低维的点到面）
    bool AlignG2OP2L(SE2& init_pose);
    bool AlignG2OP2L_2(SE2& init_pose);

    /**
     * @description: 
     * @return {*}
     */
    class EdgeSE2P2P_2 : public g2o::BaseUnaryEdge<2, Vec2d, VertexSE2> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        EdgeSE2P2P_2(const pcl::search::KdTree<Point2d>::Ptr kdtree, const Cloud2d::Ptr target_cloud, double range, double angle) : kdtree_(kdtree),  target_cloud_(target_cloud), range_(range), angle_(angle) {}

        // 判断当前激光点的最近邻集合是否为空，或者最小距离是否大于最大距离阈值
        bool isPointValid() { 
            // LOG(INFO) << "EdgeSE2P2P_2: isPointValid() ";
            auto* pose = dynamic_cast<const VertexSE2*>(_vertices[0]);
            theta_ = pose->estimate().so2().log(); // 当前位姿的角度
            // 世界系下点的坐标 p_i^W，极坐标转笛卡尔坐标公式
            pw_ = pose->estimate() * Vec2d(range_ * std::cos(angle_), range_ * std::sin(angle_));

            Point2d pt;
            pt.x = pw_.x();
            pt.y = pw_.y();

            // 在目标点云的KD树中查找一个最近邻，返回该最近邻的索引和距离
            kdtree_->nearestKSearch(pt, 1, nn_idx_, dis_);
            float max_dis2 = 0.01;
            // 判断最近邻集合是否非空，且最小距离是否小于最大距离阈值
            if (nn_idx_.size() > 0 && dis_[0] < max_dis2) {
                // 当前激光点在目标点云中的最近邻点坐标
                qw_ = Vec2d(target_cloud_->points[nn_idx_[0]].x, target_cloud_->points[nn_idx_[0]].y);   
                return true;
            }
            else 
                return false;
        }
        
        // 定义残差
        void computeError() override {
            LOG(INFO) << "EdgeSE2P2P_2: computeError() ";
            // 判断最近邻集合是否非空，且最小距离是否小于最大距离阈值
            if (isPointValid()) 
                _error =  pw_ - qw_; 
            else {
                _error = Vec2d(0, 0);
                setLevel(1);
            }
        }

        // 雅可比矩阵的解析形式
        void linearizeOplus() override {
            LOG(INFO) << "EdgeSE2P2P_2: linearizeOplus() ";
            if (isPointValid()) {
                _jacobianOplusXi <<  1, 0, 0, 1,  // de / dx， de / dy
                                    -range_ * std::sin(angle_ + theta_), range_ * std::cos(angle_ + theta_);  //  de / dtheta       
            } else {
                _jacobianOplusXi.setZero();
                setLevel(1);
            }                   
        }

        bool read(std::istream& is) override { return true; }
        bool write(std::ostream& os) const override { return true; }

    private:
        double range_ = 0;  // 距离
        double angle_ = 0;  // 角度

        // 【新增】
        double theta_ = 0;
        Vec2d pw_, qw_;
        const pcl::search::KdTree<Point2d>::Ptr kdtree_; 
        const Cloud2d::Ptr target_cloud_;
        std::vector<int> nn_idx_;    // 最近邻的索引
        std::vector<float> dis_;     // 最近邻的距离
    };

    /**
     * @description: 
     * @return {*}
     */
    class EdgeSE2P2L_2 : public g2o::BaseUnaryEdge<1, double, VertexSE2> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        EdgeSE2P2L_2(const pcl::search::KdTree<Point2d>::Ptr kdtree, const Cloud2d::Ptr target_cloud, double range, double angle) : kdtree_(kdtree),  target_cloud_(target_cloud), range_(range), angle_(angle) {}

        bool getIsLineFitSuccess() { return isLineFitSuccess_; }

        // 直线拟合是否成功
        bool isLineFitValid() { 
            auto* pose = dynamic_cast<const VertexSE2*>(_vertices[0]);
            theta_ = pose->estimate().so2().log(); // 当前位姿的角度
            // 世界系下点的坐标 p_i^W，极坐标转笛卡尔坐标公式
            pw_ = pose->estimate() * Vec2d(range_ * std::cos(angle_), range_ * std::sin(angle_));

            Point2d pt;
            pt.x = pw_.x();
            pt.y = pw_.y();

            // 在目标点云的KD树中查找一个最近邻，返回该最近邻的索引和距离
            kdtree_->nearestKSearch(pt, 5, nn_idx_, dis_);

            std::vector<Vec2d> effective_pts;  // 有效点
            float max_dis = 0.3;
            // 遍历所有五个近邻点
            for (int j = 0; j < nn_idx_.size(); ++j) {
                // 判断每个近邻点的距离是否处于最远阈值距离内
                if (dis_[j] < max_dis) 
                    // 若是，该近邻点符合要求，存储到向量中
                    effective_pts.emplace_back(Vec2d(target_cloud_->points[nn_idx_[j]].x, target_cloud_->points[nn_idx_[j]].y));
            }
            // 判断有效近邻点是否少于三个
            if (effective_pts.size() < 3) 
                // 若少于3个，则跳过当前激光点
                return false;

            
            // 利用当前点附近的几个有效近邻点，基于SVD奇异值分解，拟合出ax+by+c=0 中的最小直线系数 a,b,c，对应公式（6.11）
            if (math::FitLine2D(effective_pts, line_coeffs_)) {
                isLineFitSuccess_ = true;
                return isLineFitSuccess_;
            } else {
                isLineFitSuccess_ = false;
                return isLineFitSuccess_;
            }
        }
        
        // 定义残差
        void computeError() override {
            // 判断最近邻集合是否非空，且最小距离是否小于最大距离阈值
            if (isLineFitValid()) 
                _error[0] = line_coeffs_[0] * pw_[0] + line_coeffs_[1] * pw_[1] + line_coeffs_[2];
            else {
                _error[0] = 0.0;
                setLevel(1);
            }
        }

        // 雅可比矩阵的解析形式
        void linearizeOplus() override {
            if (isLineFitSuccess_) {
                _jacobianOplusXi << line_coeffs_[0], 
                                    line_coeffs_[1], 
                                    - line_coeffs_[0] * range_ * std::sin(angle_ + theta_) 
                                    + line_coeffs_[1] * range_ * std::cos(angle_ + theta_);        
            } else {
                _jacobianOplusXi.setZero();
                setLevel(1);
            }                   
        }

        bool read(std::istream& is) override { return true; }
        bool write(std::ostream& os) const override { return true; }

    private:
        double range_ = 0;  // 距离
        double angle_ = 0;  // 角度

        // 【新增】
        double theta_ = 0;
        Vec2d pw_, qw_;
        const pcl::search::KdTree<Point2d>::Ptr kdtree_; 
        const Cloud2d::Ptr target_cloud_;
        std::vector<int> nn_idx_;    // 最近邻的索引
        std::vector<float> dis_;     // 最近邻的距离

        Vec3d line_coeffs_;  // 拟合直线，组装J、H和误差

        bool isLineFitSuccess_ = false;
    };

private:
    // 建立目标点云的Kdtree
    void BuildTargetKdTree();

    // PCL版本的KD树实现2D点的最近邻搜索
    pcl::search::KdTree<Point2d>::Ptr kdtree_2d; 
    pcl::search::KdTree<Point3d>::Ptr kdtree_3d; 

    Cloud2d::Ptr target_cloud_2d;  // 2D目标点云
    Cloud3d::Ptr target_cloud_3d;  // 3D目标点云

    Scan2d::Ptr target_scan_ = nullptr;
    Scan2d::Ptr source_scan_ = nullptr;

    // 创建一个点云聚类对象
    pcl::EuclideanClusterExtraction<Point3d> clusterExtractor_3d;
};

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_ICP_2D_H
