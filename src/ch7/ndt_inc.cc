//
// Created by xiang on 2022/7/20.
//

#include "ch7/ndt_inc.h"
#include "common/math_utils.h"
#include "common/timer/timer.h"

#include <glog/logging.h>
#include <execution>
#include <set>

namespace sad {

/**
 * @description: 
 * 
 *      为每个体素设计一个缓冲区，当缓冲区中的点数大于一定数量时，才估计它的高斯参数。
 *      同时新加入的点也会加入这个缓冲区，达到一定数量后会重新估计
 * 
 * @param {CloudPtr} cloud_world
 * @return {*}
 */
void IncNdt3d::AddCloud(CloudPtr cloud_world) {
    // 每个体素都会有一个标记，标志它内部的高斯参数是否已被估计
    std::set<KeyType, less_vec<3>> active_voxels;  // 记录哪些voxel被更新
    for (const auto& p : cloud_world->points) {
        auto pt = ToVec3d(p);
        auto key = (pt * options_.inv_voxel_size_).cast<int>();
        auto iter = grids_.find(key);
        if (iter == grids_.end()) {
            // 栅格不存在
            // 新的点云加入缓冲区
            data_.push_front({key, {pt}});
            grids_.insert({key, data_.begin()});

            // 当缓冲区中的点数大于一定数量时，才会估计高斯参数
            if (data_.size() >= options_.capacity_) {
                // 删除一个尾部的数据
                grids_.erase(data_.back().first);
                data_.pop_back();
            }
        } else {
            // 栅格存在，添加点，更新缓存
            iter->second->second.AddPoint(pt);
            data_.splice(data_.begin(), data_, iter->second);  // 更新的那个放到最前
            iter->second = data_.begin();                      // grids时也指向最前
        }

        active_voxels.emplace(key);
    }

    // 更新active_voxels
    std::for_each(std::execution::par_unseq, 
                  active_voxels.begin(), active_voxels.end(),
                  [this](const auto& key) { 
                    UpdateVoxel(grids_[key]->second); 
                });
    flag_first_scan_ = false;
}

/**
 * @description: 获取附近栅格
 * @return {*}
 */
void IncNdt3d::GenerateNearbyGrids() {
    if (options_.nearby_type_ == NearbyType::CENTER) 
        nearby_grids_.emplace_back(KeyType::Zero());
    else if (options_.nearby_type_ == NearbyType::NEARBY6) { // 上下左右前后
        nearby_grids_ = {KeyType(0, 0, 0),  KeyType(-1, 0, 0), KeyType(1, 0, 0), KeyType(0, 1, 0),
                         KeyType(0, -1, 0), KeyType(0, 0, -1), KeyType(0, 0, 1)};
    }
}

/**
 * @description: 更新体素内部的高斯分布
 * @param {VoxelData&} v
 * @return {*}
 */
void IncNdt3d::UpdateVoxel(VoxelData& v) {
    if (flag_first_scan_) {
        if (v.pts_.size() > 1) {
            math::ComputeMeanAndCov(v.pts_, v.mu_, v.sigma_, [this](const Vec3d& p) { return p; });
            v.info_ = (v.sigma_ + Mat3d::Identity() * 1e-3).inverse();  // 避免出nan
        } else {
            v.mu_ = v.pts_[0];
            v.info_ = Mat3d::Identity() * 1e2;
        }

        v.ndt_estimated_ = true;
        v.pts_.clear();
        return;
    }

    if (v.ndt_estimated_ && v.num_pts_ > options_.max_pts_in_voxel_) 
        return;

    if (!v.ndt_estimated_ && v.pts_.size() > options_.min_pts_in_voxel_) {
        // 新增的voxel
        math::ComputeMeanAndCov(v.pts_, v.mu_, v.sigma_, [this](const Vec3d& p) { return p; });
        v.info_ = (v.sigma_ + Mat3d::Identity() * 1e-3).inverse();  // 避免出nan
        v.ndt_estimated_ = true;
        v.pts_.clear();
    } else if (v.ndt_estimated_ && v.pts_.size() > options_.min_pts_in_voxel_) {
        // 已经估计，而且还有新来的点
        Vec3d cur_mu, new_mu;
        Mat3d cur_var, new_var;
        math::ComputeMeanAndCov(v.pts_, cur_mu, cur_var, [this](const Vec3d& p) { return p; });
        math::UpdateMeanAndCov(v.num_pts_, v.pts_.size(), v.mu_, v.sigma_, cur_mu, cur_var, new_mu, new_var);

        v.mu_ = new_mu;
        v.sigma_ = new_var;
        v.num_pts_ += v.pts_.size();
        v.pts_.clear();

        // check info
        Eigen::JacobiSVD svd(v.sigma_, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Vec3d lambda = svd.singularValues();

        if (lambda[1] < lambda[0] * 1e-3) 
            lambda[1] = lambda[0] * 1e-3;
        
        if (lambda[2] < lambda[0] * 1e-3) 
            lambda[2] = lambda[0] * 1e-3;
        
        Mat3d inv_lambda = Vec3d(1.0 / lambda[0], 1.0 / lambda[1], 1.0 / lambda[2]).asDiagonal();
        v.info_ = svd.matrixV() * inv_lambda * svd.matrixU().transpose();
    }
}

/**
 * @description: 增量式NDT激光里程计估计位姿
 * @param {SE3&} init_pose
 * @return {*}
 */
bool IncNdt3d::AlignNdt(SE3& init_pose) {
    LOG(INFO) << "aligning with inc ndt, pts: " << source_->size() << ", grids: " << grids_.size();
    assert(grids_.empty() == false);

    SE3 pose = init_pose;

    // 对点的索引，预先生成
    int num_residual_per_point = 1;
    if (options_.nearby_type_ == NearbyType::NEARBY6) 
        num_residual_per_point = 7;
    
    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i) 
        index[i] = i;
    
    // 我们来写一些并发代码
    int total_size = index.size() * num_residual_per_point;

    for (int iter = 0; iter < options_.max_iteration_; ++iter) {
        std::vector<bool> effect_pts(total_size, false);
        std::vector<Eigen::Matrix<double, 3, 6>> jacobians(total_size);
        std::vector<Vec3d> errors(total_size);
        std::vector<Mat3d> infos(total_size);

        // gauss-newton 迭代
        // 最近邻，可以并发
        std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](int idx) {
            auto q = ToVec3d(source_->points[idx]);
            Vec3d qs = pose * q;  // 转换之后的q

            // 计算qs所在的栅格以及它的最近邻栅格
            Vec3i key = (qs * options_.inv_voxel_size_).cast<int>();

            for (int i = 0; i < nearby_grids_.size(); ++i) {
                Vec3i real_key = key + nearby_grids_[i];
                auto it = grids_.find(real_key);
                int real_idx = idx * num_residual_per_point + i;
                /// 这里要检查高斯分布是否已经估计
                if (it != grids_.end() && it->second->second.ndt_estimated_) {
                    auto& v = it->second->second;  // voxel
                    Vec3d e = qs - v.mu_;

                    // check chi2 th
                    double res = e.transpose() * v.info_ * e;
                    if (std::isnan(res) || res > options_.res_outlier_th_) {
                        effect_pts[real_idx] = false;
                        continue;
                    }

                    // build residual
                    Eigen::Matrix<double, 3, 6> J;
                    J.block<3, 3>(0, 0) = -pose.so3().matrix() * SO3::hat(q);   // 对R
                    J.block<3, 3>(0, 3) = Mat3d::Identity();                    // 对p

                    jacobians[real_idx] = J;
                    errors[real_idx] = e;
                    infos[real_idx] = v.info_;
                    effect_pts[real_idx] = true;
                } else 
                    effect_pts[real_idx] = false;
            }
        });

        // 累加Hessian和error,计算dx
        double total_res = 0;

        int effective_num = 0;

        Mat6d H = Mat6d::Zero();
        Vec6d err = Vec6d::Zero();

        for (int idx = 0; idx < effect_pts.size(); ++idx) {
            if (!effect_pts[idx]) 
                continue;

            // total_res += errors[idx].transpose() * infos[idx] * errors[idx];
            double e2 = errors[idx].transpose() * infos[idx] * errors[idx];
            total_res += e2;

            effective_num++;

            // TODO: Cauchy Kernel

            // 根据Huber定义的95 efficiency rule，由于NDT的残差是高斯分布，所以Cauchy核的控制参数delta取2.3849
            double delta = 2.3849; 

            double delta2 = delta * delta;
            double delta2_inv = 1.0 / delta2;
            double aux = delta2_inv * e2 + 1.0;

            Vec3d rho;
            rho[0] = delta2 * log(aux);         // Cauchy核函数
            rho[1] = 1.0 / aux;                     // Cauchy核函数的一阶导数
            rho[2] = -delta2_inv * pow(rho[1],2);   // Cauchy核函数的二阶导数

            // 加权信息矩阵，
            Mat3d weighted_infos = rho[1] * infos[idx]; // g2o中cauchy核函数的实现中，作者仅用了一阶导函数，其注释掉了二阶导函数
            // Mat3d weighted_infos = rho[1] * infos[idx] + 2 * rho[2] * errors[idx] * errors[idx].transpose(); // 二阶泰勒展开的加权信息矩阵

            H += jacobians[idx].transpose() * weighted_infos * jacobians[idx];      // 修改为加权信息矩阵
            err += -rho[1] * jacobians[idx].transpose() * infos[idx] * errors[idx]; // 乘以Cauchy核函数的一阶导

            // // 加权最小二乘的高斯牛顿解法，累加海森矩阵和残差，对应公式（7.15）
            // H += jacobians[idx].transpose() * infos[idx] * jacobians[idx];
            // err += -jacobians[idx].transpose() * infos[idx] * errors[idx];
        }

        if (effective_num < options_.min_effective_pts_) {
            LOG(WARNING) << "effective num too small: " << effective_num;
            init_pose = pose;
            return false;
        }

        Vec6d dx = H.inverse() * err;
        pose.so3() = pose.so3() * SO3::exp(dx.head<3>());
        pose.translation() += dx.tail<3>();

        // 更新
        LOG(INFO) << "iter " << iter << " total res: " << total_res << ", eff: " << effective_num
                  << ", mean res: " << total_res / effective_num << ", dxn: " << dx.norm()
                  << ", dx: " << dx.transpose();

        if (dx.norm() < options_.eps_) {
            LOG(INFO) << "converged, dx = " << dx.transpose();
            break;
        }
    }

    init_pose = pose;
    return true;
}

/**
 * @description: 计算残差和雅克比矩阵
 * @param {SE3&} input_pose
 * @param {Mat18d&} HTVH
 * @param {Vec18d&} HTVr
 * @return {*}
 */
void IncNdt3d::ComputeResidualAndJacobians(const SE3& input_pose, Mat18d& HTVH, Vec18d& HTVr) {
    assert(grids_.empty() == false);

    // 输入位姿，来自ESKF的Predict()函数预测得到的名义旋转R_、名义位移T_
    SE3 pose = input_pose;

    // 大部分流程和前面的AlignNdt()是一样的，只是会把z, H, R三者抛出去，而非自己处理
    int num_residual_per_point = 1;
    if (options_.nearby_type_ == NearbyType::NEARBY6) 
        num_residual_per_point = 7;

    // 初始化索引，0，1，2，3，4。。。
    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i) 
        index[i] = i;

    // 总的点云数目 乘以 近邻栅格数目（1或7）
    int total_size = index.size() * num_residual_per_point;
    std::vector<bool> effect_pts(total_size, false);                    // 用于标记有效点
    std::vector<Eigen::Matrix<double, 3, 18>> jacobians(total_size);    // 用于存储雅可比矩阵
    std::vector<Vec3d> errors(total_size);                              // 用于存储残差
    std::vector<Mat3d> infos(total_size);                               // 用于存储信息矩阵

    // gauss-newton 迭代
    // 最近邻，可以并发
    std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](int idx) {
        // 并发遍历到点云中的某个点，不是按顺序遍历的
        auto q = ToVec3d(source_->points[idx]);

        // 利用ESKF预测的名义位姿对该点进行转换
        // 雷达系转换到IMU系：P_I = R_IL * P_L + T_IL
        Vec3d qs = pose * q;  // R * q + t

        // 计算转换后点qs所在的栅格以及它的最近邻栅格
        Vec3i key = (qs * options_.inv_voxel_size_).cast<int>();

        // 遍历栅格，中心栅格 或者 六邻域栅格
        for (int i = 0; i < nearby_grids_.size(); ++i) {
            Vec3i real_key = key + nearby_grids_[i];
            auto it = grids_.find(real_key);
            int real_idx = idx * num_residual_per_point + i; // num_residual_per_point = 1或者7
            /// 这里要检查高斯分布是否已经估计
            if (it != grids_.end() && it->second->second.ndt_estimated_) {
                auto& v = it->second->second;  // 获取当前体素 voxel

                // NDT残差：R * q + t - mu
                Vec3d e = qs - v.mu_;

                
                // 误差的平方，info为信息矩阵，也就是所在栅格的协方差矩阵的逆 
                double res = e.transpose() * v.info_ * e;
                // check chi2 th 检查误差平方是否为nan值，或者是否大于阈值
                if (std::isnan(res) || res > options_.res_outlier_th_) {
                    effect_pts[real_idx] = false; // 若是，标记为无效点
                    continue; // 跳过当前点
                }

                // build residual 
                // 构建残差相对于旋转和平移的雅可比矩阵
                Eigen::Matrix<double, 3, 18> J; 
                J.setZero(); // 其它四项3x3的块矩阵均为零矩阵
                J.block<3, 3>(0, 0) = Mat3d::Identity();                   // NDT残差相对平移p的雅可比
                J.block<3, 3>(0, 6) = -pose.so3().matrix() * SO3::hat(q);  // NDT残差相对旋转R的雅可比

                jacobians[real_idx] = J;        // 雅可比矩阵
                errors[real_idx] = e;           // 误差
                infos[real_idx] = v.info_;      // 信息矩阵
                effect_pts[real_idx] = true;    // 标记为有效点
            } else 
                effect_pts[real_idx] = false;   // 标记为无效点
        }
    });

    // 累加Hessian和error,计算dx
    double total_res = 0;
    int effective_num = 0;

    HTVH.setZero();
    HTVr.setZero();

    // 每个点反馈的info信息矩阵因子
    // 由于NDT点数明显多于预测方程，可能导致估计结果向NDT倾斜，
    // 给信息矩阵添加一个乘积因子0.01，让更新部分更加平滑一些。
    const double info_ratio = 0.01;  

    for (int idx = 0; idx < effect_pts.size(); ++idx) {
        if (!effect_pts[idx]) 
            continue;

        total_res += errors[idx].transpose() * infos[idx] * errors[idx];
        effective_num++;

        HTVH += jacobians[idx].transpose() * infos[idx] * jacobians[idx] * info_ratio;
        HTVr += -jacobians[idx].transpose() * infos[idx] * errors[idx] * info_ratio;
    }

    LOG(INFO) << "effective: " << effective_num;
}

/**
 * @description: 构建NDT残差边
 * @param {VertexPose*} v
 * @return {*}
 */
void IncNdt3d::BuildNDTEdges(sad::VertexPose* v, std::vector<EdgeNDT*>& edges) {
    assert(grids_.empty() == false);
    SE3 pose = v->estimate();

    /// 整体流程和NDT一致，只是把查询函数放到edge内部，建立和v绑定的边
    for (const auto& pt : source_->points) {
        Vec3d q = ToVec3d(pt);
        auto edge = new EdgeNDT(v, q, 
                                [this](const Vec3d& qs, Vec3d& mu, Mat3d& info) -> bool {
                                    Vec3i key = (qs * options_.inv_voxel_size_).cast<int>();

                                    auto it = grids_.find(key);
                                    /// 这里要检查高斯分布是否已经估计
                                    if (it != grids_.end() && it->second->second.ndt_estimated_) {
                                        auto& v = it->second->second;  // voxel
                                        mu = v.mu_;
                                        info = v.info_;
                                        return true;
                                    } else 
                                        return false;
                                });

        if (edge->IsValid()) 
            edges.emplace_back(edge);
        else 
            delete edge;
    }
}

}  // namespace sad