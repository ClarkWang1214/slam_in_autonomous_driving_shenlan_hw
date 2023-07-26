//
// Created by xiang on 2022/7/14.
//

#include "ndt_3d.h"
#include "common/math_utils.h"

#include <glog/logging.h>
#include <Eigen/SVD>
#include <execution>

namespace sad {

/**
 * @description: 构建NDT体素
 * @return {*}
 */
void Ndt3d::BuildVoxels() {
    assert(target_ != nullptr); // 目标点云指针不能为空
    assert(target_->empty() == false);  // 目标点云不能为空
    grids_.clear(); // 清空体素栅格

    /// 分配体素索引
    std::vector<size_t> index(target_->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });

    // 生成体素栅格
    std::for_each(index.begin(), index.end(), [this](const size_t& idx) {
        auto pt = ToVec3d(target_->points[idx]);
        // 对目标点云中的每个点，计算它所在的体素栅格ID
        auto key = (pt * options_.inv_voxel_size_).cast<int>(); 
        // 查看该栅格是否已存在
        if (grids_.find(key) == grids_.end()) 
            grids_.insert({key, {idx}}); // 若不存在，则插入该栅格
        else 
            grids_[key].idx_.emplace_back(idx); // 若存在，则将该点的索引插入到该体素栅格中
    });

    /// 并行遍历所有体素栅格
    std::for_each(  std::execution::par_unseq, 
                    grids_.begin(), grids_.end(), 
                    [this](auto& v) { // lambda函数
                        // 判断每个体素栅格中是否多于3个点
                        if (v.second.idx_.size() > options_.min_pts_in_voxel_) { 
                            
                            // 计算每个体素栅格中点坐标的均值和协方差
                            math::ComputeMeanAndCov(v.second.idx_,                  // 当前栅格内所有点在原目标点云中的索引
                                                    v.second.mu_, v.second.sigma_,  // 返回的均值和协方差
                                                    [this](const size_t& idx) {     // lambda函数
                                                        // 获取目标点云中对应索引点的坐标，返回Eigen::Vector3d类型向量
                                                        return ToVec3d(target_->points[idx]); 
                                                    }); 
                            // 对协方差矩阵进行SVD分解
                            Eigen::JacobiSVD svd(v.second.sigma_, Eigen::ComputeFullU | Eigen::ComputeFullV);
                            Vec3d lambda = svd.singularValues(); // 获取奇异值
                            // 检查最大与最小奇异值，限制最小奇异值，避免出现奇异值过小的情况
                            if (lambda[1] < lambda[0] * 1e-3) 
                                lambda[1] = lambda[0] * 1e-3;
                            if (lambda[2] < lambda[0] * 1e-3) 
                                lambda[2] = lambda[0] * 1e-3;
                            // 计算奇异值对角矩阵的逆
                            Mat3d inv_lambda = Vec3d(1.0 / lambda[0], 1.0 / lambda[1], 1.0 / lambda[2]).asDiagonal();
                            // 恢复出协方差矩阵的逆.作为信息矩阵
                            v.second.info_ = svd.matrixV() * inv_lambda * svd.matrixU().transpose(); 
                        }
                    });

    // 遍历所有体素栅格
    for (auto iter = grids_.begin(); iter != grids_.end();) {
        // 判断栅格中点的数量是否大于阈值3个
        if (iter->second.idx_.size() > options_.min_pts_in_voxel_)
            iter++; // 若是，则继续遍历下一个栅格
        else
            iter = grids_.erase(iter);  // 若否，则删除点数不够3的的栅格
    }
}

/**
 * @description: 用NDT进行配准
 * @param {SE3&} init_pose 返回的初始位姿
 * @return {*}
 */
bool Ndt3d::AlignNdt(SE3& init_pose) {
    // LOG(INFO) << "aligning with ndt";
    assert(grids_.empty() == false);

    SE3 pose = init_pose;   // 初始位姿
    // 判断是否移除中心
    if (options_.remove_centroid_) {
        pose.translation() = target_center_ - source_center_;  // 设置平移初始值
        LOG(INFO) << "init trans set to " << pose.translation().transpose();
    }

    // 近邻栅格数：1或7
    int num_residual_per_point = 1;
    if (options_.nearby_type_ == NearbyType::NEARBY6)
        num_residual_per_point = 7;

    // 从待配准点云中随机采样，生成索引
    std::vector<int> index(source_->points.size()); 
    for (int i = 0; i < index.size(); ++i) 
        index[i] = i;

    // 我们来写一些并发代码
    int total_size = index.size() * num_residual_per_point;

    // 遍历20次
    for (int iter = 0; iter < options_.max_iteration_; ++iter) {
        std::vector<bool> effect_pts(total_size, false);                // 有效点
        std::vector<Eigen::Matrix<double, 3, 6>> jacobians(total_size); // 雅克比矩阵
        std::vector<Vec3d> errors(total_size);                          // 残差
        std::vector<Mat3d> infos(total_size);                           // 信息矩阵

        // gauss-newton 迭代
        std::for_each(  std::execution::par_unseq, // 并发
                        index.begin(), index.end(), 
                        [&](int idx) {
                            auto q = ToVec3d(source_->points[idx]); // 获取待配准点云中点的坐标
                            Vec3d qs = pose * q;  // 转换之后的q

                            // 计算qs所在的栅格以及它的最近邻栅格，体素大小为1.0
                            Vec3i key = (qs * options_.inv_voxel_size_).cast<int>();
                            // 遍历最近邻栅格，中心体素或者上下左右前后六个体素   
                            for (int i = 0; i < nearby_grids_.size(); ++i) {
                                auto key_off = key + nearby_grids_[i];  // 待配准点周围第i个近邻栅格
                            auto it = grids_.find(key_off);     // 在目标点云的体素栅格容器中查找待配准点近邻栅格所在的体素
                                // 第idx点周围第i个近邻栅格的索引
                                int real_idx = idx * num_residual_per_point + i;
                                // 判断是否找到
                                if (it != grids_.end()) {   
                                    // 若找到，
                                    auto& v = it->second;   // voxel，体素
                                    Vec3d e = qs - v.mu_;   // 残差，转换后的点坐标与体素均值作差，对应公式（7.12）

                                    // check chi2 th
                                    double res = e.transpose() * v.info_ * e; // 残差的转置乘以信息矩阵再乘以残差，得到卡方值
                                    // 卡方值大于阈值（20），则认为是异常值，将该点标记为无效点
                                    if (std::isnan(res) || res > options_.res_outlier_th_) {
                                        effect_pts[real_idx] = false;
                                        continue;
                                    }

                                    // build residual 构建残差，对应公式（7.16）
                                    Eigen::Matrix<double, 3, 6> J;
                                    J.block<3, 3>(0, 0) = -pose.so3().matrix() * SO3::hat(q);
                                    J.block<3, 3>(0, 3) = Mat3d::Identity();

                                    jacobians[real_idx] = J;
                                    errors[real_idx] = e;
                                    infos[real_idx] = v.info_;
                                    effect_pts[real_idx] = true;    // 有效点
                                } else 
                                    effect_pts[real_idx] = false;   // 无效点
                            }
                        });

        // 累加Hessian和error,计算dx
        // 原则上可以用reduce并发，写起来比较麻烦，这里写成accumulate
        double total_res = 0;
        int effective_num = 0;

        // 计算每一次迭代的NDT匹配得分累加和
        double score = 0;

        Mat6d H = Mat6d::Zero();
        Vec6d err = Vec6d::Zero();

        for (int idx = 0; idx < effect_pts.size(); ++idx) {
            if (!effect_pts[idx]) 
                continue;

            double res = errors[idx].transpose() * infos[idx] * errors[idx];
            total_res += res;

            score += std::exp(-res/2.0);
            
            // // 【新增】计算NDT匹配得分，借鉴PCL的NDT实现
            // // 对应公式（6.9）[Magnusson 2009]
            // // e^(-d_2/2 * (x_k - mu_k)^T Sigma_k^-1 (x_k - mu_k)) 
            // double e_x_cov_x = std::exp(-gauss_d2_ * res / 2);

            // LOG(INFO) << "gauss_d1_ = " << gauss_d1_;
            // LOG(INFO) << "gauss_d2_ = " << gauss_d2_;

            // // 匹配得分累加 对应公式（6.10）[Magnusson 2009]
            // score += - gauss_d1_ * e_x_cov_x;

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

        // 判断有效点数是否小于阈值10个
        if (effective_num < options_.min_effective_pts_) {
            LOG(WARNING) << "effective num too small: " << effective_num;
            return false;
        }

        Vec6d dx = H.inverse() * err; // 求解  dx = H^-1 * err
        // 令 x_k+1 = x_k + dx，更新位姿，再重新迭代，直到dx小于阈值
        pose.so3() = pose.so3() * SO3::exp(dx.head<3>());
        pose.translation() += dx.tail<3>();

        // 更新
        // LOG(INFO) << "iter " << iter << " total res: " << total_res << ", eff: " << effective_num
        //           << ", mean res: " << total_res / effective_num << ", dxn: " << dx.norm()
        //           << ", dx: " << dx.transpose();


        // 【新增】计算NDT匹配得分，借鉴PCL的NDT实现
        // ndt_matching_score_ = score / index.size();
        ndt_matching_score_ = score / index.size();

        // std::sort(chi2.begin(), chi2.end());
        // LOG(INFO) << "chi2 med: " << chi2[chi2.size() / 2] << ", .7: " << chi2[chi2.size() * 0.7]
        //           << ", .9: " << chi2[chi2.size() * 0.9] << ", max: " << chi2.back();

        if (gt_set_) {
            double pose_error = (gt_pose_.inverse() * pose).log().norm();
            LOG(INFO) << "iter " << iter << " pose error: " << pose_error;
        }

        // 当dx足够小时，收敛，高斯牛顿迭代终止
        if (dx.norm() < options_.eps_) {
            // LOG(INFO) << "converged, dx = " << dx.transpose();
            // LOG(INFO) << "NDT score: " << ndt_matching_score_;
            break;
        }
    }

    init_pose = pose;
    return true;
}

/**
 * @description: 生成附近邻域的体素栅格，未用到kdtree进行检索
 * @return {*}
 */
void Ndt3d::GenerateNearbyGrids() {
    if (options_.nearby_type_ == NearbyType::CENTER) { // 仅中心栅格
        nearby_grids_.emplace_back(KeyType::Zero()); 
    } else if (options_.nearby_type_ == NearbyType::NEARBY6) { // 上下左右前后六个栅格
        nearby_grids_ = {KeyType(0, 0, 0),  
                         KeyType(-1, 0, 0), KeyType(1, 0, 0), KeyType(0, 1, 0),
                         KeyType(0, -1, 0), KeyType(0, 0, -1), KeyType(0, 0, 1)};
    }
}

}  // namespace sad