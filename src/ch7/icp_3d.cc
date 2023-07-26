//
// Created by xiang on 2022/7/7.
//

#include "icp_3d.h"
#include "common/math_utils.h"

#include <execution>

namespace sad {

bool Icp3d::AlignP2P(SE3& init_pose) {
    LOG(INFO) << "aligning with point to point";
    assert(target_ != nullptr && source_ != nullptr);

    SE3 pose = init_pose;
    if (!options_.use_initial_translation_) {
        pose.translation() = target_center_ - source_center_;  // 设置平移初始值
    }

    // 对点的索引，预先生成
    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i) 
        index[i] = i;

    // 我们来写一些并发代码
    std::vector<bool> effect_pts(index.size(), false);
    std::vector<Eigen::Matrix<double, 3, 6>> jacobians(index.size());
    std::vector<Vec3d> errors(index.size());

    for (int iter = 0; iter < options_.max_iteration_; ++iter) {
        // gauss-newton 迭代
        // 最近邻，可以并发
        std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](int idx) {
            auto q = ToVec3d(source_->points[idx]);
            Vec3d qs = pose * q;  // 转换之后的q
            std::vector<int> nn;
            kdtree_->GetClosestPoint(ToPointType(qs), nn, 1);

            if (!nn.empty()) {
                Vec3d p = ToVec3d(target_->points[nn[0]]);
                double dis2 = (p - qs).squaredNorm();
                if (dis2 > options_.max_nn_distance_) {
                    // 点离的太远了不要
                    effect_pts[idx] = false;
                    return;
                }

                effect_pts[idx] = true;

                // build residual
                Vec3d e = p - qs;
                Eigen::Matrix<double, 3, 6> J;
                J.block<3, 3>(0, 0) = pose.so3().matrix() * SO3::hat(q);
                J.block<3, 3>(0, 3) = -Mat3d::Identity();

                jacobians[idx] = J;
                errors[idx] = e;
            } else {
                effect_pts[idx] = false;
            }
        });

        // 累加Hessian和error,计算dx
        // 原则上可以用reduce并发，写起来比较麻烦，这里写成accumulate
        double total_res = 0;
        int effective_num = 0;
        auto H_and_err = std::accumulate(
            index.begin(), index.end(), std::pair<Mat6d, Vec6d>(Mat6d::Zero(), Vec6d::Zero()),
            [&jacobians, &errors, &effect_pts, &total_res, &effective_num](const std::pair<Mat6d, Vec6d>& pre,
                                                                           int idx) -> std::pair<Mat6d, Vec6d> {
                if (!effect_pts[idx]) {
                    return pre;
                } else {                    
                    double e2 = errors[idx].dot(errors[idx]);
                    // total_res += errors[idx].dot(errors[idx]);
                    total_res += e2;

                    effective_num++;

                    // TODO: Cauchy Kernel
                    double delta =  1.0; 

                    double delta2 = delta * delta;
                    double delta2_inv = 1.0 / delta2;
                    double aux = delta2_inv * e2 + 1.0;

                    Vec3d rho;
                    rho[0] = delta2 * log(aux);         // Cauchy核函数
                    rho[1] = 1.0 / aux;                     // Cauchy核函数的一阶导数
                    rho[2] = -delta2_inv * pow(rho[1],2);   // Cauchy核函数的二阶导数

                    // Mat3d weighted_infos = rho[1] * Mat3d::Identity() + 2 * rho[2] * errors[idx] * errors[idx].transpose();
                    Mat3d weighted_infos = rho[1] * Mat3d::Identity();

                    // LOG(INFO) << "rho" << rho.transpose() << ", weighted_infos: " << weighted_infos;
                    // LOG(INFO) << "jacobians[idx].transpose() * jacobians[idx]" << jacobians[idx].transpose() * jacobians[idx];
                    // LOG(INFO) << "- jacobians[idx].transpose() * errors[idx]" << (- jacobians[idx].transpose() * errors[idx]).transpose();

                    // LOG(INFO) << "jacobians[idx].transpose() * weighted_infos * jacobians[idx]" << jacobians[idx].transpose() * weighted_infos * jacobians[idx];
                    // LOG(INFO) << "- rho[1] * jacobians[idx].transpose() * errors[idx]" << (- rho[1] * jacobians[idx].transpose() * errors[idx]).transpose();

                    // return std::pair<Mat6d, Vec6d>(pre.first + jacobians[idx].transpose() * jacobians[idx],
                    //                                pre.second - jacobians[idx].transpose() * errors[idx]);
                    return std::pair<Mat6d, Vec6d>(pre.first + jacobians[idx].transpose() * weighted_infos * jacobians[idx],
                                                   pre.second - rho[1] * jacobians[idx].transpose() * errors[idx]);
                }
            });

        if (effective_num < options_.min_effective_pts_) {
            LOG(WARNING) << "effective num too small: " << effective_num;
            return false;
        }

        Mat6d H = H_and_err.first;
        Vec6d err = H_and_err.second;

        Vec6d dx = H.inverse() * err;
        pose.so3() = pose.so3() * SO3::exp(dx.head<3>());
        pose.translation() += dx.tail<3>();

        // 更新
        LOG(INFO) << "iter " << iter << " total res: " << total_res << ", eff: " << effective_num
                  << ", mean res: " << total_res / effective_num << ", dxn: " << dx.norm();

        if (gt_set_) {
            double pose_error = (gt_pose_.inverse() * pose).log().norm();
            LOG(INFO) << "iter " << iter << " pose error: " << pose_error;
        }

        if (dx.norm() < options_.eps_) {
            LOG(INFO) << "converged, dx = " << dx.transpose();
            break;
        }
    }

    init_pose = pose;
    return true;
}

/**
 * @description: 点到面ICP
 * @param {SE3&} init_pose
 * @return {*}
 */
bool Icp3d::AlignP2Plane(SE3& init_pose) {
    LOG(INFO) << "aligning with point to plane";
    assert(target_ != nullptr && source_ != nullptr);
    // 整体流程与p2p一致，读者请关注变化部分

    SE3 pose = init_pose;
    if (!options_.use_initial_translation_) 
        pose.translation() = target_center_ - source_center_;  // 设置平移初始值：目标点云中心-源点云中心

    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i) 
        index[i] = i;

    std::vector<bool> effect_pts(index.size(), false);                  // 用于标记有效点   
    std::vector<Eigen::Matrix<double, 1, 6>> jacobians(index.size());   // 用于存储雅可比矩阵
    std::vector<double> errors(index.size());                           // 用于存储残差

    // gauss-newton 迭代20次
    for (int iter = 0; iter < options_.max_iteration_; ++iter) {
        // gauss-newton 迭代
        // 最近邻，可以并发
        std::for_each(  std::execution::par_unseq, 
                        index.begin(), index.end(), 
                        [&](int idx) {
                            auto q = ToVec3d(source_->points[idx]);
                            Vec3d qs = pose * q;  // 转换之后的q
                            std::vector<int> nn;

                            // kd树中查找转换后点的5个最近邻
                            kdtree_->GetClosestPoint(ToPointType(qs), nn, 5);  
                            // 判断查找到近邻点数是否多于3个，平面方程拟合，a*x+b*y+c*z+d=0，最少需要4个点才能拟合出平面系数
                            if (nn.size() > 3) {
                                std::vector<Vec3d> nn_eigen;
                                // 遍历近邻点集
                                for (int i = 0; i < nn.size(); ++i) 
                                    // 将近邻点转换为Vec3d类型存储
                                    nn_eigen.emplace_back(ToVec3d(target_->points[nn[i]]));
                                
                                Vec4d n;
                                // 对这几个近邻点执行平面拟合，平面系数a,b,c,d存储在四维向量n中
                                if (!math::FitPlane(nn_eigen, n)) {
                                    effect_pts[idx] = false; // 平面拟合失败，标记为无效点
                                    return;
                                }

                                // 计算点到平面的距离
                                double dis = n.head<3>().dot(qs) + n[3]; 
                                // 添加阈值检查判断拟合出的平面是否合理
                                if (fabs(dis) > options_.max_plane_distance_) {
                                    // 点离的太远了不要
                                    effect_pts[idx] = false;
                                    return;
                                }

                                effect_pts[idx] = true; // 标记为有效点

                                // 构建雅可比矩阵，对应公式（7.7）
                                Eigen::Matrix<double, 1, 6> J;
                                J.block<1, 3>(0, 0) = -n.head<3>().transpose() * pose.so3().matrix() * SO3::hat(q);
                                J.block<1, 3>(0, 3) = n.head<3>().transpose();

                                jacobians[idx] = J;
                                errors[idx] = dis;
                            } else {
                                effect_pts[idx] = false;
                            }
                        });

        // 累加Hessian和error,计算dx
        // 原则上可以用reduce并发，写起来比较麻烦，这里写成accumulate
        double total_res = 0;
        int effective_num = 0;
        auto H_and_err = std::accumulate(   index.begin(), index.end(), 
                                            std::pair<Mat6d, Vec6d>(Mat6d::Zero(), Vec6d::Zero()),
                                            [&jacobians, &errors, &effect_pts, &total_res, &effective_num]
                                            (const std::pair<Mat6d, Vec6d>& pre, int idx) -> std::pair<Mat6d, Vec6d> {
                                                if (!effect_pts[idx]) 
                                                    return pre;
                                                else {
                                                    double e2 = errors[idx] * errors[idx];
                                                    // total_res += errors[idx] * errors[idx];
                                                    total_res += e2;

                                                    effective_num++;

                                                    // TODO: Cauchy Kernel
                                                    // ICP的残差是非正态分布，需要估计残差的方差，然后对残差进行归一化
                                                    double delta =  1.0; 

                                                    double delta2 = delta * delta;
                                                    double delta2_inv = 1.0 / delta2;
                                                    double aux = delta2_inv * e2 + 1.0;

                                                    Vec3d rho;
                                                    rho[0] = delta2_inv / log(aux);         // Cauchy核函数    
                                                    rho[1] = 1.0 / aux;                     // Cauchy核函数的一阶导数
                                                    rho[2] = -delta2_inv * pow(rho[1],2);   // Cauchy核函数的二阶导数

                                                    // double weighted_infos = rho[1] + 2 * rho[2] * errors[idx] * errors[idx].transpose();
                                                    double weighted_infos = rho[1]; // g2o源码中的实现，注释掉了二阶导部分

                                                    // return std::pair<Mat6d, Vec6d>(pre.first + jacobians[idx].transpose() * jacobians[idx],
                                                    //                                pre.second - jacobians[idx].transpose() * errors[idx]);
                                                    return std::pair<Mat6d, Vec6d>(pre.first + jacobians[idx].transpose() * weighted_infos * jacobians[idx],
                                                                                   pre.second - rho[1] * jacobians[idx].transpose() * errors[idx]);
                                                }
                                            });

        if (effective_num < options_.min_effective_pts_) {
            LOG(WARNING) << "effective num too small: " << effective_num;
            return false;
        }

        Mat6d H = H_and_err.first;
        Vec6d err = H_and_err.second;

        Vec6d dx = H.inverse() * err;
        pose.so3() = pose.so3() * SO3::exp(dx.head<3>());
        pose.translation() += dx.tail<3>();

        // 更新
        LOG(INFO) << "iter " << iter << " total res: " << total_res << ", eff: " << effective_num
                  << ", mean res: " << total_res / effective_num << ", dxn: " << dx.norm();

        if (gt_set_) {
            double pose_error = (gt_pose_.inverse() * pose).log().norm();
            LOG(INFO) << "iter " << iter << " pose error: " << pose_error;
        }

        if (dx.norm() < options_.eps_) {
            LOG(INFO) << "converged, dx = " << dx.transpose();
            break;
        }
    }

    init_pose = pose;
    return true;
}

/**
 * @description: 计算残差和雅克比矩阵【新增】
 * @param {SE3&} input_pose
 * @param {Mat18d&} HTVH
 * @param {Vec18d&} HTVr
 * @return {*}
 */
void Icp3d::ComputeResidualAndJacobians_P2Plane(const SE3& input_pose, Mat18d& HTVH, Vec18d& HTVr) {
    LOG(INFO) << "aligning with point to plane";
    assert(target_ != nullptr && source_ != nullptr);

    // 大部分流程和前面的AlignP2Plane()是一样的，只是会把z, H, R三者抛出去，而非自己处理

    // 输入位姿，来自ESKF的Predict()函数预测得到的名义旋转R_、名义位移T_
    SE3 pose = input_pose;
    if (!options_.use_initial_translation_) 
        pose.translation() = target_center_ - source_center_;  // 设置平移初始值

    // 初始化索引，0，1，2，3，4。。。
    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i) 
        index[i] = i;

    std::vector<bool> effect_pts(index.size(), false);                    // 用于标记有效点
    std::vector<Eigen::Matrix<double, 1, 18>> jacobians(index.size());    // 用于存储雅可比矩阵
    std::vector<double> errors(index.size());                             // 用于存储残差
    // std::vector<Mat3d> infos(index.size());                            // 用于存储信息矩阵

    // gauss-newton 迭代
    // 最近邻，可以并发
    std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](int idx) {
        // 并发遍历到点云中的某个点，不是按顺序遍历的
        auto q = ToVec3d(source_->points[idx]);

        // 利用ESKF预测的名义位姿对该点进行转换
        // 雷达系转换到IMU系：P_I = R_IL * P_L + T_IL
        Vec3d qs = pose * q;  // R * q + t

        std::vector<int> nn;

        // kd树中查找转换后点的5个最近邻
        kdtree_->GetClosestPoint(ToPointType(qs), nn, 5);  

        // 判断查找到近邻点数是否多于3个，平面方程拟合，a*x+b*y+c*z+d=0，最少需要4个点才能拟合出平面系数
        if (nn.size() > 3) {
            std::vector<Vec3d> nn_eigen;
            // 遍历近邻点集
            for (int i = 0; i < nn.size(); ++i) 
                // 将近邻点转换为Vec3d类型存储
                nn_eigen.emplace_back(ToVec3d(target_->points[nn[i]]));
            
            Vec4d n;
            // 对这几个近邻点执行平面拟合，平面系数a,b,c,d存储在四维向量n中
            if (!math::FitPlane(nn_eigen, n)) {
                effect_pts[idx] = false; // 平面拟合失败，标记为无效点
                return;
            }

            // 计算点到平面的距离
            double dis = n.head<3>().dot(qs) + n[3]; 
            // 添加阈值检查判断拟合出的平面是否合理
            if (fabs(dis) > options_.max_plane_distance_) {
                // 点离的太远了不要
                effect_pts[idx] = false;
                return;
            }

            // 构建雅可比矩阵，对应公式（7.7）
            Eigen::Matrix<double, 1, 18> J;
            J.setZero(); // 其它四项1x3的块矩阵均为零矩阵
            J.block<1, 3>(0, 0) = n.head<3>().transpose();
            J.block<1, 3>(0, 3) = -n.head<3>().transpose() * pose.so3().matrix() * SO3::hat(q);

            jacobians[idx] = J;
            errors[idx] = dis;
            effect_pts[idx] = true; // 标记为有效点
        } else {
            effect_pts[idx] = false;
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

        total_res += errors[idx] * errors[idx];
        effective_num++;

        HTVH += jacobians[idx].transpose() * jacobians[idx] * info_ratio;    // 18x18
        HTVr += -jacobians[idx].transpose() * errors[idx] * info_ratio;      // 18x1
    }

    LOG(INFO) << "effective: " << effective_num;
}


void Icp3d::BuildTargetKdTree() {
    kdtree_ = std::make_shared<KdTree>();
    kdtree_->BuildTree(target_);
    kdtree_->SetEnableANN();
}

bool Icp3d::AlignP2Line(SE3& init_pose) {
    LOG(INFO) << "aligning with point to line";
    assert(target_ != nullptr && source_ != nullptr);
    // 点线与点面基本是完全一样的

    SE3 pose = init_pose;
    pose.translation() = target_center_ - source_center_;  // 设置平移初始值
    LOG(INFO) << "init trans set to " << pose.translation().transpose();

    std::vector<int> index(source_->points.size());
    for (int i = 0; i < index.size(); ++i) {
        index[i] = i;
    }

    std::vector<bool> effect_pts(index.size(), false);
    std::vector<Eigen::Matrix<double, 3, 6>> jacobians(index.size());
    std::vector<Vec3d> errors(index.size());

    for (int iter = 0; iter < options_.max_iteration_; ++iter) {
        // gauss-newton 迭代
        // 最近邻，可以并发
        std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](int idx) {
            auto q = ToVec3d(source_->points[idx]);
            Vec3d qs = pose * q;  // 转换之后的q
            std::vector<int> nn;
            kdtree_->GetClosestPoint(ToPointType(qs), nn, 5);  // 这里取5个最近邻
            if (nn.size() == 5) {
                // convert to eigen
                std::vector<Vec3d> nn_eigen;
                for (int i = 0; i < 5; ++i) {
                    nn_eigen.emplace_back(ToVec3d(target_->points[nn[i]]));
                }

                Vec3d d, p0;
                if (!math::FitLine(nn_eigen, p0, d, options_.max_line_distance_)) {
                    // 失败的不要
                    effect_pts[idx] = false;
                    return;
                }

                Vec3d err = SO3::hat(d) * (qs - p0);

                if (err.norm() > options_.max_line_distance_) {
                    // 点离的太远了不要
                    effect_pts[idx] = false;
                    return;
                }

                effect_pts[idx] = true;

                // build residual
                Eigen::Matrix<double, 3, 6> J;
                J.block<3, 3>(0, 0) = -SO3::hat(d) * pose.so3().matrix() * SO3::hat(q);
                J.block<3, 3>(0, 3) = SO3::hat(d);

                jacobians[idx] = J;
                errors[idx] = err;
            } else {
                effect_pts[idx] = false;
            }
        });

        // 累加Hessian和error,计算dx
        // 原则上可以用reduce并发，写起来比较麻烦，这里写成accumulate
        double total_res = 0;
        int effective_num = 0;
        auto H_and_err = std::accumulate(
            index.begin(), index.end(), std::pair<Mat6d, Vec6d>(Mat6d::Zero(), Vec6d::Zero()),
            [&jacobians, &errors, &effect_pts, &total_res, &effective_num](const std::pair<Mat6d, Vec6d>& pre,
                                                                           int idx) -> std::pair<Mat6d, Vec6d> {
                if (!effect_pts[idx]) {
                    return pre;
                } else {
                    double e2 =errors[idx].dot(errors[idx]);
                    // total_res += errors[idx] * errors[idx];
                    total_res += e2;

                    effective_num++;

                    // TODO: Cauchy Kernel
                    // ICP的残差是非正态分布，需要估计残差的方差，然后对残差进行归一化
                    double delta =  1.0; 

                    double delta2 = delta * delta;
                    double delta2_inv = 1.0 / delta2;
                    double aux = delta2_inv * e2 + 1.0;

                    Vec3d rho;
                    rho[0] = delta2_inv / log(aux);         // Cauchy核函数
                    rho[1] = 1.0 / aux;                     // Cauchy核函数的一阶导数
                    rho[2] = -delta2_inv * pow(rho[1],2);   // Cauchy核函数的二阶导数

                    // Mat3d weighted_infos = rho[1] * Mat3d::Identity() + 2 * rho[2] * errors[idx] * errors[idx].transpose();
                    Mat3d weighted_infos = rho[1] * Mat3d::Identity(); // g2o源码中的实现，注释掉了二阶导部分

                    // return std::pair<Mat6d, Vec6d>(pre.first + jacobians[idx].transpose() * jacobians[idx],
                    //                                pre.second - jacobians[idx].transpose() * errors[idx]);
                    return std::pair<Mat6d, Vec6d>(pre.first + jacobians[idx].transpose() * weighted_infos * jacobians[idx],
                                                   pre.second - rho[1] * jacobians[idx].transpose() * errors[idx]);
                }
            });

        if (effective_num < options_.min_effective_pts_) {
            LOG(WARNING) << "effective num too small: " << effective_num;
            return false;
        }

        Mat6d H = H_and_err.first;
        Vec6d err = H_and_err.second;

        Vec6d dx = H.inverse() * err;
        pose.so3() = pose.so3() * SO3::exp(dx.head<3>());
        pose.translation() += dx.tail<3>();

        if (gt_set_) {
            double pose_error = (gt_pose_.inverse() * pose).log().norm();
            LOG(INFO) << "iter " << iter << " pose error: " << pose_error;
        }

        // 更新
        LOG(INFO) << "iter " << iter << " total res: " << total_res << ", eff: " << effective_num
                  << ", mean res: " << total_res / effective_num << ", dxn: " << dx.norm();

        if (dx.norm() < options_.eps_) {
            LOG(INFO) << "converged, dx = " << dx.transpose();
            break;
        }
    }

    init_pose = pose;
    return true;
}

}  // namespace sad