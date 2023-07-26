//
// Created by xiang on 2022/7/7.
//

#include "icp_inc_3d.h"
#include "common/math_utils.h"

#include <execution>

namespace sad {

/**
 * @description: 计算残差和雅克比矩阵【新增】
 * @param {SE3&} input_pose
 * @param {Mat18d&} HTVH
 * @param {Vec18d&} HTVr
 * @return {*}
 */
void IncIcp3d::ComputeResidualAndJacobians_P2Plane(const SE3& input_pose, Mat18d& HTVH, Vec18d& HTVr) {
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


void IncIcp3d::BuildTargetKdTree() {
    kdtree_ = std::make_shared<KdTree>();
    kdtree_->BuildTree(target_);
    kdtree_->SetEnableANN();
}

}  // namespace sad