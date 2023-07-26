//
// Created by xiang on 2022/3/18.
//

#include "ch6/g2o_types.h"
#include "ch6/likelihood_filed.h"

#include <glog/logging.h>

#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/robust_kernel.h>
#include <g2o/core/robust_kernel_impl.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/dense/linear_solver_dense.h>

namespace sad {

/**
 * @description: 上一帧的Scan设置为目标点云
 * @param {Ptr} scan
 * @return {*}
 */
void LikelihoodField::SetTargetScan(Scan2d::Ptr scan) {
    target_ = scan;

    // 在target目标点云上，生成似然场
    field_ = cv::Mat(1000, 1000, CV_32F, 30.0); // 1000*1000的图像，每个像素的值设为30.0

    // 遍历目标点云中的每个点
    for (size_t i = 0; i < scan->ranges.size(); ++i) {
        // 跳过无效点
        if (scan->ranges[i] < scan->range_min || scan->ranges[i] > scan->range_max)  continue;

        double real_angle = scan->angle_min + i * scan->angle_increment; // 当前激光点的角度

        // 极坐标转换为笛卡尔坐标，加上偏移量（图像尺寸的一半），得到激光点在似然场图像中的像素坐标（x,y）
        double x = scan->ranges[i] * std::cos(real_angle) * resolution_ + 500;  
        double y = scan->ranges[i] * std::sin(real_angle) * resolution_ + 500;

        // 在像素坐标附近附近填入似然场场函数，也就是为每个激光点贴上这个模板
        // 遍历模板中的每个点，共[2*20+1, 2*20+1]个模板点
        for (auto& model_pt : model_) {
            // 激光点的像素坐标，加上模板点的偏移量
            int xx = int(x + model_pt.dx_); 
            int yy = int(y + model_pt.dy_); 
            // 判断xx, yy坐标是否在似然场图像范围内，且该似然场的值是否大于模板点的距离值，
            if (xx >= 0 && xx < field_.cols && 
                yy >= 0 && yy < field_.rows && field_.at<float>(yy, xx) > model_pt.residual_) {
                // 若是，则用该模板的的距离值去更新似然场值
                field_.at<float>(yy, xx) = model_pt.residual_; // 似然场中的值，最大30.0，然后其它都是距离值
            }
        }
    }
}

/**
 * @description: 生成似然场模板
 * @return {*}
 */
void LikelihoodField::BuildModel() {
    const int range = 20;  // 生成多少个像素的模板
    // 遍历[-range, range] * [-range, range]，(2*range+1)*(2*range+1)个像素的模板
    for (int x = -range; x <= range; ++x) {
        for (int y = -range; y <= range; ++y) 
            model_.emplace_back(x, y, std::sqrt((x * x) + (y * y))); // 模板中每个点的坐标和距离
}

/**
 * @description: 当前帧scan设置为源始点云
 * @param {Ptr} scan
 * @return {*}
 */
void LikelihoodField::SetSourceScan(Scan2d::Ptr scan) { source_ = scan; }

/**
 * @description: 使用高斯牛顿法实现似然场来配准两个scan数据
 * @param {SE2&} init_pose
 * @return {*}
 */
bool LikelihoodField::AlignGaussNewton(SE2& init_pose) {
    int iterations = 10;
    double cost = 0, lastCost = 0;
    SE2 current_pose = init_pose;
    const int min_effect_pts = 20;  // 最小有效点数
    const int image_boarder = 20;   // 预留图像边界

    has_outside_pts_ = false;  // 是否含有出了这个场的点

    // 遍历十次
    for (int iter = 0; iter < iterations; ++iter) {
        Mat3d H = Mat3d::Zero();
        Vec3d b = Vec3d::Zero();
        cost = 0;

        int effective_num = 0;  // 有效点数

        // 遍历source
        for (size_t i = 0; i < source_->ranges.size(); ++i) {
            float r = source_->ranges[i];
            if (r < source_->range_min || r > source_->range_max) 
                continue;

            float angle = source_->angle_min + i * source_->angle_increment;
            if (angle < source_->angle_min + 30 * M_PI / 180.0 || angle > source_->angle_max - 30 * M_PI / 180.0) 
                continue;

            float theta = current_pose.so2().log();
            Vec2d pw = current_pose * Vec2d(r * std::cos(angle), r * std::sin(angle));

            // 获取这个2D点在似然场图像（1000，1000）中的图像坐标，乘以分辨率，加上图像中心偏移量
            Vec2i pf = (pw * resolution_ + Vec2d(500, 500)).cast<int>();

            // 判断这个点是否在似然场图像范围内
            if (pf[0] >= image_boarder && pf[0] < field_.cols - image_boarder && 
                pf[1] >= image_boarder && pf[1] < field_.rows - image_boarder) {

                effective_num++;

                // // 图像梯度
                // float dx = 0.5 * (field_.at<float>(pf[1], pf[0] + 1) - field_.at<float>(pf[1], pf[0] - 1));
                // float dy = 0.5 * (field_.at<float>(pf[1] + 1, pf[0]) - field_.at<float>(pf[1] - 1, pf[0]));

                // // 图像梯度（双线性插值 bilinear interpolation）
                // float dx = 0.5 * (math::GetPixelValue_bilinear<float>(field_, pf[0] + 1, pf[1]) - math::GetPixelValue_bilinear<float>(field_, pf[0] - 1, pf[1]));
                // float dy = 0.5 * (math::GetPixelValue_bilinear<float>(field_, pf[0], pf[1] + 1) - math::GetPixelValue_bilinear<float>(field_, pf[0], pf[1] - 1));

                // 图像梯度（双三次插值 bicubic interpolation）
                float dx = 0.5 * (math::GetPixelValue_bicubic<float>(field_, pf[0] + 1, pf[1]) - math::GetPixelValue_bicubic<float>(field_, pf[0] - 1, pf[1]));
                float dy = 0.5 * (math::GetPixelValue_bicubic<float>(field_, pf[0], pf[1] + 1) - math::GetPixelValue_bicubic<float>(field_, pf[0], pf[1] - 1));

                Vec3d J;
                J << resolution_ * dx, resolution_ * dy,
                    -resolution_ * dx * r * std::sin(angle + theta) + resolution_ * dy * r * std::cos(angle + theta);
                H += J * J.transpose();

                // float e = field_.at<float>(pf[1], pf[0]); 
                // float e = math::GetPixelValue_bilinear<float>(field_, pf[0], pf[1]); // 双线性插值
                float e = math::GetPixelValue_bicubic<float>(field_, pf[0], pf[1]); // 双三次插值
                b += -J * e;

                cost += e * e;
            } else {
                has_outside_pts_ = true;
            }
        }

        if (effective_num < min_effect_pts) {
            return false;
        }

        // solve for dx
        Vec3d dx = H.ldlt().solve(b);
        if (isnan(dx[0])) {
            break;
        }

        cost /= effective_num;
        if (iter > 0 && cost >= lastCost) {
            break;
        }

        LOG(INFO) << "iter " << iter << " cost = " << cost << ", effect num: " << effective_num;

        current_pose.translation() += dx.head<2>();
        current_pose.so2() = current_pose.so2() * SO2::exp(dx[2]);
        lastCost = cost;
    }

    init_pose = current_pose;
    return true;
}

/**
 * @description: 获取激光点的似然场图像，每个点周围绘制一个随距离衰减的圆圈，该圆是固定的，通过BuildModel()预先生成。
 * @return {*}
 */
cv::Mat LikelihoodField::GetFieldImage() {
    cv::Mat image(field_.rows, field_.cols, CV_8UC3);
    for (int x = 0; x < field_.cols; ++x) {
        for (int y = 0; y < field_.rows; ++y) {
            // 根据似然场的值，映射到0-255的灰度值
            float r = field_.at<float>(y, x) * 255.0 / 30.0; 
            // 激光点中心为蓝色，距离越远，颜色越黑
            image.at<cv::Vec3b>(y, x) = cv::Vec3b(uchar(r), uchar(r), uchar(r)); 
        }
    }
    return image;
}

/**
 * @description: g2o优化器版本基于似然场的扫描匹配算法
 * @param {SE2&} init_pose
 * @return {*}
 */
bool LikelihoodField::AlignG2O(SE2& init_pose) {
    using BlockSolverType = g2o::BlockSolver<g2o::BlockSolverTraits<3, 1>>;
    using LinearSolverType = g2o::LinearSolverCholmod<BlockSolverType::PoseMatrixType>;
    auto* solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);

    auto* v = new VertexSE2();  
    v->setId(0);
    v->setEstimate(init_pose);
    optimizer.addVertex(v);

    const double range_th = 15.0;  // 不考虑太远的scan，不准
    const double rk_delta = 0.8;

    has_outside_pts_ = false;
    // 遍历source
    for (size_t i = 0; i < source_->ranges.size(); ++i) {
        float r = source_->ranges[i];
        if (r < source_->range_min || r > source_->range_max) 
            continue;

        if (r > range_th) 
            continue;

        float angle = source_->angle_min + i * source_->angle_increment;
        if (angle < source_->angle_min + 30 * M_PI / 180.0 || angle > source_->angle_max - 30 * M_PI / 180.0) 
            continue;

        // 似然场约束边
        auto e = new EdgeSE2LikelihoodFiled(field_, r, angle, resolution_);
        e->setVertex(0, v);

        // 判断这个点是否在似然场图像范围内
        if (e->IsOutSide()) {
            has_outside_pts_ = true;
            delete e;
            continue;
        }

        e->setInformation(Eigen::Matrix<double, 1, 1>::Identity());
        auto rk = new g2o::RobustKernelHuber;
        rk->setDelta(rk_delta);
        e->setRobustKernel(rk);
        optimizer.addEdge(e);
    }

    optimizer.setVerbose(false);
    optimizer.initializeOptimization();
    optimizer.optimize(10);

    init_pose = v->estimate();
    return true;
}

void LikelihoodField::SetFieldImageFromOccuMap(const cv::Mat& occu_map) {
    const int boarder = 25;
    field_ = cv::Mat(1000, 1000, CV_32F, 30.0);

    for (int x = boarder; x < occu_map.cols - boarder; ++x) {
        for (int y = boarder; y < occu_map.rows - boarder; ++y) {
            if (occu_map.at<uchar>(y, x) < 127) {
                // 在该点生成一个model
                for (auto& model_pt : model_) {
                    int xx = int(x + model_pt.dx_);
                    int yy = int(y + model_pt.dy_);
                    if (xx >= 0 && xx < field_.cols && yy >= 0 && yy < field_.rows &&
                        field_.at<float>(yy, xx) > model_pt.residual_) {
                        field_.at<float>(yy, xx) = model_pt.residual_;
                    }
                }
            }
        }
    }
}

}  // namespace sad