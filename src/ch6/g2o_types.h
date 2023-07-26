//
// Created by xiang on 2022/3/22.
//

#ifndef SLAM_IN_AUTO_DRIVING_G2O_TYPES_H
#define SLAM_IN_AUTO_DRIVING_G2O_TYPES_H

#include <g2o/core/base_binary_edge.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/base_vertex.h>

#include <glog/logging.h>
#include <opencv2/core.hpp>

#include "common/eigen_types.h"
#include "common/math_utils.h"

namespace sad {

class VertexSE2 : public g2o::BaseVertex<3, SE2> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    void setToOriginImpl() override { _estimate = SE2(); }

    /**
     * @description: 定义SE2位姿顶点的广义加法运算，用于位姿更新
     * @param {double*} update  x,y,theta
     * @return {*}
     */    
    void oplusImpl(const double* update) override {
        _estimate.translation()[0] += update[0];    // x,y坐标直接累加
        _estimate.translation()[1] += update[1];
        // 角度使用SO2的exp指数映射函数，将角度转换为2D旋转矩阵
        // R = [cos(theta) -sin(theta)
        //      sin(theta)  cos(theta)]
        _estimate.so2() = _estimate.so2() * SO2::exp(update[2]);    
    }

    bool read(std::istream& is) override { return true; }
    bool write(std::ostream& os) const override { return true; }
};

/**
 * @description: 似然场的一元边，定义了似然场的残差计算和雅可比矩阵计算
 * @return {*}
 */
class EdgeSE2LikelihoodFiled : public g2o::BaseUnaryEdge<1, double, VertexSE2> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EdgeSE2LikelihoodFiled(const cv::Mat& field_image, double range, double angle, float resolution = 10.0)
        : field_image_(field_image), range_(range), angle_(angle), resolution_(resolution) {}

     /**
     * @description: 判定当前激光点对应似然图像上的像素坐标是否越界
     * @return {*}
     */    
    bool IsOutSide() {
        VertexSE2* v = (VertexSE2*)_vertices[0];
        SE2 pose = v->estimate();
        Vec2d pw = pose * Vec2d(range_ * std::cos(angle_), range_ * std::sin(angle_));
        // 世界系点坐标到似然场图像坐标的转换关系 p_i^f，对应公式（6.20）
        Vec2i pf = (pw * resolution_ + Vec2d(field_image_.rows / 2, field_image_.cols / 2)).cast<int>();  // 图像坐标

        if (pf[0] >= image_boarder_ && pf[0] < field_image_.cols - image_boarder_ && pf[1] >= image_boarder_ &&
            pf[1] < field_image_.rows - image_boarder_) 
            return false;
        else 
            return true;
    }

    /**
     * @description: 计算残差
     * @return {*}
     */    
    void computeError() override {
        VertexSE2* v = (VertexSE2*)_vertices[0];
        SE2 pose = v->estimate();
        // 世界系下点的坐标 p_i^W，极坐标转笛卡尔坐标公式
        Vec2d pw = pose * Vec2d(range_ * std::cos(angle_), range_ * std::sin(angle_));
        // 世界系点坐标到似然场图像坐标的转换关系 p_i^f，对应公式（6.20）
        Vec2d pf = pw * resolution_ + Vec2d(field_image_.rows / 2, field_image_.cols / 2) - Vec2d(0.5, 0.5);  // 图像坐标

        if (pf[0] >= image_boarder_ && pf[0] < field_image_.cols - image_boarder_ && 
            pf[1] >= image_boarder_ && pf[1] < field_image_.rows - image_boarder_) {
            // 对似然场图像进行双线性插值，得到似然场的值
            _error[0] = math::GetPixelValue_bilinear<float>(field_image_, pf[0], pf[1]); // 双线性插值
            // _error[0] = math::GetPixelValue_bicubic<float>(field_image_, pf[0], pf[1]); // 双三次插值
            // _error[0] = field_image_.at<float>(pf[1], pf[0]);
        } else {
            _error[0] = 0;
            setLevel(1);
        }
    }

    /**
     * @description: 定义雅可比矩阵的解析形式
     * @return {*}
     */   
    void linearizeOplus() override {
        VertexSE2* v = (VertexSE2*)_vertices[0];
        SE2 pose = v->estimate();
        float theta = pose.so2().log();
        Vec2d pw = pose * Vec2d(range_ * std::cos(angle_), range_ * std::sin(angle_));
        Vec2d pf = pw * resolution_ + Vec2d(field_image_.rows / 2, field_image_.cols / 2) - Vec2d(0.5, 0.5);  // 图像坐标

        if (pf[0] >= image_boarder_ && pf[0] < field_image_.cols - image_boarder_ && pf[1] >= image_boarder_ &&
            pf[1] < field_image_.rows - image_boarder_) {
            // 图像梯度（双线性插值 bilinear interpolation）
            float dx = 0.5 * (math::GetPixelValue_bilinear<float>(field_image_, pf[0] + 1, pf[1]) -
                              math::GetPixelValue_bilinear<float>(field_image_, pf[0] - 1, pf[1]));
            float dy = 0.5 * (math::GetPixelValue_bilinear<float>(field_image_, pf[0], pf[1] + 1) -
                              math::GetPixelValue_bilinear<float>(field_image_, pf[0], pf[1] - 1));

            // // 图像梯度（双三次插值 bicubic interpolation）
            // float dx = 0.5 * (math::GetPixelValue_bicubic<float>(field_image_, pf[0] + 1, pf[1]) -
            //                   math::GetPixelValue_bicubic<float>(field_image_, pf[0] - 1, pf[1]));
            // float dy = 0.5 * (math::GetPixelValue_bicubic<float>(field_image_, pf[0], pf[1] + 1) -
            //                   math::GetPixelValue_bicubic<float>(field_image_, pf[0], pf[1] - 1));

            _jacobianOplusXi <<  resolution_ * dx, resolution_ * dy,
                                -resolution_ * dx * range_ * std::sin(angle_ + theta) +
                                 resolution_ * dy * range_ * std::cos(angle_ + theta);
        } else {
            _jacobianOplusXi.setZero();
            setLevel(1);
        }
    }

    bool read(std::istream& is) override { return true; }
    bool write(std::ostream& os) const override { return true; }

private:
    const cv::Mat& field_image_;
    double range_ = 0;
    double angle_ = 0;
    float resolution_ = 10.0;
    inline static const int image_boarder_ = 10;
};

/**
 * SE2 pose graph使用
 * error = v1.inv * v2 * meas.inv
 */
class EdgeSE2 : public g2o::BaseBinaryEdge<3, SE2, VertexSE2, VertexSE2> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EdgeSE2() {}

    /**
     * @description: 计算子地图的位姿图的残差
     * @return {*}
     */    
    void computeError() override {
        VertexSE2* v1 = (VertexSE2*)_vertices[0];
        VertexSE2* v2 = (VertexSE2*)_vertices[1];
        // SE2位姿顶点v1的逆乘以SE2位姿顶点v2乘以测量值的逆
        // _error = (v1->estimate().inverse() * v2->estimate() * _measurement.inverse()).log(); // (T_1^-1 * T_2) * T_12^-1
        // _error = (_measurement.inverse() * v1->estimate().inverse() * v2->estimate()).log(); // T_12^-1 * (T_1^-1 * T_2)

        // 【新增】
        const SO2 R1 = v1->estimate().so2(); 
        const SO2 R2 = v2->estimate().so2(); 
        const SO2 R12 = _measurement.so2();
        const Vec2d p1 = v1->estimate().translation();
        const Vec2d p2 = v2->estimate().translation();
        const Vec2d p12 = _measurement.translation();

        double delta_theta = R1.log() - R2.log() - R12.log();
        // Vec2d delta_p = R12.inverse() * (R1.inverse() * (p2 - p1) - p12);
        Vec2d delta_p = R1.inverse() * (p2 - p1) - p12;
        LOG(INFO) << "delta_theta: " << delta_theta << ", delta_p: " << delta_p.transpose();
        _error[0] = delta_p.x(); // 2x1
        _error[1] = delta_p.y();
        _error[2] = delta_theta; // 1x1
    }

    // TODO jacobian
    // 子地图的位姿图优化包含两种观测：相邻位姿图之间的相对位姿观测；回环检测计算出来的两个子地图之间的相对位姿关系。
    // computeError中定义好的残差项的雅可比矩阵计算比较繁琐，这里交给g2o自动求导来完成，也就是说，这里的linearizeOplus函数不需要实现。

    void linearizeOplus() {
        VertexSE2* v1 = (VertexSE2*)_vertices[0];
        VertexSE2* v2 = (VertexSE2*)_vertices[1];

        const SO2 R1 = v1->estimate().so2(); 
        const SO2 R12 = _measurement.so2();
        const Vec2d p1 = v1->estimate().translation();
        const Vec2d p2 = v2->estimate().translation();
        const Vec2d p12 = _measurement.translation();
        
        // R1的导数
        Mat2d R1Derivative;
        double theta1 = R1.log();
        R1Derivative << -sin(theta1), cos(theta1), -cos(theta1), -sin(theta1); // 2x2 这是列主序

        _jacobianOplusXj.setZero();
        // _jacobianOplusXi.block<2, 2>(0, 0) = -R12.inverse().matrix() * R1.inverse().matrix();
        _jacobianOplusXi.block<2, 2>(0, 0) = -R1.inverse().matrix();
        // _jacobianOplusXi.block<2, 1>(0, 2) = R12.inverse().matrix() * R1Derivative * (p2 - p1);
        _jacobianOplusXi.block<2, 1>(0, 2) = R1Derivative * (p2 - p1);
        _jacobianOplusXi.block<1, 3>(2, 0) = Vec3d(0, 0, -1);

        _jacobianOplusXj.setIdentity();
        // _jacobianOplusXj.block<2, 2>(0, 0) = R12.inverse().matrix() * R1.inverse().matrix();
        _jacobianOplusXj.block<2, 2>(0, 0) = R1.inverse().matrix();
    }

    bool read(std::istream& is) override { return true; }
    bool write(std::ostream& os) const override { return true; }

private:
};

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_G2O_TYPES_H
