//
// Created by xiang on 22-9-16.
//

#ifndef SLAM_IN_AUTO_DRIVING_IESKF_HPP
#define SLAM_IN_AUTO_DRIVING_IESKF_HPP

#include "common/eigen_types.h"
#include "common/imu.h"
#include "common/math_utils.h"
#include "common/nav_state.h"

namespace sad {

/**
 * 迭代版的ESKF，运动方程与第3章保持一致
 *
 * @tparam S
 */
template <typename S>
class IESKF {
   public:
    using SO3 = Sophus::SO3<S>;                     // 旋转变量类型
    using VecT = Eigen::Matrix<S, 3, 1>;            // 向量类型
    using Vec18T = Eigen::Matrix<S, 18, 1>;         // 18维向量类型
    using Mat3T = Eigen::Matrix<S, 3, 3>;           // 3x3矩阵类型
    using MotionNoiseT = Eigen::Matrix<S, 18, 18>;  // 运动噪声类型
    using OdomNoiseT = Eigen::Matrix<S, 3, 3>;      // 里程计噪声类型
    using GnssNoiseT = Eigen::Matrix<S, 6, 6>;      // GNSS噪声类型
    using Mat18T = Eigen::Matrix<S, 18, 18>;        // 18维方差类型
    using NavStateT = NavState<S>;                  // 整体名义状态变量类型

    struct Options {
        Options() = default;
        /// IEKF配置
        int num_iterations_ = 3;  // 迭代次数
        double quit_eps_ = 1e-3;  // 终止迭代的dx大小

        /// IMU 测量与零偏参数
        double imu_dt_ = 0.01;         // IMU测量间隔
        double gyro_var_ = 1e-5;       // 陀螺测量标准差
        double acce_var_ = 1e-2;       // 加计测量标准差
        double bias_gyro_var_ = 1e-6;  // 陀螺零偏游走标准差
        double bias_acce_var_ = 1e-4;  // 加计零偏游走标准差

        /// RTK 观测参数
        double gnss_pos_noise_ = 0.1;                   // GNSS位置噪声
        double gnss_height_noise_ = 0.1;                // GNSS高度噪声
        double gnss_ang_noise_ = 1.0 * math::kDEG2RAD;  // GNSS旋转噪声

        /// 其他配置
        bool update_bias_gyro_ = true;  // 是否更新bias
        bool update_bias_acce_ = true;  // 是否更新bias
    };

    /**
     * 初始零偏取零
     */
    IESKF(Options option = Options()) : options_(option) { BuildNoise(option); }

    /**
     * 由外部指定初始零偏
     * @param init_bg
     * @param init_ba
     * @param gravity
     */
    IESKF(Options options, const VecT& init_bg, const VecT& init_ba, const VecT& gravity = VecT(0, 0, -9.8))
        : options_(options) {
        BuildNoise(options);
        bg_ = init_bg;
        ba_ = init_ba;
        g_ = gravity;
    }

    /// 设置初始条件
    void SetInitialConditions(Options options, const VecT& init_bg, const VecT& init_ba,
                              const VecT& gravity = VecT(0, 0, -9.8)) {
        BuildNoise(options);
        options_ = options;
        bg_ = init_bg;
        ba_ = init_ba;
        g_ = gravity;

        cov_ = 1e-4 * Mat18T::Identity();
        cov_.template block<3, 3>(6, 6) = 0.1 * math::kDEG2RAD * Mat3T::Identity();
    }

    /// 使用IMU递推
    bool Predict(const IMU& imu);

    /**
     * NDT观测函数，输入一个SE3 Pose, 返回本书(8.10)中的几个项
     * HT V^{-1} H
     * H^T V{-1} r
     * 二者都可以用求和的形式来做
     */
    using CustomObsFunc = std::function<void(const SE3& input_pose, Eigen::Matrix<S, 18, 18>& HT_Vinv_H,
                                             Eigen::Matrix<S, 18, 1>& HT_Vinv_r)>;

    /// 使用自定义观测函数更新滤波器
    bool UpdateUsingCustomObserve(CustomObsFunc obs);

    /// accessors
    /// 全量状态
    NavStateT GetNominalState() const { return NavStateT(current_time_, R_, p_, v_, bg_, ba_); }

    /// SE3 状态
    SE3 GetNominalSE3() const { return SE3(R_, p_); }

    /**
     * @description: 设置状态
     * @param {NavStated&} x
     * @return {*}
     */    
    void SetX(const NavStated& x) {
        current_time_ = x.timestamp_;
        R_ = x.R_;
        p_ = x.p_;
        v_ = x.v_;
        bg_ = x.bg_;
        ba_ = x.ba_;
    }

    // 设置协方差矩阵
    void SetCov(const Mat18T& cov) { cov_ = cov; }

    // 获取重力矢量
    Vec3d GetGravity() const { return g_; }

private:
   /**
     * @description: 构建噪声协方差矩阵 Q
     * @param {Options&} options
     * @return {*}
     */ 
    void BuildNoise(const Options& options) {
        // 噪声项不参与递推，单独归入噪声部分中。
        // 连续时间的噪声项可以看作随机过程的能力谱密度
        // 离散时间的噪声变量就是日常看到的随机变量
        double ev = options.acce_var_;      // 加计噪声
        double et = options.gyro_var_;      // 陀螺噪声
        double eg = options.bias_gyro_var_; // 陀螺零偏
        double ea = options.bias_acce_var_; // 加计零偏

        double ev2 = ev;  // * ev;
        double et2 = et;  // * et;
        double eg2 = eg;  // * eg;
        double ea2 = ea;  // * ea;

        // 设置协方差矩阵Q的对角线上的元素，对应公式（3.45）
        Q_.diagonal() << 0, 0, 0, 
                         ev2, ev2, ev2,     // 加计噪声
                         et2, et2, et2,     // 陀螺噪声
                         eg2, eg2, eg2,     // 陀螺零偏
                         ea2, ea2, ea2,     // 加计零偏
                         0, 0, 0;

        double gp2 = options.gnss_pos_noise_ * options.gnss_pos_noise_;
        double gh2 = options.gnss_height_noise_ * options.gnss_height_noise_;
        double ga2 = options.gnss_ang_noise_ * options.gnss_ang_noise_;
        gnss_noise_.diagonal() << gp2, gp2, gh2, ga2, ga2, ga2;
    }

    /**
     * @description: 更新名义状态变量，重置误差状态量error state
     * @return {*}
     */  
    void Update() {
        // 名义状态加上误差状态
        p_ += dx_.template block<3, 1>(0, 0);
        v_ += dx_.template block<3, 1>(3, 0);
        R_ = R_ * SO3::exp(dx_.template block<3, 1>(6, 0));

        if (options_.update_bias_gyro_) 
            bg_ += dx_.template block<3, 1>(9, 0);

        if (options_.update_bias_acce_) 
            ba_ += dx_.template block<3, 1>(12, 0);

        g_ += dx_.template block<3, 1>(15, 0);  // 更新重力矢量
    }

    double current_time_ = 0.0;

    // nominal state
    SO3 R_;
    VecT p_ = VecT::Zero();
    VecT v_ = VecT::Zero();
    VecT bg_ = VecT::Zero();
    VecT ba_ = VecT::Zero();
    VecT g_{0, 0, -9.8};

    // error state
    Vec18T dx_ = Vec18T::Zero();

    // covariance
    Mat18T cov_ = Mat18T::Identity();

    // noise
    MotionNoiseT Q_ = MotionNoiseT::Zero();
    GnssNoiseT gnss_noise_ = GnssNoiseT::Zero();

    Options options_;
};

using IESKFD = IESKF<double>;
using IESKFF = IESKF<float>;

template <typename S>
bool IESKF<S>::Predict(const IMU& imu) {
    /// Predict 部分与ESKF完全一样，不再解释
    assert(imu.timestamp_ >= current_time_);

    double dt = imu.timestamp_ - current_time_;
    if (dt > (5 * options_.imu_dt_) || dt < 0) {
        LOG(INFO) << "skip this imu because dt_ = " << dt;
        current_time_ = imu.timestamp_;
        return false;
    }

    VecT new_p = p_ + v_ * dt + 0.5 * (R_ * (imu.acce_ - ba_)) * dt * dt + 0.5 * g_ * dt * dt;
    VecT new_v = v_ + R_ * (imu.acce_ - ba_) * dt + g_ * dt;
    SO3 new_R = R_ * SO3::exp((imu.gyro_ - bg_) * dt);

    R_ = new_R;
    v_ = new_v;
    p_ = new_p;

    // 协方差矩阵F由误差运动学的线性形式给出，对应公式（8.4）
    Mat18T F = Mat18T::Identity();                                                  // 18x18的单位矩阵
    F.template block<3, 3>(0, 3) = Mat3T::Identity() * dt;                          // 公式（8.2a）相对于速度v的雅克比
    F.template block<3, 3>(3, 6) = -R_.matrix() * SO3::hat(imu.acce_ - ba_) * dt;   // 公式（8.2b）相对于旋转R的雅克比
    F.template block<3, 3>(3, 12) = -R_.matrix() * dt;                              // 公式（8.2b）相对于加计零偏ba的雅克比
    F.template block<3, 3>(3, 15) = Mat3T::Identity() * dt;                         // 公式（8.2b）相对于重力g的雅克比
    F.template block<3, 3>(6, 6) = SO3::exp(-(imu.gyro_ - bg_) * dt).matrix();      // 公式（8.2c）相对于旋转R的雅克比
    F.template block<3, 3>(6, 9) = -Mat3T::Identity() * dt;                         // 公式（8.2c）相对于陀螺零偏bg的雅克比

    // 更新协方差矩阵，对应公式（8.3）
    cov_ = F * cov_ * F.transpose() + Q_;
    current_time_ = imu.timestamp_;
    return true;
}

/**
 * @description: 使用自定义观测函数更新滤波器
 * @return {*}
 */
template <typename S>
bool IESKF<S>::UpdateUsingCustomObserve(IESKF::CustomObsFunc obs) {
    // H阵由用户给定

    SO3 start_R = R_;
    Eigen::Matrix<S, 18, 1> HTVr;
    Eigen::Matrix<S, 18, 18> HTVH;
    Eigen::Matrix<S, 18, Eigen::Dynamic> K;
    Mat18T Pk, Qk;

    for (int iter = 0; iter < options_.num_iterations_; ++iter) {
        // 调用obs function
        // 让外部算法计算对应的H^T V^{-1} H 和 H^T V^{-1} r，然后代入IESKF中得到增量dx
        obs(GetNominalSE3(), // 获取预测predict()得到名义旋转R_、名义位移T_
            HTVH, HTVr); 

        // 投影P
        Mat18T J = Mat18T::Identity();
        //  与公式（8.5）不一致？？  start_R.inverse() * R_ ？
        // J.template block<3, 3>(6, 6) = Mat3T::Identity() - 0.5 * SO3::hat((R_.inverse() * start_R).log()); 
        J.template block<3, 3>(6, 6) = Mat3T::Identity() - 0.5 * SO3::hat((start_R.inverse() * R_).log()); 

        // 在每次迭代，未收敛前，cov_不变
        // 使用的都是Predict()函数预测得到的 协方差矩阵 P_{k+1}，对应公式（8.3）
        // P_pred = F * Pk * F.transpose() + Qk;
        Pk = J * cov_ * J.transpose();

        // 卡尔曼更新
        Qk = (Pk.inverse() + HTVH).inverse();   // 这个记作中间变量，最后更新时可以用
        dx_ = Qk * HTVr;                        // 对应公式（8.22）
        // LOG(INFO) << "iter " << iter << " dx = " << dx_.transpose() << ", dxn: " << dx_.norm();

        // 误差状态dx 合入 名义变量中
        Update(); // R_ p_ v_ bg_ ba_ g_ 都更新了

        // 判断18维的误差状态（属于 \mathbb R^18 向量空间）的模长是否小于阈值
        if (dx_.norm() < options_.quit_eps_) 
            break; // 若误差状态迭代到足够小了，说明收敛了，退出迭代
    }

    // 更新协方差矩阵P_{k+1}，使用最后一次迭代的Qk，Pk，对应公式（8.23）
    cov_ = (Mat18T::Identity() - Qk * HTVH) * Pk;

    // IESKF结束迭代后，将此时的协方差矩阵 P_{k+1} 投影至结束时刻的切空间中，保持整个IESKF的一致性
    Mat18T J = Mat18T::Identity();
    Vec3d dtheta = (start_R.inverse() * R_).log();
    J.template block<3, 3>(6, 6) = Mat3T::Identity() - 0.5 * SO3::hat(dtheta);
    cov_ = J * cov_ * J.inverse();

    dx_.setZero();
    return true;
}

}  // namespace sad
#endif  // SLAM_IN_AUTO_DRIVING_IEKF_HPP
