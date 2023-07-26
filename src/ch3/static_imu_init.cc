//
// Created by xiang on 2021/11/11.
//

#include "ch3/static_imu_init.h"
#include "common/math_utils.h"

#include <glog/logging.h>

namespace sad {

/**
 * @description: 为了让ESKF跑起来，需要知道一些初始条件（初始零偏、初始重力方向）
 *               静止初始化，将IMU放在某个地方静止一段时间，由于物体本身没有任何运动，可以简单认为IMU的陀螺仪只测到了零偏，而加速度计测到了零偏与重力的和
 *               可以设置一个静止初始化流程来获取这些量：
 *                  1. 将IMU静止一段时间（10秒），静止检查由轮速计判定，当两轮的轮速均小于阈值时，认为车辆静止。在没有轮速测量的场合，也可以直接认为车辆静止
 *                     来测定相关变量。
 *                  2. 统计静止时间内的陀螺仪与加速度计的均值 \bar{d_g}, \bar{d_a}
 *                  3. 由于车辆并未发生转动，这段时间的陀螺均值可以取 b_g = \bar{d_g}
 *                  4. 加速度计的测量方程为   \tilde{a} = R^T * (a - g) + b_a + n_a
 *                     当车辆静止时，实际加速度a = 0，旋转R看成是单位阵I，加速度计实际测到 b_a - g，其中零偏b_a为小量，g的长度为固定值。在这些情况下，我们取
 *                     方向为 -\bar{d_a}，大小为9.8的矢量作为重力矢量。这一步确定重力的朝向。
 *                  5. 将这段时间的加速度计读数去掉重力，重新计算加速度计的均值 \bar{d_a}
 *                  6. 取加速度零偏 b_a = \bar{d_a}
 *                  7. 同时，认为零偏不动，估计陀螺与加速度计的测量方差。该方差可以用于ESKF的噪声参数。
 * @param {IMU&} imu
 * @return {*}
 */
bool StaticIMUInit::AddIMU(const IMU& imu) {
    // 判断是否初始化成功
    if (init_success_) 
        return true;

    // 判断是否使用odom来判断车辆静止 且 当前车辆非静止
    if (options_.use_speed_for_static_checking_ && !is_static_) {
        LOG(WARNING) << "等待车辆静止";      // 提示等待车辆静止
        init_imu_deque_.clear();            // IMU队列清空
        return false;
    }

    // 判断IMU队列是否为空
    if (init_imu_deque_.empty()) {
        // 记录初始静止时间
        init_start_time_ = imu.timestamp_; // IMU时间戳
    }

    // 记入初始化队列
    init_imu_deque_.push_back(imu);

    double init_time = imu.timestamp_ - init_start_time_;  // 初始化经过时间
    // 判断初始化经过时间是否大于10秒阈值
    if (init_time > options_.init_time_seconds_) 
        // 若大于10秒阈值，则尝试初始化逻辑，静止不够10秒，则不进行初始化
        TryInit();

    // 判断IMU初始化队列长度是否大于最大长度阈值2000
    while (init_imu_deque_.size() > options_.init_imu_queue_max_size_) 
        // 若超过阈值，弹出最旧的IMU数据，维持初始化队列长度
        init_imu_deque_.pop_front();

    current_time_ = imu.timestamp_; // 更新当前时间
    return false;
}

/**
 * @description: 由轮速计判断车辆是否静止
 * @param {Odom&} odom
 * @return {*}
 */
bool StaticIMUInit::AddOdom(const Odom& odom) {
    // 判断车辆是否静止
    if (init_success_) 
        return true;

    // 判断【左轮脉冲数】和【右轮脉冲数】是否都小于【静止情况下的里程计脉冲阈值（默认：5）】
    if (odom.left_pulse_ < options_.static_odom_pulse_ && odom.right_pulse_ < options_.static_odom_pulse_) {
        is_static_ = true;  // 若都小于阈值，则认为车辆静止
    } else {
        is_static_ = false; // 否则，认为车辆非静止
    }

    // 更新当前时间为接收到的轮速计的时间戳
    current_time_ = odom.timestamp_;
    return true;
}

/**
 * @description: 此时，IMU静止时间已超10秒阈值，可对系统进行初始化
 * @return {*}
 */
bool StaticIMUInit::TryInit() {
    // 判定IMU初始化队列长度是否小于10
    if (init_imu_deque_.size() < 10) 
        // 说明IMU数据太少，初始化失败
        return false;

    // 计算均值和方差
    Vec3d mean_gyro, mean_acce; // 陀螺仪均值，加速度计均值
    math::ComputeMeanAndCovDiag(init_imu_deque_, mean_gyro, cov_gyro_, 
                                [](const IMU& imu) { return imu.gyro_; }); // lambda表达式，返回IMU的陀螺仪读数
    math::ComputeMeanAndCovDiag(init_imu_deque_, mean_acce, cov_acce_, 
                                [this](const IMU& imu) { return imu.acce_; });  // lambda表达式，返回IMU的加速度计读数

    // 以acce均值为方向，取9.8长度为重力
    LOG(INFO) << "mean acce: " << mean_acce.transpose();

    // 取上面估计的（归一化后的加速度读数均值的负值）作为方向，乘以重力大小9.81，得到重力矢量
    gravity_ = -mean_acce / mean_acce.norm() * options_.gravity_norm_;

    // 上面确定了初始重力矢量后，加速度计读数减去该重力矢量，重新计算加速度计的均值和协方差
    math::ComputeMeanAndCovDiag(init_imu_deque_, mean_acce, cov_acce_,
                                [this](const IMU& imu) { return imu.acce_ + gravity_; });

    // 检查估计的陀螺仪方差是否大于最大静止陀螺方差
    if (cov_gyro_.norm() > options_.max_static_gyro_var) {
        // 若大于最大静止陀螺方差，估计的测量方差太大了，则初始化失败
        LOG(ERROR) << "陀螺仪测量噪声太大" << cov_gyro_.norm() << " > " << options_.max_static_gyro_var;
        return false;
    }

    // 检查估计的加计协方差是否大于最大静止加计方差
    if (cov_acce_.norm() > options_.max_static_acce_var) {
        // 若大于最大静止加计方差，估计的测量方差太大了，则初始化失败
        LOG(ERROR) << "加计测量噪声太大" << cov_acce_.norm() << " > " << options_.max_static_acce_var;
        return false;
    }

    // 估计测量噪声和零偏
    init_bg_ = mean_gyro;   // 初始陀螺零偏取这段时间的陀螺均值
    init_ba_ = mean_acce;   // 初始加计零偏取这段时间的加计均值（去掉重力后）

    LOG(INFO) << "IMU 初始化成功，初始化时间= " << current_time_ - init_start_time_ << ", bg = " << init_bg_.transpose()
              << ", ba = " << init_ba_.transpose() << ", gyro sq = " << cov_gyro_.transpose()
              << ", acce sq = " << cov_acce_.transpose() << ", grav = " << gravity_.transpose()
              << ", norm: " << gravity_.norm();
    LOG(INFO) << "mean gyro: " << mean_gyro.transpose() << " acce: " << mean_acce.transpose();
    init_success_ = true;   // 初始化成功
    return true;
}

}  // namespace sad
