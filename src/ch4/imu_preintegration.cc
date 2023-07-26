#include "imu_preintegration.h"
#include <glog/logging.h>

namespace sad {

/**
 * @description: 构造函数
 * @param {Options} options
 * @return {*}
 */
IMUPreintegration::IMUPreintegration(Options options) {
    bg_ = options.init_bg_;
    ba_ = options.init_ba_;
    const float ng2 = options.noise_gyro_ * options.noise_gyro_;
    const float na2 = options.noise_acce_ * options.noise_acce_;
    noise_gyro_acce_.diagonal() << ng2, ng2, ng2, na2, na2, na2;
}

/**
 * @description: 预积分更新过程
 *                  1. 更新位置和速度的累加测量值dp_{i,j}，dv_{i,j}
 *                  2. 更新运动模型的噪声系数矩阵A和B，用的是上一时刻的累乘量dR_{i,j-1}
 *                  3. 更新观测量pvq对两个零偏ba，bg的五个雅可比矩阵
 *                  4. 更新旋转的累乘测量值dR_{i,j}
 *                  5. 更新积分累计时间dt
 * @param {IMU} &imu
 * @param {double} dt
 * @return {*}
 */
void IMUPreintegration::Integrate(const IMU &imu, double dt) {
    // 去掉零偏的测量
    Vec3d gyr = imu.gyro_ - bg_;  // 陀螺
    Vec3d acc = imu.acce_ - ba_;  // 加计

    dp_ = dp_ + dv_ * dt + 0.5f * dR_.matrix() * acc * dt * dt; // 累加更新为j时刻的位移观测dp，公式（4.16）
    dv_ = dv_ + dR_ * acc * dt;                                 // 累加更新为j时刻的速度观测dv，公式（4.13）

    // 【dR先不累乘更新为j时刻的dR，因为A, B系数阵还需要j-1时刻的dR】

    // 运动方程雅可比矩阵系数，A,B阵，见(4.29)
    // A,B阵的左上角的两项在计算得到j-1到j时刻的旋转矩阵后再填入
    Eigen::Matrix<double, 9, 9> A;
    A.setIdentity();
    Eigen::Matrix<double, 9, 6> B;
    B.setZero();

    Mat3d acc_hat = SO3::hat(acc);  // 反对称矩阵
    double dt2 = dt * dt;

    // NOTE A, B左上角块与公式稍有不同                            // 公式（4.30）
    // 系数A与B中的 dR_ 表示的是 delta R_{i,j-1}，代表的是j-1时刻的累乘量，因此在这一步之前不可以更新dR_！！！
    A.block<3, 3>(3, 0) = -dR_.matrix() * dt * acc_hat;         // A矩阵的第3行，第1列的块矩阵
    A.block<3, 3>(6, 0) = -0.5f * dR_.matrix() * acc_hat * dt2; // A矩阵的第3行，第1列的块矩阵
    A.block<3, 3>(6, 3) = dt * Mat3d::Identity();               // A矩阵的第3行，第2列的块矩阵

    B.block<3, 3>(3, 3) = dR_.matrix() * dt;                    // B矩阵的第2行，第2列的块矩阵
    B.block<3, 3>(6, 3) = 0.5f * dR_.matrix() * dt2;            // B矩阵的第3行，第2列的块矩阵

    // 更新各雅可比，见式(4.39)
    // 因为j时刻平移观测相对于零偏的雅各比递推公式中，用到了j-1时刻速度观测与旋转观测相对于零偏雅各比，
    // 且j时刻速度相对陀螺零偏的雅各比中，也用到了j-1时刻旋转相对陀螺零偏的雅各比
    // 在实现中需要：
    // (1) 先利用【j-1时刻的旋转和速度雅各比】更新【j时刻位移雅各比】
    // (2) 再利用【j-1时刻的旋转观测雅各比】更新【j时刻速度观测雅各比】
    // (3) 最后再去更新j时刻旋转相对于陀螺零偏的雅各比
    dP_dba_ = dP_dba_ + dV_dba_ * dt - 0.5f * dR_.matrix() * dt2;                      // (4.39d)
    dP_dbg_ = dP_dbg_ + dV_dbg_ * dt - 0.5f * dR_.matrix() * dt2 * acc_hat * dR_dbg_;  // (4.39e)
    dV_dba_ = dV_dba_ - dR_.matrix() * dt;                                             // (4.39b)
    dV_dbg_ = dV_dbg_ - dR_.matrix() * dt * acc_hat * dR_dbg_;                         // (4.39c)

    // 旋转部分
    Vec3d omega = gyr * dt;         // 转动量    j-1到j时刻的旋转量
    Mat3d rightJ = SO3::jr(omega);  // 右雅可比   对应的右雅可比
    SO3 deltaR = SO3::exp(omega);   // exp后     取指数映射，得到j-1到j时刻的旋转矩阵 delta R_{j-1,j}
    dR_ = dR_ * deltaR;             // (4.7a)    旋转矩阵累乘

    // 使用【j-1到j时刻】的旋转矩阵delta R_{j-1,j}更新系数矩阵A的左上角块矩阵
    A.block<3, 3>(0, 0) = deltaR.matrix().transpose();      
    // 使用【j-1到j时刻】的右雅可比更新系数矩阵的左上角块矩阵
    B.block<3, 3>(0, 0) = rightJ * dt;  

    // 更新噪声项，见式(4.31)
    cov_ = A * cov_ * A.transpose() + B * noise_gyro_acce_ * B.transpose();

    // 最后再来更新旋转相对陀螺零偏的雅各比dR_dbg
    dR_dbg_ = deltaR.matrix().transpose() * dR_dbg_ - rightJ * dt;  // (4.39a)

    // 增量积分时间
    dt_ += dt;
}

/**
 * @description: 旋转观测只和陀螺零偏bg相关
 *               当陀螺零偏发生变化时，旋转观测需要进行修正，对应公式（4.32a）
 * @param {Vec3d} &bg 当前时刻估计的陀螺零偏，bg_ 为初始零偏
 * @return {*} 修正后的旋转观测
 */
SO3 IMUPreintegration::GetDeltaRotation(const Vec3d &bg) { 
    return dR_ * SO3::exp(dR_dbg_ * (bg - bg_)); 
}

/**
 * @description: 速度观测与陀螺零偏bg和加计零偏ba都相关
 *               当陀螺和加计零偏发生变化时，速度观测需要进行修正，对应公式（4.32b）
 * @param {Vec3d} &bg 当前时刻估计的陀螺零偏，bg_ 为初始零偏
 * @param {Vec3d} &ba 当前时刻估计的加计零偏，ba_ 为初始零偏
 * @return {*} 修正后的速度观测
 */
Vec3d IMUPreintegration::GetDeltaVelocity(const Vec3d &bg, const Vec3d &ba) {
    return dv_ + dV_dbg_ * (bg - bg_) + dV_dba_ * (ba - ba_);
}

/**
 * @description: 位移观测与陀螺零偏bg和加计零偏ba都相关
 *               当陀螺和加计零偏发生变化时，平移观测需要进行修正，对应公式（4.32c）
 * @param {Vec3d} &bg 当前时刻估计的陀螺零偏，bg_ 为初始零偏
 * @param {Vec3d} &ba 当前时刻估计的加计零偏，ba_ 为初始零偏
 * @return {*} 修正后的位移观测
 */
Vec3d IMUPreintegration::GetDeltaPosition(const Vec3d &bg, const Vec3d &ba) {
    return dp_ + dP_dbg_ * (bg - bg_) + dP_dba_ * (ba - ba_);
}

/**
 * @description:  类似于ESKF中的名义状态运动学方程，用于预测下一时刻的状态
 * @param {NavStated} &start 前一时刻的状态
 * @param {Vec3d} &grav 重力矢量
 * @return {*} dt_时间后预测得到的状态
 */
NavStated IMUPreintegration::Predict(const sad::NavStated &start, const Vec3d &grav) const {
    SO3 Rj = start.R_ * dR_;
    Vec3d vj = start.R_ * dv_ + start.v_ + grav * dt_;
    Vec3d pj = start.R_ * dp_ + start.p_ + start.v_ * dt_ + 0.5f * grav * dt_ * dt_;

    auto state = NavStated(start.timestamp_ + dt_, Rj, pj, vj);
    state.bg_ = bg_;
    state.ba_ = ba_;
    return state;
}

}  // namespace sad
