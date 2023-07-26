//
// Created by xiang on 2021/7/19.
//

#include "ch4/gins_pre_integ.h"
#include "ch4/g2o_types.h"
#include "common/g2o_types.h"

#include <glog/logging.h>

#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/robust_kernel.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>

namespace sad {

/**
 * @description: 添加IMU数据
 * @param {IMU&} imu
 * @return {*}
 */
void GinsPreInteg::AddImu(const IMU& imu) {
    // 判断是否接收到第一个GNSS数据以及第一个IMU数据
    // 预积分需要从第一个接收到GNSS数据之后的IMU数据开始累计
    if (first_gnss_received_ && first_imu_received_) 
        pre_integ_->Integrate(imu, imu.timestamp_ - last_imu_.timestamp_);

    first_imu_received_ = true;     // 标记已收到第一个IMU数据
    last_imu_ = imu;                // 更新上一帧IMU数据
    current_time_ = imu.timestamp_; // 更新当前时间为IMU时间戳
    // LOG(INFO) << std::fixed << std::setprecision(9)<< "imu.timestamp_: " << imu.timestamp_;
}

/**
 * @description: 配置项
 * @param {Options} options
 * @return {*}
 */
void GinsPreInteg::SetOptions(sad::GinsPreInteg::Options options) {

    // 【新增】设置配置项
    options_ = options; 

    double bg_rw2 = 1.0 / (options_.bias_gyro_var_ * options_.bias_gyro_var_);
    options_.bg_rw_info_.diagonal() << bg_rw2, bg_rw2, bg_rw2; // 陀螺零偏随机游走信息矩阵，用于g2o优化
    double ba_rw2 = 1.0 / (options_.bias_acce_var_ * options_.bias_acce_var_);
    options_.ba_rw_info_.diagonal() << ba_rw2, ba_rw2, ba_rw2; // 加计零偏随机游走信息矩阵，用于g2o优化

    // GNSS位移、高度、角度方差
    double gp2 = options_.gnss_pos_noise_ * options_.gnss_pos_noise_; 
    double gh2 = options_.gnss_height_noise_ * options_.gnss_height_noise_;
    double ga2 = options_.gnss_ang_noise_ * options_.gnss_ang_noise_;

    // GNSS信息矩阵，方差倒数
    options_.gnss_info_.diagonal() << 1.0 / ga2, 1.0 / ga2, 1.0 / ga2, 1.0 / gp2, 1.0 / gp2, 1.0 / gh2; 
    
    // 预积分类初始化
    pre_integ_ = std::make_shared<IMUPreintegration>(options_.preinteg_options_);

    double o2 = 1.0 / (options_.odom_var_ * options_.odom_var_); // 轮速计方差的倒数
    options_.odom_info_.diagonal() << o2, o2, o2; // 轮速计信息矩阵

    // 用于g2o优化中先验边的信息矩阵设置
    prior_info_.block<6, 6>(9, 9) = Mat6d ::Identity() * 1e6;

    if (this_frame_) {
        this_frame_->bg_ = options_.preinteg_options_.init_bg_;
        this_frame_->ba_ = options_.preinteg_options_.init_ba_;
    }
}

/**
 * @description: 添加GNSS数据
 * @param {GNSS&} gnss
 * @return {*}
 */
void GinsPreInteg::AddGnss(const GNSS& gnss) {
    this_frame_ = std::make_shared<NavStated>(current_time_);   // 接收到GNSS数据
    this_gnss_ = gnss;

    last_gnss_set_ = true;

    // 判断接收到的是否为首个GNSS数据
    if (!first_gnss_received_) {
        // 判断首个GNSS数据航向是否有效
        if (!this_gnss_.heading_valid_)  // 默认为false
            return;  // 要求【首个GNSS必须有航向】，若没有航向，直接返回

        // 此时是转换成功且有航向数据的有效GNSS数据

        // 首个gnss信号，将初始pose设置为该gnss信号
        // 将第一个GNSS数据的位姿作为当前帧状态的初始位姿
        this_frame_->timestamp_ = gnss.unix_time_;                  // 取GNSS时间戳作为当前帧时间戳
        this_frame_->p_ = gnss.utm_pose_.translation();             // 来自首个带姿态的GNSS信号的平移
        this_frame_->R_ = gnss.utm_pose_.so3();                     // 来自首个带姿态的GNSS信号的旋转
        this_frame_->v_.setZero();                                  // 速度设为0
        this_frame_->bg_ = options_.preinteg_options_.init_bg_;     // 静止初始化得到的陀螺零偏
        this_frame_->ba_ = options_.preinteg_options_.init_ba_;     // 静止初始化得到的加计零偏

        // 初始化预积分器实例对象
        pre_integ_ = std::make_shared<IMUPreintegration>(options_.preinteg_options_);

        last_frame_ = this_frame_;
        last_gnss_ = this_gnss_;
        first_gnss_received_ = true;
        current_time_ = gnss.unix_time_;
        return; // 第一次接收到GNSS数据后，不做预积分，直接返回，去完成IMU数据的累计，等待下一次GNSS数据到来
    }

    // 预积分前一个GNSS到当前GNSS时刻之间的IMU数据
    // 第一次的时候，last_imu_中的加速度和角速度都是0，时间戳也是0.0，此时还没有IMU数据，预积分函数中的量都是0
    pre_integ_->Integrate(last_imu_, gnss.unix_time_ - current_time_); // 当前GNSS数据时间戳减去上一帧GNSS数据时间戳

    current_time_ = gnss.unix_time_; // 更新当前时间为接收到GNSS数据的时间戳

    // 根据上一帧的GNSS状态数据预测当前帧GNSS的状态
    // 接收到首个GNSS的时候，this_frame就是last_frame，因为预积分量都是0
    *this_frame_ = pre_integ_->Predict(*last_frame_, options_.gravity_);
    
    Optimize(); // 触发一次g2o优化
    // LOG(INFO) << std::fixed << std::setprecision(9)<< "gnss.unix_time_: " << gnss.unix_time_;

    last_frame_ = this_frame_;  // 更新上一帧GNSS数据为当前帧GNSS数据
    last_gnss_ = this_gnss_;    // 更新上一帧GNSS数据为当前帧GNSS数据
}

/**
 * @description: Odom轮速数据到来时也触发一次优化，没有j时刻的GNSS约束，只有j时刻的Odom的速度约束，看优化函数应该这么改
 *               如果Odom不触发，在RTK中GNSS数据来的时候触发一次优化，主要是R,t,v部分的更新
 * @param {Odom&} odom
 * @return {*}
 */
void GinsPreInteg::AddOdom(const sad::Odom& odom) {
    last_odom_ = odom;
    last_odom_set_ = true;

    // 累计IMU到当前Odom时刻的预积分量
    pre_integ_->Integrate(last_imu_, odom.timestamp_ - current_time_);

    // 更新当前时间为接收到的轮速计的时间戳
    current_time_ = odom.timestamp_; 

    *this_frame_ = pre_integ_->Predict(*last_frame_, options_.gravity_);
    // LOG(INFO) << "odom.timestamp_: " << odom.timestamp_;
    // LOG(INFO) << std::fixed << std::setprecision(9)<< "odom.timestamp_: " << odom.timestamp_;

    // // g2o优化
    Optimize();
    
    last_frame_ = this_frame_;  // 更新状态
}

/**
 * @description: g2o优化
 * @return {*}
 */
void GinsPreInteg::Optimize() {
    if (pre_integ_->dt_ < 1e-2)  // 1624426913.27996874 - 1624426913.27186489 = 0.00810385 比1e-3 = 0.01大
        // 未得到积分             // 1624426913.9799695 - 1624426913.9718504 = 0.0081191 比1e-3 = 0.01大
        return;                  //
    

    //////////////////////////////// 构建图优化问题 ////////////////////////////////

    using BlockSolverType = g2o::BlockSolverX; 
    using LinearSolverType = g2o::LinearSolverEigen<BlockSolverType::PoseMatrixType>; // 线性求解器

    // 梯度下降算法：列文伯格马夸尔特LM算法（可选：GN、LM、DogLeg）
    auto* solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer; // 稀疏优化器
    optimizer.setAlgorithm(solver); // 设置求解器

    //////////////////////////////// 添加顶点 ////////////////////////////////

    // 上一时刻的顶点， pose, v, bg, ba
    auto v0_pose = new VertexPose();                // 前一时刻的位姿顶点
    v0_pose->setId(0);                              // 设置id           
    v0_pose->setEstimate(last_frame_->GetSE3());    // 设置估计值，上一帧的位姿（R, p）
    optimizer.addVertex(v0_pose);                   // 添加位姿顶点，注意添加顺序

    auto v0_vel = new VertexVelocity();             // 前一时刻的速度顶点
    v0_vel->setId(1);                               // 设置id
    v0_vel->setEstimate(last_frame_->v_);           // 设置估计值，上一帧的速度
    optimizer.addVertex(v0_vel);                    // 添加速度顶点，注意添加顺序

    auto v0_bg = new VertexGyroBias();              // 前一时刻的陀螺零偏顶点
    v0_bg->setId(2);                                // 设置id
    v0_bg->setEstimate(last_frame_->bg_);           // 设置估计值，上一帧的陀螺零偏
    optimizer.addVertex(v0_bg);                     // 添加陀螺零偏顶点，注意添加顺序

    auto v0_ba = new VertexAccBias();               // 前一时刻的加计零偏顶点
    v0_ba->setId(3);                                // 设置id
    v0_ba->setEstimate(last_frame_->ba_);           // 设置估计值，上一帧的加计零偏
    optimizer.addVertex(v0_ba);                     // 添加加计零偏顶点，注意添加顺序

    // 本时刻顶点，pose, v, bg, ba
    auto v1_pose = new VertexPose();                // 当前时刻的位姿顶点
    v1_pose->setId(4);                              // 设置id
    v1_pose->setEstimate(this_frame_->GetSE3());    // 设置估计值，当前帧的位姿（R, p）
    optimizer.addVertex(v1_pose);                   // 添加位姿顶点，注意添加顺序

    auto v1_vel = new VertexVelocity();             // 当前时刻的速度顶点
    v1_vel->setId(5);                               // 设置id
    v1_vel->setEstimate(this_frame_->v_);           // 设置估计值，当前帧的速度
    optimizer.addVertex(v1_vel);                    // 添加速度顶点，注意添加顺序

    auto v1_bg = new VertexGyroBias();              // 当前时刻的陀螺零偏顶点
    v1_bg->setId(6);                                // 设置id
    v1_bg->setEstimate(this_frame_->bg_);           // 设置估计值，当前帧的陀螺零偏
    optimizer.addVertex(v1_bg);                     // 添加陀螺零偏顶点，注意添加顺序

    auto v1_ba = new VertexAccBias();               // 当前时刻的加计零偏顶点
    v1_ba->setId(7);                                // 设置id
    v1_ba->setEstimate(this_frame_->ba_);           // 设置估计值，当前帧的加计零偏
    optimizer.addVertex(v1_ba);                     // 添加加计零偏顶点，注意添加顺序

    //////////////////////////////// 添加边 ////////////////////////////////

    // i和j之间的预积分边类对象
    auto edge_inertial = new EdgeInertial(pre_integ_, options_.gravity_); // 设置预积分噪声协方差矩阵的逆作为信息矩阵
    edge_inertial->setVertex(0, v0_pose);           // 前一时刻的位姿（R0, p0）
    edge_inertial->setVertex(1, v0_vel);            // 前一时刻的速度 v0
    edge_inertial->setVertex(2, v0_bg);             // 前一时刻的陀螺零偏 bg0
    edge_inertial->setVertex(3, v0_ba);             // 前一时刻的加计零偏 ba0
    edge_inertial->setVertex(4, v1_pose);           // 当前时刻的位姿（R1, p1）
    edge_inertial->setVertex(5, v1_vel);            // 当前时刻的速度 v1
    auto* rk = new g2o::RobustKernelHuber();    
    rk->setDelta(200.0);
    edge_inertial->setRobustKernel(rk);             // 设置鲁棒核函数
    optimizer.addEdge(edge_inertial);               // 添加预积分边到图优化器中

    // i和j之间的陀螺零偏随机游走边
    auto* edge_gyro_rw = new EdgeGyroRW();
    edge_gyro_rw->setVertex(0, v0_bg);                  // 前一时刻的陀螺零偏 bg0
    edge_gyro_rw->setVertex(1, v1_bg);                  // 当前时刻的陀螺零偏 bg1
    edge_gyro_rw->setInformation(options_.bg_rw_info_); // 设置信息矩阵
    optimizer.addEdge(edge_gyro_rw);                    // 往图优化器中添加边

    // i和j之间的加计零偏随机游走边
    auto* edge_acc_rw = new EdgeAccRW();
    edge_acc_rw->setVertex(0, v0_ba);                   // 前一时刻的加计零偏 ba0
    edge_acc_rw->setVertex(1, v1_ba);                   // 当前时刻的加计零偏 ba1
    edge_acc_rw->setInformation(options_.ba_rw_info_);  // 设置信息矩阵
    optimizer.addEdge(edge_acc_rw);                     // 往图优化器中添加边

    // i时刻的先验边
    auto* edge_prior = new EdgePriorPoseNavState(*last_frame_, prior_info_);
    edge_prior->setVertex(0, v0_pose);                  // 前一时刻的位姿（R0, p0）
    edge_prior->setVertex(1, v0_vel);                   // 前一时刻的速度 v0
    edge_prior->setVertex(2, v0_bg);                    // 前一时刻的陀螺零偏 bg0
    edge_prior->setVertex(3, v0_ba);                    // 前一时刻的加计零偏 ba0
    optimizer.addEdge(edge_prior);                      // 往图优化器中添加边

    // i时刻的GNSS边
    auto edge_gnss0 = new EdgeGNSS(v0_pose, last_gnss_.utm_pose_); // 添加前一时刻的位姿顶点到GNSS边中
    edge_gnss0->setInformation(options_.gnss_info_); // 设置信息矩阵：协方差矩阵的逆
    optimizer.addEdge(edge_gnss0);   // 往图优化器中添加边

    // Odom边  当前里程计边
    EdgeEncoder3D* edge_odom = nullptr;
    Vec3d vel_world = Vec3d::Zero();
    Vec3d vel_odom = Vec3d::Zero();
    EdgeGNSS* edge_gnss1 = nullptr;

    if(last_gnss_set_){
        // j时刻的GNSS边
        edge_gnss1 = new EdgeGNSS(v1_pose, this_gnss_.utm_pose_);
        edge_gnss1->setInformation(options_.gnss_info_); // 设置信息矩阵：协方差矩阵的逆
        optimizer.addEdge(edge_gnss1);  // 往图优化器中添加边
        last_gnss_set_ = false;
    }
    
    // // 判断是否有最新的Odom里程计数据到达
    // if (last_odom_set_) {
        // 若已到达，记录此时速度观测velocity obs，并保留读数
        // v_wheel = 2 * pi * r * pulse / circle_pulse / odom_span
        // 其中，r为轮子半径，pulse为编码器读数，circle_pulse为编码器每圈脉冲数，odom_span为里程计测量间隔
        double velo_l = options_.wheel_radius_ * last_odom_.left_pulse_ / options_.circle_pulse_ * 2 * M_PI / options_.odom_span_; // 公式（3.76）
        double velo_r = options_.wheel_radius_ * last_odom_.right_pulse_ / options_.circle_pulse_ * 2 * M_PI / options_.odom_span_;
        double average_vel = 0.5 * (velo_l + velo_r); // 计算平均速度
        vel_odom = Vec3d(average_vel, 0.0, 0.0);    // 里程计测量的速度，只有x方向，假设y、z两个轴上没有速度测量
        vel_world = this_frame_->R_ * vel_odom;     // 本体系转换到世界系

        // 添加j时刻的速度顶点到Odom边中，权重为轮速计观测到的世界系下左右轮的平均速度
        edge_odom = new EdgeEncoder3D(v1_vel, vel_world);   
        edge_odom->setInformation(options_.odom_info_);     // 设置里Odom的信息矩阵
        optimizer.addEdge(edge_odom);                       // 往图优化器中添加Odom观测边

        // 重置odom数据到达标志位，等待最新的odom轮速计数据
        last_odom_set_ = false;
    // }

    //////////////////////////////// 优化 ////////////////////////////////

    // optimizer.setVerbose(options_.verbose_);    // 设置是否输出调试信息
    optimizer.setVerbose(false);    // 设置是否输出调试信息
    optimizer.initializeOptimization();         // 初始化优化器
    optimizer.optimize(20);                     // 进行20次优化

    // 判断是否输出调试信息
    // if (options_.verbose_) {
    if (false) {
        // 获取结果，统计各类误差
        LOG(INFO) << "chi2/error: ";
        LOG(INFO) << "preintegration: " << edge_inertial->chi2() << "/" << edge_inertial->error().transpose();
        LOG(INFO) << "gnss0: " << edge_gnss0->chi2() << ", " << edge_gnss0->error().transpose();
        // LOG(INFO) << "gnss1: " << edge_gnss1->chi2() << ", " << edge_gnss1->error().transpose();
        LOG(INFO) << "bias: " << edge_gyro_rw->chi2() << "/" << edge_acc_rw->error().transpose();
        
        // LOG(INFO) << "p_1 estimated: " << v1_pose->estimate().translation().transpose();
        // LOG(INFO) << "v_1 estimated: " << v1_vel->estimate().transpose();
        // LOG(INFO) << "bg_1 estimated: " << v1_bg->estimate().transpose();
        // LOG(INFO) << "ba_1 estimated: " << v1_ba->estimate().transpose();

        LOG(INFO) << "prior: " << edge_prior->chi2() << "/" << edge_prior->error().transpose();
        if (edge_odom) {
            LOG(INFO) << "body vel: " << (v1_pose->estimate().so3().inverse() * v1_vel->estimate()).transpose();
            LOG(INFO) << "meas: " << vel_odom.transpose();
            LOG(INFO) << "odom: " << edge_odom->chi2() << "/" << edge_odom->error().transpose();
        }
    }

    // 保存g2o图优化估计的前一帧状态量
    last_frame_->R_ = v0_pose->estimate().so3();            // 旋转 
    last_frame_->p_ = v0_pose->estimate().translation();    // 平移
    last_frame_->v_ = v0_vel->estimate();                   // 速度
    last_frame_->bg_ = v0_bg->estimate();                   // 陀螺零偏
    last_frame_->ba_ = v0_ba->estimate();                   // 加计零偏

    // 保存g2o图优化估计的当前帧状态量
    this_frame_->R_ = v1_pose->estimate().so3();            // 旋转
    this_frame_->p_ = v1_pose->estimate().translation();    // 平移
    this_frame_->v_ = v1_vel->estimate();                   // 速度
    this_frame_->bg_ = v1_bg->estimate();                   // 陀螺零偏
    this_frame_->ba_ = v1_ba->estimate();                   // 加计零偏

    // 使用当前帧估计的陀螺零偏、加速度计零偏作为下一次预积分的初始零偏
    options_.preinteg_options_.init_bg_ = this_frame_->bg_; 
    options_.preinteg_options_.init_ba_ = this_frame_->ba_; 

    // 每次只对i到j时刻之间的IMU数据进行预积分，g2o图优化结束后，重置预积分器的状态，方便下次触发优化时使用
    pre_integ_ = std::make_shared<IMUPreintegration>(options_.preinteg_options_);
}

/**
 * @description: 获取预计分的预测状态
 * @return {*}
 */
NavStated GinsPreInteg::GetState() const {
    // 判断当前帧是否为空
    if (this_frame_ == nullptr) return {};
    
    // 判断预积分器对象是否为空，若是，直接返回当前帧状态
    if (pre_integ_ == nullptr)  return *this_frame_;

    // 根据当前帧的状态，预测dt_后的状态
    return pre_integ_->Predict(*this_frame_, options_.gravity_);
}
}  // namespace sad