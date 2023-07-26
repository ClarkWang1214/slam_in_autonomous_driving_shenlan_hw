//
// Created by xiang on 23-1-19.
//

#include "ch4/g2o_types.h"
#include "common/g2o_types.h"

namespace sad {

/**
 * @description: 构造函数中需要指定预积分类对象
 * @param {shared_ptr<IMUPreintegration>} preinteg 预积分对象指针
 * @param {Vec3d&} gravity 重力矢量
 * @param {double} weight 权重
 * @return {*}
 */
EdgeInertial::EdgeInertial(std::shared_ptr<IMUPreintegration> preinteg, const Vec3d& gravity, double weight) : preint_(preinteg), dt_(preinteg->dt_) {
    resize(6);  // 6个关联顶点（位姿顶点i、速度顶点i、陀螺零偏顶点i、加计零偏顶点i、位姿顶点j、速度顶点j）
    grav_ = gravity; // 重力矢量

    // 使用预积分中累加的噪声协方差矩阵的逆作为预积分边的信息矩阵
    setInformation(preinteg->cov_.inverse() * weight); // 默认权重weight为1.0
}

/**
 * @description: 取出边所连接的顶点的当前估计值，与零偏更新后的预积分观测值进行比较，计算三个残差
 * @return {*}
 */
void EdgeInertial::computeError() {
    auto* p1 = dynamic_cast<const VertexPose*>(_vertices[0]);       // i时刻的位姿 Ri, pi
    auto* v1 = dynamic_cast<const VertexVelocity*>(_vertices[1]);   // i时刻的速度 vi
    auto* bg1 = dynamic_cast<const VertexGyroBias*>(_vertices[2]);  // i时刻的陀螺零偏 bgi
    auto* ba1 = dynamic_cast<const VertexAccBias*>(_vertices[3]);   // i时刻的加计零偏 bai
    auto* p2 = dynamic_cast<const VertexPose*>(_vertices[4]);       // j时刻的位姿 Rj, pj
    auto* v2 = dynamic_cast<const VertexVelocity*>(_vertices[5]);   // j时刻的速度 vj

    const SO3 R1 = p1->estimate().so3();                // i时刻估计的旋转
    const SO3 R2 = p2->estimate().so3();                // j时刻估计的旋转
    const Vec3d vi = v1->estimate();                    // i时刻估计的速度
    const Vec3d vj = v2->estimate();                    // j时刻估计的速度
    const Vec3d pi = p1->estimate().translation();      // i时刻估计的位置
    const Vec3d pj = p2->estimate().translation();      // j时刻估计的位置

    Vec3d bg = bg1->estimate(); // 估计的陀螺零偏
    Vec3d ba = ba1->estimate(); // 估计的加计零偏

    const SO3 dR = preint_->GetDeltaRotation(bg);
    const Vec3d dv = preint_->GetDeltaVelocity(bg, ba);
    const Vec3d dp = preint_->GetDeltaPosition(bg, ba);

    /// 预积分误差项
    // 旋转残差，见式(4.41a)，Log(Delta_R_ij^T * R_i^T * R_j)
    const Vec3d er = (dR.inverse() * R1.inverse() * R2).log();
    // 计算i时刻旋转矩阵的转置，保存为中间变量RiT
    const Mat3d RiT = R1.inverse().matrix();
    // 速度残差，见式(4.41b)，R_i^T * (v_j - v_i - g * dt_ij) - Delta_v_ij
    const Vec3d ev = RiT * (vj - vi - grav_ * dt_) - dv;
    // 位置残差，见式(4.41c)，R_i^T * (p_j - p_i - v_i * dt_ij - 1/2 * g * dt_ij^2) - Delta_p_ij
    const Vec3d ep = RiT * (pj - pi - vi * dt_ - grav_ * dt_ * dt_ / 2) - dp;

    _error << er, ev, ep;
}

/**
 * @description: 计算三条边（旋转残差、速度残差、平移残差）相对于六个顶点（SE3位姿、速度、陀螺零偏、加计零偏）的雅可比。
 * @return {*}
 */
void EdgeInertial::linearizeOplus() {
    auto* p1 = dynamic_cast<const VertexPose*>(_vertices[0]);       // i时刻的位姿 Ri, pi
    auto* v1 = dynamic_cast<const VertexVelocity*>(_vertices[1]);   // i时刻的速度 vi
    auto* bg1 = dynamic_cast<const VertexGyroBias*>(_vertices[2]);  // i时刻的陀螺零偏 bgi
    auto* ba1 = dynamic_cast<const VertexAccBias*>(_vertices[3]);   // i时刻的加计零偏 bai
    auto* p2 = dynamic_cast<const VertexPose*>(_vertices[4]);       // j时刻的位姿 Rj, pj
    auto* v2 = dynamic_cast<const VertexVelocity*>(_vertices[5]);   // j时刻的速度 vj

    // double my_update_[] = {1e-9,1e-9,1e-9};
    // // p1->oplus(my_update_);
    // g2o::dynamic_aligned_buffer<number_t> buffer{ 12 };
    const number_t delta = g2o::cst(1e-9);
    const number_t scalar = 1 / (2*delta);

    const int p1_dim = p1->dimension();     // 6
    const int v1_dim = v1->dimension();     // 3
    const int bg1_dim = bg1->dimension();   // 3
    const int ba1_dim = ba1->dimension();   // 3
    const int p2_dim = p2->dimension();     // 6
    const int v2_dim = v2->dimension();     // 3

    // 一些中间符号
    const SO3 R1 = p1->estimate().so3();                // i时刻估计的旋转
    const SO3 R2 = p2->estimate().so3();                // j时刻估计的旋转
    const Vec3d vi = v1->estimate();                    // i时刻估计的速度
    const Vec3d vj = v2->estimate();                    // j时刻估计的速度
    const Vec3d pi = p1->estimate().translation();      // i时刻估计的位置
    const Vec3d pj = p2->estimate().translation();      // j时刻估计的位置
    
    const SO3 R1T = R1.inverse();                       // i时刻估计的旋转的转置

    const Vec3d bg = bg1->estimate();                   // 估计的陀螺零偏
    const Vec3d ba = ba1->estimate();                   // 估计的加计零偏
    const Vec3d dbg = bg - preint_->bg_;                // 当前估计的陀螺零偏减去初始零偏，得到陀螺零偏的变化量

    // 预积分观测量相对于IMU零偏的雅可比
    auto dR_dbg = preint_->dR_dbg_;                     // 预积分【旋转】观测量相对于【陀螺】零偏的雅可比
    auto dv_dbg = preint_->dV_dbg_;                     // 预积分【速度】观测量相对于【陀螺】零偏的雅可比
    auto dp_dbg = preint_->dP_dbg_;                     // 预积分【平移】观测量相对于【陀螺】零偏的雅可比
    auto dv_dba = preint_->dV_dba_;                     // 预积分【速度】观测量相对于【加计】零偏的雅可比
    auto dp_dba = preint_->dP_dba_;                     // 预积分【平移】观测量相对于【加计】零偏的雅可比

    const SO3 dR = preint_->GetDeltaRotation(bg);       // 根据估计出的陀螺零偏修正预积分的旋转观测量
    const SO3 eR = SO3(dR).inverse() * R1T * R2;        // 公式(4.41a)中的旋转残差 \Delta R_ij^T * R_i^T * R_j
    const Vec3d er = eR.log();                          // 转为旋转向量 
    const Mat3d invJr = SO3::jr_inv(eR);                // 旋转矩阵的右雅可比矩阵的逆矩阵

    /// 雅可比矩阵
    /// 注意有3个index, 顶点的，自己误差的，顶点内部变量的
    /// 变量顺序：pose1(R1,p1), v1, bg1, ba1, pose2(R2,p2), v2
    /// 残差顺序：eR, ev, ep，残差顺序为行，变量顺序为列

    //       | R1 | p1 | v1 | bg1 | ba1 | R2 | p2 | v2 |
    //  vert | 0       | 1  | 2   | 3   | 4       | 5  |
    //  col  | 0    3  | 0  | 0   | 0   | 0    3  | 0  |
    //    row
    //  eR 0 |
    //  ev 3 |
    //  ep 6 |

    /// 残差对R1, 9x3
    _jacobianOplus[0].setZero();
    // 旋转残差相对于R1的雅可比dR/dR1，见公式（4.42）
    _jacobianOplus[0].block<3, 3>(0, 0) = -invJr * (R2.inverse() * R1).matrix();
    // 速度残差相对于R1的雅可比dv/dR1，见公式（4.47）
    _jacobianOplus[0].block<3, 3>(3, 0) = SO3::hat(R1T * (vj - vi - grav_ * dt_));
    // 平移残差相对于R1的雅可比dp/dR1，见公式（4.48）
    _jacobianOplus[0].block<3, 3>(6, 0) = SO3::hat(R1T * (pj - pi - vi * dt_ - 0.5 * grav_ * dt_ * dt_));

    /// 残差对p1, 9x3
    // dp/dp1, 4.48a
    // 平移残差相对于p1的雅可比dp/dp1
    _jacobianOplus[0].block<3, 3>(6, 3) = -R1T.matrix();

    // const Mat3d RiT_clark = -R1T.matrix();
    // LOG(INFO) << "-R1T.matrix(): " << RiT_clark;
    
    /// 残差对v1, 9x3
    _jacobianOplus[1].setZero();
    // 速度残差相对于v1的雅可比dv/dv1
    _jacobianOplus[1].block<3, 3>(3, 0) = -R1T.matrix();
    // 平移残差相对于v1的雅可比dp/dv1, 4.48c
    _jacobianOplus[1].block<3, 3>(6, 0) = -R1T.matrix() * dt_;

    /// 残差对bg1
    _jacobianOplus[2].setZero();
    // 旋转残差相对于陀螺零偏dbg1的雅可比dR/dbg1, 4.45
    _jacobianOplus[2].block<3, 3>(0, 0) = -invJr * eR.inverse().matrix() * SO3::jr((dR_dbg * dbg).eval()) * dR_dbg;
    // 速度残差相对于陀螺零偏dbg1的雅可比dv/dbg1
    _jacobianOplus[2].block<3, 3>(3, 0) = -dv_dbg;
    // 平移残差相对于陀螺零偏dbg1的雅可比dp/dbg1
    _jacobianOplus[2].block<3, 3>(6, 0) = -dp_dbg;

    /// 残差对ba1
    _jacobianOplus[3].setZero();
    // 速度残差相对于加计零偏dba1的雅可比dv/dba1
    _jacobianOplus[3].block<3, 3>(3, 0) = -dv_dba;
    // 平移残差相对于加计零偏dba1的雅可比dp/dba1
    _jacobianOplus[3].block<3, 3>(6, 0) = -dp_dba;

    /// 残差对pose2
    _jacobianOplus[4].setZero();
    // 旋转残差相对于R2的雅可比dR/dR2, 4.43
    _jacobianOplus[4].block<3, 3>(0, 0) = invJr;
    // 平移残差相对于P2的雅可比dp/dp2, 4.48b
    _jacobianOplus[4].block<3, 3>(6, 3) = R1T.matrix();

    /// 残差对v2
    _jacobianOplus[5].setZero();
    // 速度残差相对于v2的雅可比dv/dv2, 4,46b
    _jacobianOplus[5].block<3, 3>(3, 0) = R1T.matrix();  // OK

    // 【第三题数值求导】
    double my_update_[] = {0,0,1e-9};
    // SO3 R1 = R1 * SO3::exp(Eigen::Map<const Vec3d>(&my_update_[0])); // 给R1加一个微小扰动
    // SO3 R2_perturb = R2 * SO3::exp(Eigen::Map<const Vec3d>(&my_update_[0])); // 给R1加一个微小扰动
    Vec3d pi_perturb = pi + Eigen::Map<const Vec3d>(&my_update_[0]);  
    // Vec3d pj_perturb = pj + Eigen::Map<const Vec3d>(&my_update_[0]);   
    // Vec3d vi_perturb = vi + Eigen::Map<const Vec3d>(&my_update_[0]);  
    // Vec3d vj_perturb = vj + Eigen::Map<const Vec3d>(&my_update_[0]);   
    // Vec3d bg_perturb = bg + Eigen::Map<const Vec3d>(&my_update_[0]);  
    // Vec3d ba_perturb = ba + Eigen::Map<const Vec3d>(&my_update_[0]);        

    // 将当前估计的零偏，与初始零偏相减，得到零偏的增量 delta bg, delta ba
    SO3 dR1 = preint_->GetDeltaRotation(bg);// 修正后的预积分旋转观测量，对应式(4.32a)
    Vec3d dv = preint_->GetDeltaVelocity(bg, ba); // 修正后的预积分速度观测量，对应式(4.32b)
    Vec3d dp = preint_->GetDeltaPosition(bg, ba); // 修正后的预积分平移观测量，对应式(4.32c)

    /// 预积分误差项
    // 旋转残差，见式(4.41a)，Log(Delta_R_ij^T * R_i^T * R_j)
    Vec3d er_perturb = (dR1.inverse() * R1.inverse() * R2).log();
    // 计算i时刻旋转矩阵的转置，保存为中间变量RiT
    Mat3d RiT = R1.inverse().matrix();
    // 速度残差，见式(4.41b)，R_i^T * (v_j - v_i - g * dt_ij) - Delta_v_ij
    Vec3d ev_perturb = RiT * (vj - vi - grav_ * dt_) - dv;
    // 位置残差，见式(4.41c)，R_i^T * (p_j - p_i - v_i * dt_ij - 1/2 * g * dt_ij^2) - Delta_p_ij
    Vec3d ep_perturb = RiT * (pj - pi_perturb - vi * dt_ - grav_ * dt_ * dt_ / 2) - dp; 

    ErrorVector errorBak;
    errorBak << er_perturb, ev_perturb, ep_perturb;

    double my_update_2[] = {0,0,-1e-9};
    // R1_perturb = R1 * SO3::exp(Eigen::Map<const Vec3d>(&my_update_2[0])); // 给R1加一个微小扰动
    // R2_perturb = R2 * SO3::exp(Eigen::Map<const Vec3d>(&my_update_2[0])); // 给R1加一个微小扰动
    pi_perturb = pi + Eigen::Map<const Vec3d>(&my_update_2[0]);  
    // pj_perturb = pj + Eigen::Map<const Vec3d>(&my_update_2[0]);   
    // vi_perturb = vi + Eigen::Map<const Vec3d>(&my_update_2[0]);  
    // vj_perturb = vj + Eigen::Map<const Vec3d>(&my_update_2[0]);   
    // bg_perturb = bg + Eigen::Map<const Vec3d>(&my_update_2[0]);  
    // ba_perturb = ba + Eigen::Map<const Vec3d>(&my_update_2[0]);        

    // RiT_perturb = R1_perturb.inverse().matrix();

    // 将当前估计的零偏，与初始零偏相减，得到零偏的增量 delta bg, delta ba
    // dR_perturb = preint_->GetDeltaRotation(bg_perturb);// 修正后的预积分旋转观测量，对应式(4.32a)
    // dv_perturb = preint_->GetDeltaVelocity(bg_perturb, ba_perturb); // 修正后的预积分速度观测量，对应式(4.32b)
    // dp_perturb = preint_->GetDeltaPosition(bg_perturb, ba_perturb); // 修正后的预积分平移观测量，对应式(4.32c)

    /// 预积分误差项
    // 旋转残差，见式(4.41a)，Log(Delta_R_ij^T * R_i^T * R_j)
    er_perturb = (dR1.inverse() * R1.inverse() * R2).log();
    // 计算i时刻旋转矩阵的转置，保存为中间变量RiT
    RiT = R1.inverse().matrix();
    // 速度残差，见式(4.41b)，R_i^T * (v_j - v_i - g * dt_ij) - Delta_v_ij
    ev_perturb = RiT * (vj - vi - grav_ * dt_) - dv;
    // 位置残差，见式(4.41c)，R_i^T * (p_j - p_i - v_i * dt_ij - 1/2 * g * dt_ij^2) - Delta_p_ij
    ep_perturb = RiT * (pj - pi_perturb - vi * dt_ - grav_ * dt_ * dt_ / 2) - dp; 

    ErrorVector errorBak2;
    errorBak2 << er_perturb, ev_perturb, ep_perturb;

    errorBak = errorBak - errorBak2;

    errorBak = scalar * errorBak; 

    // LOG(INFO) << "errorBak: " << errorBak.transpose();
}

}  // namespace sad