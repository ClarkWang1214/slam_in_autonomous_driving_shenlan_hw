//
// Created by xiang on 22-12-29.
//

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "common/eigen_types.h"
#include "common/math_utils.h"
#include "tools/ui/pangolin_window.h"

/// 本节程序演示一个正在作圆周运动的车辆
/// 车辆的角速度与线速度可以在flags中设置

DEFINE_double(angular_velocity, 10.0, "角速度（角度）制");
DEFINE_double(linear_velocity, 5.0, "车辆前进线速度 m/s");
DEFINE_bool(use_quaternion, false, "是否使用四元数计算");

// 【新增】重力加速度大小
DEFINE_double(gravity, 9.81, "重力加速度 m/s^2");

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);

    /// 可视化
    sad::ui::PangolinWindow ui;
    if (ui.Init() == false) {
        return -1;
    }

    double angular_velocity_rad = FLAGS_angular_velocity * sad::math::kDEG2RAD;  // 弧度制角速度
    SE3 pose;                                                                    // TWB表示的位姿
    Vec3d omega(0, 0, angular_velocity_rad);                                     // 角速度矢量

    // v_body本体系速度，是世界系速度矢量转换到车辆坐标系下的结果，可以被各种传感器测量到
    Vec3d v_body(FLAGS_linear_velocity, 0, 0);                                   

    // 【新增】世界系下的加速度
    Vec3d a_world(0, 0, -FLAGS_gravity);                                         
    const double dt = 0.05;                                                      // 每次更新的时间

    Vec3d v_world, a_body;

    while (ui.ShouldQuit() == false) {
        // 本体系下的速度转换到世界系下 v_w = R_wb * v_b
        // Vec3d v_world = pose.so3() * v_body;
        // pose.translation() += v_world * dt;

        // 世界系加速度对平移的影响
        // p_w = p_w + v_w * dt + 0.5 * a_w * dt * dt
        v_world = pose.so3() * v_body;                  // 将本体系速度转换到世界系下
        pose.translation() += v_world * dt              // 速度乘以时间，得到位移
                           + 0.5 * a_world * dt * dt;   // 世界系加速度对平移的影响

        // 世界系下的加速度转换到本体系 
        // a_b = R_wb^T * a_w = R_bw * a_w
        Vec3d a_body = pose.so3().inverse() * a_world;

        //更新本体系下的速度
        // v_b = v_b + a_b * dt
        v_body += a_body * dt; 

        // 更新自身旋转
        if (FLAGS_use_quaternion) {
            // 四元数指数映射 q = exp(1/2*[0, \vec{omega}]^T) = Exp(omega)
            // 
            //      q_omega    =    [cos(0.5 * omega * dt), sin(0.5 * \vec{omega} * dt)]^T
            // omega比较小时  约等于 [1, 0.5 * \vec{omega} * dt]^T
            //
            // 对应于公式（3.15a）
            Quatd q = pose.unit_quaternion() * Quatd(1, 0.5 * omega[0] * dt, 0.5 * omega[1] * dt, 0.5 * omega[2] * dt);
            q.normalize();
            pose.so3() = SO3(q); // 四元数转换为旋转矩阵
        } else {
            // pose.so3为Rwb，表示本体系到世界系的旋转
            pose.so3() = pose.so3() * SO3::exp(omega * dt);
        }

        LOG(INFO) << "a_body: " << a_body.transpose();
        LOG(INFO) << "v_world: " << v_world.transpose();
        LOG(INFO) << "v_body: " << v_body.transpose();
        LOG(INFO) << "pose: " << pose.translation().transpose();
        ui.UpdateNavState(sad::NavStated(0, pose, v_world));

        usleep(dt * 1e6);
    }

    ui.Quit();
    return 0;
}