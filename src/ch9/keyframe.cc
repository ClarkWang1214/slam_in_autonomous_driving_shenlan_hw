//
// Created by xiang on 22-12-7.
//

#include "keyframe.h"
#include <glog/logging.h>
#include <pcl/io/pcd_io.h>
#include <iomanip>
#include "common/point_cloud_utils.h"

namespace sad {

/**
 * @description: 将当前关键帧的点云保存到本地，并从内存中清除
 * @param {string} &path
 * @return {*}
 */
void Keyframe::SaveAndUnloadScan(const std::string &path) {
    if (cloud_) {
        sad::SaveCloudToFile(path + "/" + std::to_string(id_) + ".pcd", *cloud_);
        cloud_ = nullptr;   // 保存后清除内存，指针置空
    }
}

/**
 * @description: 从pcd文件中加载关键帧点云数据
 * @param {string} &path
 * @return {*}
 */
void Keyframe::LoadScan(const std::string &path) {
    cloud_.reset(new PointCloudType);
    pcl::io::loadPCDFile(path + "/" + std::to_string(id_) + ".pcd", *cloud_);
}

/**
 * @description: 将关键帧数据保存到文件流中
 * @param {ostream} &os
 * @return {*}
 */
void Keyframe::Save(std::ostream &os) {
    // lambda表达式，保存7维数据到文件流中：平移3维，旋转四元数4维，返回SE3位姿
    auto save_SE3 = [](std::ostream &f, SE3 pose) {
                        auto q = pose.so3().unit_quaternion();  // SO3李群转单位四元数（实部在后）
                        Vec3d t = pose.translation();           // 平移向量
                        f << t[0] << " " << t[1] << " " << t[2] << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << " ";
                    };
    // 先保存该行数据的前5个：关键帧id，时间戳，RTK头部有效标记rtk_heading_valid，RTK数据有效位rtk_valid，RTK内点标记位rtk_inlier
    os << id_ << " " << std::setprecision(18) << timestamp_ << " " 
       << rtk_heading_valid_ << " " << rtk_valid_ << " " << rtk_inlier_ << " ";
    save_SE3(os, lidar_pose_);  // 保存雷达位姿
    save_SE3(os, rtk_pose_);    // 保存RTK位姿
    save_SE3(os, opti_pose_1_); // 保存第一轮优化位姿
    save_SE3(os, opti_pose_2_); // 保存第二轮优化位姿
    os << std::endl;
}

/**
 * @description: 从keyframes.txt文件中的一行读取关键帧中的信息
 * @param {istream} &is
 * @return {*}
 */
void Keyframe::Load(std::istream &is) {
    // 先读取该行数据的前5个：关键帧id，时间戳，RTK有效标记rtk_heading_valid，rtk_valid，rtk_inlier
    is >> id_ >> timestamp_ >> rtk_heading_valid_ >> rtk_valid_ >> rtk_inlier_;

    // lambda表达式，从文件流中读取7维数据：平移3维，旋转四元数4维，返回SE3位姿
    auto load_SE3 = [](std::istream &f) -> SE3 {
                        SE3 ret;
                        double q[4];
                        double t[3];
                        f >> t[0] >> t[1] >> t[2] >> q[0] >> q[1] >> q[2] >> q[3];
                        return SE3(Quatd(q[3], q[0], q[1], q[2]), Vec3d(t[0], t[1], t[2]));
                    };
    lidar_pose_ = load_SE3(is); // 读取雷达位姿
    rtk_pose_ = load_SE3(is);   // 读取RTK位姿
    opti_pose_1_ = load_SE3(is);// 第一轮优化位姿
    opti_pose_2_ = load_SE3(is);// 第二轮优化位姿
}

/**
 * @description: 从文件中读取关键帧
 * @param {string} &path  ./data/ch9/keyframes.txt
 * @param {map<IdType, std::shared_ptr<Keyframe>>} &keyframes
 * @return {*}
 */
bool LoadKeyFrames(const std::string &path, std::map<IdType, std::shared_ptr<Keyframe>> &keyframes) {
    std::ifstream fin(path);    
    if (!fin) // 判断文件打开是否失败
        return false;

    // 从文件中读取一行，直到文件末尾
    while (!fin.eof()) {
        std::string line;
        std::getline(fin, line);    // 读取一行数据

        // 如果读取到空行，说明文件末尾，退出循环
        if (line.empty()) 
            break;


        std::stringstream ss;
        ss << line; // 将一行数据转换为字符串流
        // 创建一个关键帧对象，从字符串流中读取数据
        auto kf = std::make_shared<Keyframe>();
        // 从字符串流中读取关键帧数据：id，时间戳，RTK有效标记，RTK数据有效位，RTK内点标记位，
        //                           雷达位姿，RTK位姿，第一轮优化位姿，第二轮优化位姿
        kf->Load(ss);   
        // 将关键帧对象加入到map中，以关键帧id为索引
        keyframes.emplace(kf->id_, kf);
    }

    LOG(INFO) << "Loaded kfs: " << keyframes.size();
    return true;
}
}  // namespace sad
