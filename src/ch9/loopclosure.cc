//
// Created by xiang on 22-12-9.
//

#include "loopclosure.h"

#include <glog/logging.h>
#include <pcl/common/transforms.h>
#include <pcl/io/pcd_io.h>
#include <pcl/registration/ndt.h>
#include <yaml-cpp/yaml.h>
#include <execution>

#include "common/lidar_utils.h"
#include "common/point_cloud_utils.h"

namespace sad {

LoopClosure::LoopClosure(const std::string& config_yaml) : yaml_(config_yaml) {}

/**
 * @description: 初始化回环检测
 * @return {*}
 */
bool LoopClosure::Init() {
    if (!LoadKeyFrames("./data/ch9/keyframes.txt", keyframes_)) {
        LOG(ERROR) << "cannot load keyframes";
        return false;
    }
    LOG(INFO) << "keyframes: " << keyframes_.size();

    // 从yaml文件中加载回环检测相关的参数值
    auto yaml = YAML::LoadFile(yaml_);
    min_id_interval_ = yaml["loop_closing"]["min_id_interval"].as<int>();   // 最小时间间隔
    min_distance_ = yaml["loop_closing"]["min_distance"].as<double>();      // 最小距离
    skip_id_ = yaml["loop_closing"]["skip_id"].as<int>();                   // ID间隔（默认为5个），每隔多少个关键帧做一次回环检测
    ndt_score_th_ = yaml["loop_closing"]["ndt_score_th"].as<double>();      // NDT分值阈值，作为回环有效性的判定依据
    use_pcl_ndt_ = yaml["use_pcl_ndt"].as<bool>();          // 是否使用PCL NDT库
    return true;
}

/**
 * @description: 回环检测算法是简单的先检测，后计算的方式
 * @return {*}
 */
void LoopClosure::Run() {
    // 检测回环
    DetectLoopCandidates();
    
    // 计算回环
    ComputeLoopCandidates();

    // 保存结果
    SaveResults();
}

/**
 * @description: 检测出候选关键帧对
 *          每两个在空间上相隔较近，但时间上存在一定距离的关键帧，进行一次回环检测
 *          每隔5个关键帧做一次检查点关键帧的更新判定
 * @return {*}
 */
void LoopClosure::DetectLoopCandidates() {
    // 每两个在空间上相隔较近，但时间上存在一定距离的关键帧，进行一次回环检测，
    // 我们称这样的一对关键帧为检查点
    KFPtr check_first = nullptr; 
    KFPtr check_second = nullptr;

    // 从第一轮优化中的位姿中检测回环候选
    LOG(INFO) << "detecting loop candidates from pose in stage 1";

    // 本质上是两重循环
    // 第一次遍历第一轮优化轨迹中的关键帧。
    // 对每两个在空间上相隔较近，但时间是存在一定距离的关键帧进行一次回环检测
    for (auto iter_first = keyframes_.begin(); iter_first != keyframes_.end(); ++iter_first) {
        auto kf_first = iter_first->second; // 获取第一个关键帧

        // 判定第一个检查点关键帧是否非空，且，与第一个非空检查关键帧之间的ID差值是否小于等于ID间隔阈值（默认间隔5个关键帧）
        if (check_first != nullptr && abs(int(kf_first->id_) - int(check_first->id_)) <= skip_id_) 
            // 两个关键帧之间ID相差太近，跳过当前循环，直到间隔大于5继续
            continue; 

        // 第二次遍历关键帧，从上面的第一个关键帧开始，遍历到最后
        for (auto iter_second = iter_first; iter_second != keyframes_.end(); ++iter_second) {
            auto kf_second = iter_second->second;

            // 判定第二个检查点关键帧是否非空，且，与第二个非空检查关键帧之间的ID差值是否小于等于ID间隔阈值（默认间隔5个关键帧）
            if (check_second != nullptr && abs(int(kf_second->id_) - int(check_second->id_)) <= skip_id_) 
                // 两个关键帧之前ID太近
                continue;
            
            // 判定第一个关键帧与第二个关键帧之间的id之差的绝对值是否小于最小间隔50
            if (abs(int(kf_first->id_) - int(kf_second->id_)) < min_id_interval_) 
                /// 在同一条轨迹中，如果间隔太近（小于50个关键帧），就不考虑回环
                continue;
            
            // 此时，第一个关键帧与第二个关键帧之间的id索引相差50个

            // 计算这两个关键帧之间第一轮优化后的位姿之间的位移差
            Vec3d dt = kf_first->opti_pose_1_.translation() - kf_second->opti_pose_1_.translation();
            double t2d = dt.head<2>().norm();   // 计算xy平面上的距离
            double range_th = min_distance_;    // 距离阈值设置为最小30米

            // 如果两个关键帧之间的距离小于30米，就认为是一个候选关键帧对
            if (t2d < range_th) {
                LoopCandidate c(kf_first->id_,  // 第一个关键帧的id
                                kf_second->id_, // 第二个关键帧的id
                                kf_first->opti_pose_1_.inverse() * kf_second->opti_pose_1_);    // 位姿差
                loop_candiates_.emplace_back(c);    // 将【候选关键帧对】保存到loop_candiates_中
                check_first = kf_first;     // 将第一个关键帧赋值给check_first，作为第一个检查点关键帧
                check_second = kf_second;   // 将第二个关键帧赋值给check_second，作为第二个检查点关键帧
            }
        }
    }
    LOG(INFO) << "detected candidates: " << loop_candiates_.size();
}

/**
 * @description: 对每个关键帧对，使用scan to map的配准方式，在关键帧附近抽取一定数量的点云作为子地图
 *               然后对scan和子地图进行配准。这种方式可以避免单个扫描数据量太少，点云纹理不充分的问题，缺点是计算量较大。
 *               由于每个检查点的计算都是独立的，使用并发编程让上面步骤能在机器上并发执行，加快计算速度。
 * 
 *               配准的实际执行由NDT完成，采用NDT分值作为回环有效性的判定依据。与里程计方法不同，在回环检测配准过程中，
 *               经常要面对初值位姿估计很差的情况，希望不太依赖于给定的位姿初值。因此，给NDT方法增加由粗到精的配准过程，也就是多分辨率NDT
 * @return {*}
 */
void LoopClosure::ComputeLoopCandidates() {
    // 执行计算
    // 并发遍历所有候选关键帧对，对每对关键帧计算
    std::for_each(std::execution::par_unseq, 
                    loop_candiates_.begin(), loop_candiates_.end(),
                    [this](LoopCandidate& c) { 
                        ComputeForCandidate(c); 
                    });

    // 保存成功的回环候选关键帧对
    std::vector<LoopCandidate> succ_candidates;

    // 选择最优的匹配结果
    auto lc_max_score = std::max_element(loop_candiates_.begin(), loop_candiates_.end(),
                                    [](const auto& lc1, const auto& lc2) { return lc1.ndt_score_ < lc2.ndt_score_; });

    auto lc_min_score = std::min_element(loop_candiates_.begin(), loop_candiates_.end(),
                                    [](const auto& lc1, const auto& lc2) { return lc1.ndt_score_ < lc2.ndt_score_; });

    LOG(INFO) << "max_score: " << lc_max_score->ndt_score_;  
    LOG(INFO) << "min_score: " << lc_min_score->ndt_score_;   

    // 遍历所有候选关键帧对
    for (const auto& lc : loop_candiates_) {
        // 计算每个候选关键帧对的NDT分值，判断NDT匹配的好坏
        // 参考PCL NDT的transprobability值，设计一个匹配度评估指标
        if (lc.ndt_score_ > ndt_score_th_) 
            succ_candidates.emplace_back(lc);
        // if(use_pcl_ndt_) {
        //     if (lc.ndt_score_ > ndt_score_th_) 
        //         succ_candidates.emplace_back(lc);
        // } 
        // else {
        //     if (lc.ndt_score_ < ndt_score_th_) 
        //         succ_candidates.emplace_back(lc);
        // }
    }
    LOG(INFO) << "success: " << succ_candidates.size() << "/" << loop_candiates_.size();

    // 交换两个vector容器中的内容
    loop_candiates_.swap(succ_candidates); // 此时，loop_candiates_中保存的是所有成功的回环候选关键帧对
}

/**
 * @description: 
 * @param {LoopCandidate&} c
 * @return {*}
 */
void LoopClosure::ComputeForCandidate(sad::LoopCandidate& c) {
    // 每一对候选关键帧之间进行NDT配准
    // LOG(INFO) << "aligning " << c.idx1_ << " with " << c.idx2_;
    const int submap_idx_range = 40;

    // 获取候选关键帧对中的这两个关键帧
    KFPtr kf1 = keyframes_.at(c.idx1_), 
          kf2 = keyframes_.at(c.idx2_);

    // lambda表达式，用于构建当前候选关键帧附近的子地图点云
    auto build_submap = [this](int given_id, bool build_in_world) -> CloudPtr {
                            // 用关键帧周围的点云构成子地图
                            CloudPtr submap(new PointCloudType); 
                            // 在当前关键帧索引前后各取40/4=10个关键帧，构成子地图
                            // idx - 40, idx - 36, idx - 32, ..., idx - 4, idx, idx + 4, ..., idx + 32, idx + 36, idx + 40
                            for (int idx = -submap_idx_range; idx < submap_idx_range; idx += 4) {
                                int id = idx + given_id; 
                                if (id < 0) continue; // 跳过负索引
                                
                                auto iter = keyframes_.find(id);
                                if (iter == keyframes_.end())  continue; // 跳过不存在的关键帧

                                auto kf = iter->second; // 获取找到的关键帧

                                CloudPtr cloud(new PointCloudType); 
                                // 读取以索引命名当前关键帧点云，存于cloud中
                                pcl::io::loadPCDFile("./data/ch9/" + std::to_string(id) + ".pcd", *cloud);

                                // 剔除地面点云，保留高度大于1.0的点云即可
                                sad::RemoveGround(cloud, 0.1);

                                if (cloud->empty())  continue; // 若非地面点云为空，跳过

                                // 转到世界系下
                                SE3 Twb = kf->opti_pose_1_; // 获取第一轮优化得到的位姿

                                // 是否在世界系下构建子地图
                                if (!build_in_world) 
                                    // 若为false，则转换到当前关键帧kf1所在坐标系下 
                                    //  Twb1^-1 * Twb = Tb1w * Twb = Tb1b 
                                    Twb = keyframes_.at(given_id)->opti_pose_1_.inverse() * Twb;

                                CloudPtr cloud_trans(new PointCloudType); // 用于存储转换后的点云
                                // 将当前关键帧周围的其它关键帧点云全部转换到 世界系 或者 当前关键帧所在坐标系下
                                pcl::transformPointCloud(*cloud, *cloud_trans, Twb.matrix());

                                // 将转换到世界系下的所有附近关键帧点云拼接在一起
                                *submap += *cloud_trans; 
                            }
                            return submap; // 返回当前关键帧周围的子地图点云
                        };

    // 用于构建当前候选关键帧附近的子地图点云，true表示在世界系下构建子地图
    auto submap_kf1 = build_submap(kf1->id_, true);

    CloudPtr submap_kf2(new PointCloudType);
    // 获取另一个候选关键帧点云，直接作为第二个子地图的点云
    pcl::io::loadPCDFile("./data/ch9/" + std::to_string(kf2->id_) + ".pcd", *submap_kf2);

    // 判断关键帧1周围的子地图 或者 关键帧2点云 是否为空
    if (submap_kf1->empty() || submap_kf2->empty()) {
        c.ndt_score_ = 0; // 若至少一个为空，则NDT得分为0，直接返回
        return;
    }

    // 创建一个NDT高斯分布变换对象
    pcl::NormalDistributionsTransform<PointType, PointType> ndt_pcl;

    ndt_pcl.setTransformationEpsilon(0.05);     
    ndt_pcl.setStepSize(0.7);   // 设置线搜索允许的最大步长
    ndt_pcl.setMaximumIterations(40);

    Ndt3d ndt;  // 【新增】自己实现的3D NDT

    // 第一轮优化得到的SE3位姿转换为Eigen::Matrix4f类型
    Mat4f Tw2 = kf2->opti_pose_1_.matrix().cast<float>();
    SE3 Tw2_se3 = kf2->opti_pose_1_;

    /// 不同分辨率下的匹配
    CloudPtr output(new PointCloudType);
    std::vector<double> res{10.0, 5.0, 4.0, 3.0}; // 分辨率由大到小
    // 多分辨率NDT，由粗到精的配准过程
    for (auto& r : res) {
        // 指定不同分辨率，对点云进行体素滤波voxel filter
        auto rough_map1 = VoxelCloud(submap_kf1, r * 0.1); // 滤波后的子地图1
        auto rough_map2 = VoxelCloud(submap_kf2, r * 0.1); // 滤波后的子地图2

        if(use_pcl_ndt_) {
            // LOG(INFO) << "pcl ndt: ";
            ndt_pcl.setResolution(r);   // 设置内部NDT网格结构的体素分辨率
            ndt_pcl.setInputTarget(rough_map1);
            ndt_pcl.setInputSource(rough_map2);
            ndt_pcl.align(*output, Tw2);
            Tw2 = ndt_pcl.getFinalTransformation();
        } else {
            // LOG(INFO) << "7th chapter ndt: ";
            ndt.SetResolution(r);       // 【新增】设置自己实现的3D NDT的分辨率
            ndt.SetTarget(rough_map1);  // 【新增】设置自己实现的3D NDT的目标点云
            ndt.SetSource(rough_map2);  // 【新增】设置自己实现的3D NDT的源点云
            ndt.AlignNdt(Tw2_se3);
        }        
    }

    if(use_pcl_ndt_) {
        Mat4d T = Tw2.cast<double>();
        Quatd q(T.block<3, 3>(0, 0));
        q.normalize();
        Vec3d t = T.block<3, 1>(0, 3);
        c.Tij_ = kf1->opti_pose_1_.inverse() * SE3(q, t);
        c.ndt_score_ = ndt_pcl.getTransformationProbability();
        LOG(INFO) << "ndt_pcl_score_: "<<ndt_pcl.getTransformationProbability();
    }
    else {
        Mat4d T = Tw2_se3.matrix().cast<double>();
        Quatd q(T.block<3, 3>(0, 0));
        q.normalize();
        Vec3d t = T.block<3, 1>(0, 3);
        c.Tij_ = kf1->opti_pose_1_.inverse() * SE3(q, t);
        c.ndt_score_ = ndt.GetNdtMatchingScore();  // 【新增】
        LOG(INFO) << "ndt_score_: "<<ndt.GetNdtMatchingScore();
    }
        
}

/**
 * @description: 保存所有成功的回环，在下次调用优化算法时读取回环检测的结果
 * @return {*}
 */
void LoopClosure::SaveResults() {
    auto save_SE3 = [](std::ostream& f, SE3 pose) {
        auto q = pose.so3().unit_quaternion();
        Vec3d t = pose.translation();
        f << t[0] << " " << t[1] << " " << t[2] << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << " ";
    };

    std::ofstream fout("./data/ch9/loops.txt");
    for (const auto& lc : loop_candiates_) {
        fout << lc.idx1_ << " " << lc.idx2_ << " " << lc.ndt_score_ << " ";
        // save_SE3(fout, lc.Tij_);
        auto q = lc.Tij_.so3().unit_quaternion();
        Vec3d t = lc.Tij_.translation();
        fout << t[0] << " " << t[1] << " " << t[2] << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << " ";
        fout << std::endl;
    }
}

}  // namespace sad