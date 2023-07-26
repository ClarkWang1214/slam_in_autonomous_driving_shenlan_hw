//
// Created by xiang on 22-12-8.
//

#include "optimization.h"
#include "common/math_utils.h"

#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/robust_kernel.h>
#include <g2o/core/sparse_block_matrix.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <glog/logging.h>
#include <yaml-cpp/yaml.h>
#include <boost/format.hpp>

namespace sad {

/// 打印优化信息
template <typename T>
std::string print_info(const std::vector<T>& edges, double th) {
    std::vector<double> chi2;
    // 遍历所有约束边
    for (auto& edge : edges) {
        // 判断边的优化级别是否为0， 1表示不参与优化
        if (edge->level() == 0) {
            // 对参与优化的边计算残差卡方值
            edge->computeError();
            chi2.push_back(edge->chi2());   // 将当前边的卡方值保存
        }
    }

    // 对残差卡方值进行排序
    std::sort(chi2.begin(), chi2.end());
    // 计算平均卡方值
    double ave_chi2 = std::accumulate(chi2.begin(), chi2.end(), 0.0) / chi2.size();
    // 格式化输出
    boost::format fmt("数量: %d, 均值: %f, 中位数: %f, 0.1分位: %f, 0.9分位: %f, 0.95分位：%f, 最大值: %f, 阈值: %f\n");
    // 判断卡方值集合是否为空
    if (!chi2.empty()) {
        // 返回格式化后的字符串
        std::string str = ( fmt % 
                            chi2.size() %                   // 数量
                            ave_chi2 %                      // 均值
                            chi2[chi2.size() / 2] %         // 中位数
                            chi2[int(chi2.size() * 0.1)] %  // 0.1分位
                            chi2[int(chi2.size() * 0.9)] %  // 0.9分位
                            chi2[int(chi2.size() * 0.95)] % // 0.95分位
                            chi2.back() %                   // 最大值
                            th)                             // 阈值
                            .str();
        return str;
    }
    return std::string();
}

/**
 * @description: 构造函数
 * @param {string&} yaml
 * @return {*}
 */
Optimization::Optimization(const std::string& yaml) { yaml_ = yaml; }

/**
 * @description: 初始化
 * @param {int} stage
 * @return {*}
 */
bool Optimization::Init(int stage) {
    stage_ = stage;

    // 加载关键帧信息
    if (!LoadKeyFrames("./data/ch9/keyframes.txt", keyframes_)) {
        LOG(ERROR) << "cannot load keyframes.txt";
        return false;
    }

    LOG(INFO) << "keyframes: " << keyframes_.size();

    // 读参数 mapping.yaml
    auto yaml = YAML::LoadFile(yaml_); 
    rtk_outlier_th_ = yaml["rtk_outlier_th"].as<double>();          // 获取RTK的outlier阈值
    lidar_continuous_num_ = yaml["lidar_continuous_num"].as<int>(); // 获取雷达连续帧数
    rtk_has_rot_ = yaml["rtk_has_rot"].as<bool>();                  // RTK是否有旋转信息

    rtk_pos_noise_ = yaml["rtk_pos_noise"].as<double>();                    // RTK位置噪声
    rtk_ang_noise_ = yaml["rtk_ang_noise"].as<double>() * math::kDEG2RAD;   // RTK旋转噪声
    rtk_height_noise_ratio_ = yaml["rtk_height_noise_ratio"].as<double>();  // RTK高度噪声比例

    std::vector<double> rtk_ext_t = yaml["rtk_ext"]["t"].as<std::vector<double>>(); // RTK外参
    TBG_ = SE3(SO3(), Vec3d(rtk_ext_t[0], rtk_ext_t[1], rtk_ext_t[2]));     // uclk数据集，单天线方案，仅有位移无旋转
    LOG(INFO) << "TBG = \n" << TBG_.matrix();

    if (stage_ == 2) 
        // 若是第二轮，则加载回环候选边
        LoadLoopCandidates();
    
    return true;
}

/**
 * @description: 开始后端位姿图优化g2o
 * @param {*}
 * @return {*}
 */
void Optimization::Run() {
    LOG(INFO) << "running optimization on stage " << stage_; // 输出所在优化阶段
    if (!rtk_has_rot_ && stage_ == 1) 
        // 【1】如果是第一轮优化，且RTK没有旋转，则进行RTK和雷达的初始对齐，基于SVD分解的ICP
        InitialAlign(); 

    // 【2】构建优化问题
    BuildProblem();  

    // 保存优化前的位姿图
    SaveG2O("./data/ch9/before.g2o");
    LOG(INFO) << "RTK 误差：" << print_info(gnss_edge_, rtk_outlier_th_);
    LOG(INFO) << "RTK 平移误差：" << print_info(gnss_trans_edge_, rtk_outlier_th_);
    LOG(INFO) << "lidar 误差：" << print_info(lidar_edge_, 0);

    // 【3】带着RK鲁棒核函数求解第一遍，迭代100次
    Solve();           
    // 【4】移除outlier异常值
    RemoveOutliers();  
    // 【5】再求解一遍，迭代100次
    Solve();           

    // 保存优化后的位姿图
    SaveG2O("./data/ch9/after.g2o");

    // 保存结果到keyframes.txt中
    SaveResults();  
    LOG(INFO) << "done";
}

/**
 * @description: 保存位姿结果到.g2o文件中
 * @param {string&} file_name
 * @return {*}
 */
void Optimization::SaveG2O(const std::string& file_name) {
    std::ofstream fout(file_name);
    // 遍历所有位姿顶点
    for (auto& v : vertices_) 
        v.second->write(fout);
    // 遍历所有雷达约束边 EdgeRelativeMotion
    for (auto& e : lidar_edge_) 
        e->write(fout);
    // 遍历所有回环约束边 EdgeRelativeMotion
    for (auto& e : loop_edge_) 
        e->write(fout);
    
    fout.close();
}

/**
 * @description: 构建优化问题
 * @return {*}
 */
void Optimization::BuildProblem() {
    using BlockSolverType = g2o::BlockSolverX;
    using LinearSolverType = g2o::LinearSolverEigen<BlockSolverType::PoseMatrixType>;
    auto* solver = new g2o::OptimizationAlgorithmLevenberg(
        g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));

    optimizer_.setAlgorithm(solver);

    AddVertices();      // 添加顶点
    AddRTKEdges();      // 添加RTK边
    AddLidarEdges();    // 添加lidar边
    AddLoopEdges();     // 添加loop边
}

/**
 * @description: 添加顶点
 * @return {*}
 */
void Optimization::AddVertices() {
    // 遍历所有关键帧
    for (auto& kfp : keyframes_) {
        auto kf = kfp.second;

        // make g2o vertex for this kf
        // 创建位姿顶点
        auto v = new VertexPose();
        v->setId(kf->id_);  // 位姿顶点的索引设为关键帧id

        // 判断是第一轮优化，还是第二轮优化
        if (stage_ == 1) 
            v->setEstimate(kf->lidar_pose_);    // 设置第一轮的估计值为雷达位姿
        else 
            v->setEstimate(kf->opti_pose_1_);   // 设置第二轮的估计值为第一轮优化位姿
        
        optimizer_.addVertex(v);        // 将当前关键帧的位姿顶点添加到图优化器中
        vertices_.emplace(kf->id_, v);  // 添加位姿顶点到map容器中，按关键帧id排序
    }
    LOG(INFO) << "vertex: " << vertices_.size();
}

/**
 * @description: 添加RTK边
 * @return {*}
 */
void Optimization::AddRTKEdges() {
    /// RTK 噪声设置
    Mat3d info_pos = Mat3d::Identity() * 1.0 / (rtk_pos_noise_ * rtk_pos_noise_);   // 位移噪声协方差的倒数
    info_pos(2, 2) = 1.0 / (rtk_height_noise_ratio_ * rtk_pos_noise_ * rtk_height_noise_ratio_ * rtk_pos_noise_);
    Mat3d info_ang = Mat3d::Identity() * 1.0 / (rtk_ang_noise_ * rtk_ang_noise_);   // 角度噪声协方差的倒数
    Mat6d info_all = Mat6d::Identity();

    // info_all.block<3, 3>(0, 0) = info_pos; // 写反了，角度在前，平移在后！！！
    // info_all.block<3, 3>(3, 3) = info_ang;
    info_all.block<3, 3>(0, 0) = info_ang; 
    info_all.block<3, 3>(3, 3) = info_pos;

    LOG(INFO) << "Info of rtk trans: " << info_pos.diagonal().transpose();

    // 是否为第2轮优化
    if (stage_ == 2) {
        // 若是，信息矩阵权重再减小一点
        info_pos *= 0.01;
        info_all *= 0.01;
    }

    for (auto& kfp : keyframes_) {
        auto kf = kfp.second;

        // RTK数据是否有效
        if (!kf->rtk_valid_) 
            continue; // 跳过无效RTK数据

        // RTK航向角是否有效
        if (kf->rtk_heading_valid_) {
            // 若有效，说明是双天线方案，既有旋转，也有平移信息
            auto edge = new EdgeGNSS(vertices_.at(kf->id_), kf->rtk_pose_);  // 6自由度RTK约束边
            edge->setInformation(info_all);         // 设置信息矩阵 6x6
            auto rk = new g2o::RobustKernelHuber(); // Huber鲁棒核函数
            rk->setDelta(rtk_outlier_th_);          // 设置outlier外点判定阈值
            edge->setRobustKernel(rk);              // 给约束边设置鲁棒核
            optimizer_.addEdge(edge);
            gnss_edge_.emplace_back(edge);
        } else {
            // 否则，说明是单天线方案，仅有平移，无旋转信息，需给定RTK外参
            auto edge = new EdgeGNSSTransOnly(vertices_.at(kf->id_), kf->rtk_pose_.translation(), TBG_); // 3自由度约束边
            edge->setInformation(info_pos);             // 设置位置信息矩阵 3x3
            auto rk = new g2o::RobustKernelCauchy();    // Cauchy鲁棒核函数
            rk->setDelta(rtk_outlier_th_);              // 设置outlier外点判定阈值
            edge->setRobustKernel(rk);                  // 给约束边设置鲁棒核
            optimizer_.addEdge(edge);   
            gnss_trans_edge_.emplace_back(edge);
        }
    }

    LOG(INFO) << "gnss edges: " << gnss_edge_.size() << ", " << gnss_trans_edge_.size();
}

/**
 * @description: 添加lidar边
 * @return {*}
 */
void Optimization::AddLidarEdges() {
    const double lidar_pos_noise = 0.01, 
                 lidar_ang_noise = 0.1 * math::kDEG2RAD;  // RTK 观测的噪声
    Mat3d info_pos = Mat3d::Identity() * 1.0 / (lidar_pos_noise * lidar_pos_noise);
    Mat3d info_ang = Mat3d::Identity() * 1.0 / (lidar_ang_noise * lidar_ang_noise);
    Mat6d info_all = Mat6d::Identity();
    info_all.block<3, 3>(0, 0) = info_pos;
    info_all.block<3, 3>(3, 3) = info_ang;

    for (auto iter = keyframes_.begin(); iter != keyframes_.end(); ++iter) {
        auto iter_next = iter;
        for (int i = 0; i < lidar_continuous_num_; ++i) {
            iter_next++;

            if (iter_next == keyframes_.end()) 
                break;

            // 添加iter和iter_next之间的相邻运动
            auto edge = new EdgeRelativeMotion( vertices_.at(iter->second->id_), 
                                                vertices_.at(iter_next->second->id_),
                                                iter->second->lidar_pose_.inverse() * iter_next->second->lidar_pose_);
            edge->setInformation(info_all);
            optimizer_.addEdge(edge);
            lidar_edge_.emplace_back(edge);
        }
    }

    LOG(INFO) << "lidar edges: " << lidar_edge_.size();
}

/**
 * @description:  添加回环边
 * @return {*}
 */
void Optimization::AddLoopEdges() {
    // 第一轮优化不添加loop回环边
    if (stage_ == 1) 
        return; 

    const double loop_pos_noise = 0.1,                   // 回环位移噪声
                 loop_ang_noise = 0.5 * math::kDEG2RAD;  // 回环角度噪声
    // 位置和角度的信息矩阵
    Mat3d info_pos = Mat3d::Identity() * 1.0 / (loop_pos_noise * loop_pos_noise);
    Mat3d info_ang = Mat3d::Identity() * 1.0 / (loop_ang_noise * loop_ang_noise);
    Mat6d info_all = Mat6d::Identity();
    info_all.block<3, 3>(0, 0) = info_pos;
    info_all.block<3, 3>(3, 3) = info_ang;

    const double loop_rk_th = 5.2; // 鲁棒核函数的阈值

    // 遍历所有回环候选边
    for (const auto& lc : loop_candidates_) {
        auto edge = new EdgeRelativeMotion(vertices_.at(lc.idx1_), 
                                           vertices_.at(lc.idx2_), 
                                           lc.Tij_);
        edge->setInformation(info_all); // 设置信息矩阵6x6

        auto rk = new g2o::RobustKernelCauchy();    // 创建鲁棒Cauchy核函数
        rk->setDelta(loop_rk_th);                   // 设置阈值
        edge->setRobustKernel(rk);                  // 设置鲁棒核函数
        optimizer_.addEdge(edge);                   // 添加边
        loop_edge_.emplace_back(edge);              // 添加到回环边容器中
    }
}

/**
 * @description: 求解优化问题
 * @return {*}
 */
void Optimization::Solve() {
    optimizer_.setVerbose(true);            // 设置输出优化信息
    optimizer_.initializeOptimization(0);   // 初始化优化器
    optimizer_.optimize(100);               // 进行100次优化

    LOG(INFO) << "RTK 误差：" << print_info(gnss_edge_, rtk_outlier_th_);
    LOG(INFO) << "RTK 平移误差：" << print_info(gnss_trans_edge_, rtk_outlier_th_);
    LOG(INFO) << "lidar 误差：" << print_info(lidar_edge_, 0);
    LOG(INFO) << "loop 误差：" << print_info(loop_edge_, 0);
}

/**
 * @description: 移除异常值
 * @return {*}
 */
void Optimization::RemoveOutliers() {
    // 主要用于移除GNSS的异常值
    int cnt_outlier_removed = 0;
    // lambda表达式，用于移除异常值
    auto remove_outlier = [&cnt_outlier_removed](g2o::OptimizableGraph::Edge* e) {
                            // 判断残差的卡方值: _error.dot(information()*_error) 是否大于鲁棒核函数设置的阈值delta
                            if (e->chi2() > e->robustKernel()->delta()) {
                                // 若超过了，说明是outlier外点
                                e->setLevel(1); // 将该边的优化级别设为1，即不参与优化
                                cnt_outlier_removed++;  // 外点计数加一
                            } else 
                                // 否则，说明第一轮中该边为inlier内点
                                // 此时，去掉鲁棒核中的delta值
                                e->setRobustKernel(nullptr); 
                        };

    // 遍历每个RTK位姿约束边
    std::for_each(gnss_edge_.begin(), gnss_edge_.end(), remove_outlier);

    // 遍历每个RTK仅位移约束边
    std::for_each(gnss_trans_edge_.begin(), gnss_trans_edge_.end(), remove_outlier);
    LOG(INFO) << "gnss outlier: " << cnt_outlier_removed << "/" << gnss_edge_.size() + gnss_trans_edge_.size();

    cnt_outlier_removed = 0;
    std::for_each(loop_edge_.begin(), loop_edge_.end(), remove_outlier);
    LOG(INFO) << "loop outlier: " << cnt_outlier_removed << "/" << loop_edge_.size();
}

/**
 * @description: 保存优化结果
 * @return {*}
 */
void Optimization::SaveResults() {
    for (auto& v : vertices_) {
        // 根据优化的阶段，选择不同的位姿
        if (stage_ == 1) 
            // 保存第一轮优化的估计位姿
            keyframes_.at(v.first)->opti_pose_1_ = v.second->estimate();
        else 
            // 保存第二轮优化的估计位姿
            keyframes_.at(v.first)->opti_pose_2_ = v.second->estimate();
    }

    // 比较优化pose和rtk pose
    std::vector<double> rtk_trans_error; // RTK位移误差
    // 遍历所有关键帧
    for (auto& kfp : keyframes_) {
        auto kf = kfp.second;   // 当前关键帧
        Vec3d tWG = kf->rtk_pose_.translation();  // RTK位姿的平移部分
        Vec3d t_opti = (kf->opti_pose_1_ * TBG_).translation(); // 优化位姿的平移部分
        double n = (tWG - t_opti).head<2>().norm(); // 计算RTK位移与优化位移之间距离前两维的模长
        rtk_trans_error.emplace_back(n);
    }

    // 对RTK位移误差进行排序
    std::sort(rtk_trans_error.begin(), rtk_trans_error.end());
    LOG(INFO) << "med error: " << rtk_trans_error[rtk_trans_error.size() / 2];

    // 写入文件
    system("rm ./data/ch9/keyframes.txt");
    std::ofstream fout("./data/ch9/keyframes.txt");
    for (auto& kfp : keyframes_) 
        kfp.second->Save(fout);
    
    fout.close();
}

/**
 * @description: 如果RTK全程不带旋转，那么先对激光和RTK做一次ICP来对齐整条轨迹
 * @return {*}
 */
void Optimization::InitialAlign() {
    // should be p1 = R*p2 + t
    std::vector<Vec3d> pts1, pts2;
    // 遍历所有关键帧
    for (auto& kfp : keyframes_) {
        // 分别存储关键帧的RTK位姿和雷达位姿中的平移部分到pts1和pts2中
        pts1.emplace_back(kfp.second->rtk_pose_.translation());  // RTK点坐标
        pts2.emplace_back(kfp.second->lidar_pose_.translation());// 雷达点坐标
    }

    // 计算两组点的质心
    Vec3d p1, p2;  // center of mass
    int N = pts1.size();
    for (int i = 0; i < N; i++) {
        p1 += pts1[i];
        p2 += pts2[i];
    }
    p1 = p1 / N;
    p2 = p2 / N;

    LOG(INFO) << "p1: " << p1.transpose() << ", p2: " << p2.transpose();

    //  remove the center 得到去质心坐标
    std::vector<Vec3d> q1(N), q2(N);  
    for (int i = 0; i < N; i++) {
        q1[i] = pts1[i] - p1;
        q2[i] = pts2[i] - p2;
    }

    // compute q1*q2^T
    // 计算雷达与RTK位移之间的
    Mat3d W = Mat3d::Zero(); // 3x3的零矩阵
    for (int i = 0; i < N; i++) 
        W += q1[i] * q2[i].transpose(); // 对应十四讲公式（7.57）

    // SVD on W
    // 对W矩阵进行SVD分解
    // 详见十四讲7.9.1节中基于SVD方法的ICP算法，公式（7.58）
    Eigen::JacobiSVD<Mat3d> svd(W, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Mat3d U = svd.matrixU();
    Mat3d V = svd.matrixV();
    Mat3d R = U * (V.transpose());  // 得到旋转矩阵，对应十四讲公式（7.59）则取-R作为最优值
    
    if (R.determinant()<0)
        R = -R;

    // 再根据公式（7.54），利用计算得到的R来计算平移矩阵t
    Vec3d t = p1 - R * p2;

    // change lidar pose 改变雷达位姿
    SE3 T(R, t);
    LOG(INFO) << "initial trans: \n" << T.matrix();
    for (auto& kfp : keyframes_) 
        kfp.second->lidar_pose_ = T * kfp.second->lidar_pose_;
}

/**
 * @description: 加载回环候选边
 * @return {*}
 */
void Optimization::LoadLoopCandidates() {
    std::ifstream fin("./data/ch9/loops.txt");
    if (!fin) {
        LOG(WARNING) << "cannot load file: ./data/ch9/loops.txt";
        return;
    }

    auto load_SE3 = [](std::istream& f) -> SE3 {
        SE3 ret;
        double q[4];
        double t[3];
        f >> t[0] >> t[1] >> t[2] >> q[0] >> q[1] >> q[2] >> q[3];
        return SE3(Quatd(q[3], q[0], q[1], q[2]), Vec3d(t[0], t[1], t[2]));
    };

    while (fin.eof() == false) {
        std::string line;
        std::getline(fin, line);
        if (line.empty()) 
            break;

        std::stringstream ss;
        ss << line;

        LoopCandidate lc;
        ss >> lc.idx1_ >> lc.idx2_ >> lc.ndt_score_;
        lc.Tij_ = load_SE3(ss);
        loop_candidates_.emplace_back(lc);
    }

    LOG(INFO) << "loaded loops: " << loop_candidates_.size();
}

}  // namespace sad