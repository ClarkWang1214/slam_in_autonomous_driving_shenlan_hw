//
// Created by xiang on 2022/3/15.
//

#include "ch6/icp_2d.h"
#include "common/math_utils.h"

#include <glog/logging.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/impl/kdtree.hpp>

namespace sad {

/**
 * @description: 基于高斯牛顿迭代的2D ICP算法
 * 
 *      点到点ICP方法中，在某个目标点与另一个点云中的最近邻之间计算欧式距离误差。
 *      这个误差会随着这两个点的距离平方增长，最后累加之后形成问题的目标函数。
 *      可以想象成每个点与它的最近邻之间安装了一个弹簧。这些弹簧产生的拉力最终会把
 *      点云拉至能量最小的位置上。然而，在ICP方法中，每迭代一次，就必须将这些弹簧重新置零安装一遍，比较耗时
 * 
 * @param {SE2&} init_pose 初始位姿
 * @return {*}
 */
bool Icp2d::AlignGaussNewton(SE2& init_pose) {
    int iterations = 10;                // 迭代次数
    double cost = 0, lastCost = 0;      // 代价函数，上一次的代价函数
    SE2 current_pose = init_pose;       // 当前位姿
    float max_dis2 = 0.01;        // 最近邻时的最远距离（平方）
    int min_effect_pts = 20;      // 最小有效点数

    // LOG(INFO) << "initial pose: " << init_pose.translation().transpose() << ", theta: " << init_pose.so2().log();

    // 迭代10次
    for (int iter = 0; iter < iterations; ++iter) {
        Mat3d H = Mat3d::Zero();    // 用于存储所有点累加后的 H = H + J * J^T
        Vec3d b = Vec3d::Zero();    // 用于存储所有点累加后的 b = b - J * e
        cost = 0;

        int effective_num = 0;  // 有效点数

        // 遍历源始点云
        for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
            double r = source_scan_->ranges[i]; // 源始点云的距离
            // 判断每个点的距离是否越界
            if (r < source_scan_->range_min || r > source_scan_->range_max) continue;
            
            // 根据最小角度和角分辨率计算每个点的角度
            double angle = source_scan_->angle_min + i * source_scan_->angle_increment;
            double theta = current_pose.so2().log(); // 当前位姿的角度
            // p_i^W = T_WB * p_i^B，将机器人坐标系下的点（极坐标转笛卡尔坐标）转换到世界坐标系下
            Vec2d pw = current_pose * Vec2d(r * std::cos(angle), r * std::sin(angle));

            // Vec2d translation = current_pose.translation(); // 当前位姿的位移
            // double x = translation.x();
            // double y = translation.y();
            // Vec2d pw2 = Vec2d(x + r * std::cos(angle+theta), y + r * std::sin(angle+theta));
            // LOG(INFO) << std::fixed << std::setprecision(6) << "pw = " << pw.transpose();
            // LOG(INFO) << std::fixed << std::setprecision(6) << "pw2 = " << pw2.transpose();

            Point2d pt;
            pt.x = pw.x();
            pt.y = pw.y();

            // 最近邻
            std::vector<int> nn_idx;    // 最近邻的索引
            std::vector<float> dis;     // 最近邻的距离
            // 在目标点云的KD树中查找一个最近邻，返回该最近邻的索引和距离
            kdtree_2d->nearestKSearch(pt, 1, nn_idx, dis);

            if (nn_idx.size() > 0 && dis[0] < max_dis2) {
                effective_num++;
                Mat32d J;   // 用于存储当前点的Jacobian矩阵，3x2
                // 残差公式（6.6）对于各状态变量（x,y,theta）的偏导数，对应公式（6.8）、（6.9）
                // 2D位姿的明显优势是不必再使用流形中的符号，直接使用拆散了的x,y,theta即可
                J << 1, 0,                                                      // d e / d x
                     0, 1,                                                      // d e / d y
                     -r * std::sin(angle + theta), r * std::cos(angle + theta); // d e / d theta

                // 高斯牛顿GN法的增量方程 H = J * J^T，对应《十四讲》公式（6.32）
                H += J * J.transpose();

                // 将当前激光点与目标点云的最近邻点作差，得到误差，对应公式（6.6）
                Vec2d e(pt.x - target_cloud_2d->points[nn_idx[0]].x, pt.y - target_cloud_2d->points[nn_idx[0]].y);
                b += -J * e; // 高斯牛顿GN法的增量方程 b = b - J * e，对应《十四讲》公式（6.32）

                cost += e.dot(e); // 计算残差的平方和，对应公式（6.7）
            }
        }

        // 判断有效点数是否小于最小有效点数
        if (effective_num < min_effect_pts) return false;
        
        // solve for dx：  H * dx = b   ->   dx = H.ldlt().solve(b)
        // 求解线性方程dx，对应公式（6.7）
        Vec3d dx = H.ldlt().solve(b);

        // 判断dx是否为nan值
        if (isnan(dx[0])) break;
        
        // 计算平均残差
        cost /= effective_num;

        // 判断迭代次数是否多于1次，且平均残差是否大于上一次的平均残差
        if (iter > 0 && cost >= lastCost)  break; // 若是，说明上一次迭代的平均残差已经是最小的了，不需再迭代，直接退出
        
        // LOG(INFO) << "iter " << iter << " cost = " << cost << ", effect num: " << effective_num;

        // 对应于VertexSE2中的oplusImpl()中广义加法运算的实现，用于位姿更新
        current_pose.translation() += dx.head<2>(); // 更新位姿的平移部分
        current_pose.so2() = current_pose.so2() * SO2::exp(dx[2]);  // 更新位姿的旋转部分
        LOG(INFO) << "current_pose: " << current_pose.translation().transpose() << ", theta: " << current_pose.so2().log();
        lastCost = cost;    // 更新上一次的平均残差
    }

    init_pose = current_pose;   // 将此时的位姿赋值作为下一次的初始位姿
    // LOG(INFO) << "estimated pose: " << current_pose.translation().transpose() << ", theta: " << current_pose.so2().log();

    return true;
}

/**
 * @description: 【新增】基于G2O优化器的2D 点到点ICP算法
 * @param {SE2&} init_pose 初始位姿
 * @return {*}
 */
bool Icp2d::AlignG2OP2P(SE2& init_pose) {
    int iterations = 10;                // 迭代次数
    double rk_delta = 0.8;
    float max_dis2 = 0.01;        // 最近邻时的最远距离（平方）
    int min_effect_pts = 20;      // 最小有效点数
    SE2 current_pose = init_pose;   // 当前位姿
    for (int iter = 0; iter < iterations; ++iter) {
        using BlockSolverType = g2o::BlockSolver<g2o::BlockSolverTraits<3, 1>>;
        using LinearSolverType = g2o::LinearSolverCholmod<BlockSolverType::PoseMatrixType>;
        auto* solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
        g2o::SparseOptimizer optimizer;
        optimizer.setAlgorithm(solver);
        
        auto* v = new VertexSE2();      // 新建SE2位姿顶点
        v->setId(0);                    // 设置顶点的id
        v->setEstimate(current_pose);   // 设置顶点的估计值为初始位姿
        optimizer.addVertex(v);         // 将顶点添加到优化器中
        int effective_num = 0;  // 有效点数
        // 遍历源始点云
        for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
            double range = source_scan_->ranges[i]; // 源始点云的距离
            // 判断每个点的距离是否越界
            if (range < source_scan_->range_min || range > source_scan_->range_max) 
                continue;

            // 根据最小角度和角分辨率计算每个点的角度
            double angle = source_scan_->angle_min + i * source_scan_->angle_increment;
            double theta = current_pose.so2().log(); // 当前位姿的角度

            // 世界系下点的坐标 p_i^W，极坐标转笛卡尔坐标公式
            Vec2d pw = current_pose * Vec2d(range * std::cos(angle), range * std::sin(angle));

            Point2d pt;
            pt.x = pw.x();
            pt.y = pw.y();

            // 最近邻
            std::vector<int> nn_idx;    // 最近邻的索引
            std::vector<float> dis;     // 最近邻的距离
            // 在目标点云的KD树中查找一个最近邻，返回该最近邻的索引和距离
            kdtree_2d->nearestKSearch(pt, 1, nn_idx, dis);

            // 判断最近邻集合是否非空，且最小距离是否小于最大距离阈值
            if (nn_idx.size() > 0 && dis[0] < max_dis2) {
                effective_num++;    // 有效点数自增一
                Vec2d qw = Vec2d(target_cloud_2d->points[nn_idx[0]].x, target_cloud_2d->points[nn_idx[0]].y);   // 当前激光点在目标点云中的最近邻点坐标
                auto *edge = new EdgeSE2P2P(range, angle, qw, theta);   // 构建约束边，参数为：激光点的距离、角度、近邻点坐标、当前旋转角度
                edge->setVertex(0, v);                  // 设置边的第一个顶点为SE2位姿顶点
                edge->setInformation(Mat2d::Identity());// 观测为2维点坐标，因此信息矩阵需设为2x2单位矩阵
                auto rk = new g2o::RobustKernelHuber;   // Huber鲁棒核函数
                rk->setDelta(rk_delta);                 // 设置阈值
                edge->setRobustKernel(rk);              // 为边设置鲁棒核函数
                optimizer.addEdge(edge);                // 将约束边添加到优化器中
            } 
        }

        // 判断有效激光点数是否少于最小有效点数阈值
        if (effective_num < min_effect_pts) 
            return false;

        optimizer.setVerbose(false);        // 不输出优化过程
        optimizer.initializeOptimization(); // 初始化优化器
        optimizer.optimize(1);              // g2o内部仅非线性优化求解一次

        // 取出优化后的SE2位姿，更新当前位姿，用于下一次迭代
        current_pose = v->estimate();
    }
    init_pose = current_pose;
    LOG(INFO) << "estimated pose: " << current_pose.translation().transpose() << ", theta: " << current_pose.so2().log();
    // LOG(INFO) << "g2o: estimated pose: " << init_pose.translation().transpose() << ", theta: " << init_pose.so2().log();
    return true;
}

/**
 * @description: 【新增】基于G2O优化器的2D 点到点ICP算法
 * 
 *  for循环遍历scan点云，kdtree与target_cloud均传入位姿边构造函数中，g2o内部采用每次迭代估计的位姿v->estimate()转换得到pw点，执行近邻搜索得到qw，optimize(10)
 *  
 *  
 * 
 * @param {SE2&} init_pose 初始位姿
 * @return {*}
 */
bool Icp2d::AlignG2OP2P_2(SE2& init_pose) {
    int iterations = 10;                // 迭代次数
    double rk_delta = 0.8;
    float max_dis2 = 0.01;        // 最近邻时的最远距离（平方）
    int min_effect_pts = 20;      // 最小有效点数

    using BlockSolverType = g2o::BlockSolver<g2o::BlockSolverTraits<3, 1>>;
    using LinearSolverType = g2o::LinearSolverCholmod<BlockSolverType::PoseMatrixType>;
    auto* solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);
    
    auto* v = new VertexSE2();      // 新建SE2位姿顶点
    v->setId(0);                    // 设置顶点的id
    v->setEstimate(init_pose);   // 设置顶点的估计值为初始位姿
    optimizer.addVertex(v);         // 将顶点添加到优化器中
    int effective_num = 0;  // 有效点数
    // 遍历源始点云
    for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
        double range = source_scan_->ranges[i]; // 源始点云的距离
        // 判断每个点的距离是否越界
        if (range < source_scan_->range_min || range > source_scan_->range_max) 
            continue;

        // 根据最小角度和角分辨率计算每个点的角度
        double angle = source_scan_->angle_min + i * source_scan_->angle_increment;
        
        auto *edge = new EdgeSE2P2P_2(kdtree_2d, target_cloud_2d, range, angle);   

        edge->setVertex(0, v);                  // 设置边的第一个顶点为SE2位姿顶点
        if (edge->isPointValid()){
            effective_num++; 
            edge->setInformation(Mat2d::Identity());// 观测为2维点坐标，信息矩阵需设为2x2单位矩阵
            auto rk = new g2o::RobustKernelHuber;   // Huber鲁棒核函数
            rk->setDelta(rk_delta);                 // 设置阈值
            edge->setRobustKernel(rk);              // 为边设置鲁棒核函数    
            optimizer.addEdge(edge);                // 将约束边添加到优化器中
        }
    }

    // 判断有效激光点数是否少于最小有效点数阈值
    if (effective_num < min_effect_pts) 
        return false;

    optimizer.setVerbose(false);        // 不输出优化过程
    optimizer.initializeOptimization(); // 初始化优化器
    optimizer.optimize(10);              // g2o内部仅非线性优化求解一次

    // 取出优化后的SE2位姿，更新当前位姿，用于下一次迭代
    init_pose = v->estimate();
    LOG(INFO) << "estimated pose: " << v->estimate().translation().transpose() << ", theta: " << v->estimate().so2().log();
    // LOG(INFO) << "g2o: estimated pose: " << init_pose.translation().transpose() << ", theta: " << init_pose.so2().log();
    return true;
}

/**
 * @description: 【新增】基于G2O优化器的2D 点到点ICP算法
 *  
 *  for循环遍历scan点云，每个点都使用initial_pose转换的pw点，去进行近邻搜索得到qw，optimize(10)
 * 
 * @param {SE2&} init_pose 初始位姿
 * @return {*}
 */
bool Icp2d::AlignG2OP2P_3(SE2& init_pose) {
    int iterations = 10;                // 迭代次数
    double rk_delta = 0.8;
    float max_dis2 = 0.01;        // 最近邻时的最远距离（平方）
    int min_effect_pts = 20;      // 最小有效点数
    SE2 current_pose = init_pose;   // 当前位姿

    using BlockSolverType = g2o::BlockSolver<g2o::BlockSolverTraits<3, 1>>;
    using LinearSolverType = g2o::LinearSolverCholmod<BlockSolverType::PoseMatrixType>;
    auto* solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);
    
    auto* v = new VertexSE2();      // 新建SE2位姿顶点
    v->setId(0);                    // 设置顶点的id
    v->setEstimate(current_pose);   // 设置顶点的估计值为初始位姿
    optimizer.addVertex(v);         // 将顶点添加到优化器中
    int effective_num = 0;  // 有效点数
    // 遍历源始点云
    for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
        double range = source_scan_->ranges[i]; // 源始点云的距离
        // 判断每个点的距离是否越界
        if (range < source_scan_->range_min || range > source_scan_->range_max) 
            continue;

        // 根据最小角度和角分辨率计算每个点的角度
        double angle = source_scan_->angle_min + i * source_scan_->angle_increment;
        double theta = current_pose.so2().log(); // 当前位姿的角度

        // 世界系下点的坐标 p_i^W，极坐标转笛卡尔坐标公式
        Vec2d pw = current_pose * Vec2d(range * std::cos(angle), range * std::sin(angle));

        Point2d pt;
        pt.x = pw.x();
        pt.y = pw.y();

        // 最近邻
        std::vector<int> nn_idx;    // 最近邻的索引
        std::vector<float> dis;     // 最近邻的距离
        // 在目标点云的KD树中查找一个最近邻，返回该最近邻的索引和距离
        kdtree_2d->nearestKSearch(pt, 1, nn_idx, dis);

        // 判断最近邻集合是否非空，且最小距离是否小于最大距离阈值
        if (nn_idx.size() > 0 && dis[0] < max_dis2) {
            effective_num++;    // 有效点数自增一
            Vec2d qw = Vec2d(target_cloud_2d->points[nn_idx[0]].x, target_cloud_2d->points[nn_idx[0]].y);   // 当前激光点在目标点云中的最近邻点坐标
            auto *edge = new EdgeSE2P2P(range, angle, qw, theta);   // 构建约束边，参数为：激光点的距离、角度、近邻点坐标、当前旋转角度
            edge->setVertex(0, v);                  // 设置边的第一个顶点为SE2位姿顶点
            edge->setInformation(Mat2d::Identity());// 观测为2维点坐标，因此信息矩阵需设为2x2单位矩阵
            auto rk = new g2o::RobustKernelHuber;   // Huber鲁棒核函数
            rk->setDelta(rk_delta);                 // 设置阈值
            edge->setRobustKernel(rk);              // 为边设置鲁棒核函数
            optimizer.addEdge(edge);                // 将约束边添加到优化器中
        } 
    }

    // 判断有效激光点数是否少于最小有效点数阈值
    if (effective_num < min_effect_pts) 
        return false;

    optimizer.setVerbose(false);        // 不输出优化过程
    optimizer.initializeOptimization(); // 初始化优化器
    optimizer.optimize(10);              // g2o内部仅非线性优化求解一次

    // 取出优化后的SE2位姿，更新当前位姿，用于下一次迭代
    init_pose = v->estimate();
    LOG(INFO) << "estimated pose: " << v->estimate().translation().transpose() << ", theta: " << v->estimate().so2().log();
    // LOG(INFO) << "g2o: estimated pose: " << init_pose.translation().transpose() << ", theta: " << init_pose.so2().log();
    return true;
}

/**
 * @description: 点-线ICP算法
 * @param {SE2&} init_pose
 * @return {*}
 */
bool Icp2d::AlignGaussNewtonPoint2Plane(SE2& init_pose) {
    int iterations = 10;
    double cost = 0, lastCost = 0;
    float max_dis = 0.3;      // 最近邻时的最远距离
    int min_effect_pts = 20;  // 最小有效点数
    SE2 current_pose = init_pose;   // 初始化当前位姿

    // 遍历十次
    for (int iter = 0; iter < iterations; ++iter) {
        Mat3d H = Mat3d::Zero();
        Vec3d b = Vec3d::Zero();
        cost = 0;

        int effective_num = 0;  // 有效点数

        // 遍历当前scan中的所有2D激光点
        for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
            double r = source_scan_->ranges[i];  // 激光点的距离
            if (r < source_scan_->range_min || r > source_scan_->range_max) 
                continue;

            // 当前激光点的角度
            double angle = source_scan_->angle_min + i * source_scan_->angle_increment;
            // 从上一次迭代得到的位姿 T_wb 的2x2旋转矩阵中，利用对数映射获取对应的旋转角度
            double theta = current_pose.so2().log();
            // 机器人坐标系下的极坐标转换为笛卡尔坐标，并转为世界坐标系下的坐标 p_i^W，
            Vec2d pw = current_pose * Vec2d(r * std::cos(angle), r * std::sin(angle));
            Point2d pt;
            pt.x = pw.x();
            pt.y = pw.y();

            // 查找5个最近邻
            std::vector<int> nn_idx;
            std::vector<float> dis;
            kdtree_2d->nearestKSearch(pt, 5, nn_idx, dis);

            std::vector<Vec2d> effective_pts;  // 有效点
            // 遍历所有五个近邻点
            for (int j = 0; j < nn_idx.size(); ++j) {
                // 判断每个近邻点的距离是否处于最远阈值距离内
                if (dis[j] < max_dis) 
                    // 若是，该近邻点符合要求，存储到向量中
                    effective_pts.emplace_back(Vec2d(target_cloud_2d->points[nn_idx[j]].x, target_cloud_2d->points[nn_idx[j]].y));
            }
            // 判断有效近邻点是否少于三个
            if (effective_pts.size() < 3) 
                // 若少于3个，则跳过当前激光点
                continue;
            
            // 拟合直线，组装J、H和误差
            Vec3d line_coeffs;
            // 利用当前点附近的几个有效近邻点，基于SVD奇异值分解，拟合出ax+by+c=0 中的最小直线系数 a,b,c，对应公式（6.11）
            if (math::FitLine2D(effective_pts, line_coeffs)) {
                effective_num++; // 有效点数
                Vec3d J;   
                // 雅可比矩阵，对应公式（6.17），残差对SE2位姿的偏导数，
                // 由链式法则得到：de/dx = de/dp * dp/dx
                //  de/dp：残差相对于近邻点2D坐标的偏导
                //  dp/dx：2D点坐标相对于SE2位姿的偏导
                J << line_coeffs[0], line_coeffs[1],
                    -line_coeffs[0] * r * std::sin(angle + theta) + line_coeffs[1] * r * std::cos(angle + theta);

                // 高斯牛顿法的 海森矩阵 H = J * J^T，     H * dx = b = -J * f(x) = -J * e
                H += J * J.transpose();

                // 对应公式（6.14）  e = ax + by + c
                double e = line_coeffs[0] * pw[0] + line_coeffs[1] * pw[1] + line_coeffs[2];
                b += -J * e;

                cost += e * e;
            }
        }

        if (effective_num < min_effect_pts) 
            return false;

        // solve for dx
        Vec3d dx = H.ldlt().solve(b);
        if (isnan(dx[0])) 
            break;

        cost /= effective_num;
        if (iter > 0 && cost >= lastCost) 
            break;

        LOG(INFO) << "iter " << iter << " cost = " << cost << ", effect num: " << effective_num;

        current_pose.translation() += dx.head<2>();
        current_pose.so2() = current_pose.so2() * SO2::exp(dx[2]);
        lastCost = cost;
    }

    init_pose = current_pose;
    LOG(INFO) << "estimated pose: " << current_pose.translation().transpose()
              << ", theta: " << current_pose.so2().log();

    return true;
}

/**
 * @description: 基于G2O优化器的2D 点到线 ICP算法
 * @param {SE2&} init_pose 初始位姿
 * @return {*}
 */
bool Icp2d::AlignG2OP2L(SE2& init_pose) {
    int iterations = 10;        // 迭代次数
    double rk_delta = 0.8;
    float max_dis = 0.3;       // 最近邻时的最远距离（平方）
    int min_effect_pts = 20;    // 最小有效点数
    
    SE2 current_pose = init_pose;   // 当前位姿
    for (int iter = 0; iter < iterations; ++iter) {
        using BlockSolverType = g2o::BlockSolver<g2o::BlockSolverTraits<3, 1>>;
        using LinearSolverType = g2o::LinearSolverCholmod<BlockSolverType::PoseMatrixType>;
        auto* solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
        g2o::SparseOptimizer optimizer;
        optimizer.setAlgorithm(solver);

        auto* v = new VertexSE2();      // 新建SE2位姿顶点
        v->setId(0);                    // 设置顶点的id
        v->setEstimate(current_pose);   // 设置顶点的估计值为初始位姿
        optimizer.addVertex(v);         // 将顶点添加到优化器中
        int effective_num = 0;  // 有效点数
        // 遍历源始点云
        for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
            double range = source_scan_->ranges[i]; // 源始点云的距离
            // 判断每个点的距离是否越界
            if (range < source_scan_->range_min || range > source_scan_->range_max) 
                continue;

            // 当前激光点的角度
            double angle = source_scan_->angle_min + i * source_scan_->angle_increment;
            // 从上一次迭代得到的位姿 T_wb 的2x2旋转矩阵中，利用对数映射获取对应的旋转角度
            double theta = current_pose.so2().log();
            // 机器人坐标系下的极坐标转换为笛卡尔坐标，并转为世界坐标系下的坐标 p_i^W，
            Vec2d pw = current_pose * Vec2d(range * std::cos(angle), range * std::sin(angle));
            Point2d pt;
            pt.x = pw.x();
            pt.y = pw.y();

            // 查找5个最近邻
            std::vector<int> nn_idx;    // 最近邻的索引
            std::vector<float> dis;     // 最近邻的距离
            kdtree_2d->nearestKSearch(pt, 5, nn_idx, dis);

            std::vector<Vec2d> effective_pts;  // 有效点
            // 遍历所有五个近邻点
            for (int j = 0; j < nn_idx.size(); ++j) {
                // 判断每个近邻点的距离是否处于最远阈值距离内
                if (dis[j] < max_dis) 
                    // 若是，该近邻点符合要求，存储到向量中
                    effective_pts.emplace_back(Vec2d(target_cloud_2d->points[nn_idx[j]].x, target_cloud_2d->points[nn_idx[j]].y));
            }
            // 判断有效近邻点是否少于三个
            if (effective_pts.size() < 3) 
                // 若少于3个，则跳过当前激光点
                continue;

            // 拟合直线，组装J、H和误差
            Vec3d line_coeffs;
            // 利用当前点附近的几个有效近邻点，基于SVD奇异值分解，拟合出ax+by+c=0 中的最小直线系数 a,b,c，对应公式（6.11）
            if (math::FitLine2D(effective_pts, line_coeffs)) {
                effective_num++; // 有效点数
                auto *edge = new EdgeSE2P2L(range, angle, line_coeffs);
                edge->setVertex(0, v);                  // 设置边的第一个顶点为SE2位姿顶点
                edge->setInformation(Eigen::Matrix<double, 1, 1>::Identity());// 观测为2维点坐标，因此信息矩阵需设为2x2单位矩阵
                auto rk = new g2o::RobustKernelHuber;   // Huber鲁棒核函数
                rk->setDelta(rk_delta);                 // 设置阈值
                edge->setRobustKernel(rk);              // 为边设置鲁棒核函数
                optimizer.addEdge(edge);                // 将约束边添加到优化器中
            }
        }

        // 判断有效激光点数是否少于最小有效点数阈值
        if (effective_num < min_effect_pts) 
            return false;

        optimizer.setVerbose(false);        // 不输出优化过程
        optimizer.initializeOptimization(); // 初始化优化器
        optimizer.optimize(1);              // g2o内部仅非线性优化求解一次

        // 取出优化后的SE2位姿，更新当前位姿，用于下一次迭代
        current_pose = v->estimate();
    }
    init_pose = current_pose;
    LOG(INFO) << "estimated pose: " << current_pose.translation().transpose() << ", theta: " << current_pose.so2().log();
    // LOG(INFO) << "g2o: estimated pose: " << init_pose.translation().transpose() << ", theta: " << init_pose.so2().log();
    return true;
}

/**
 * @description: 基于G2O优化器的2D 点到线 ICP算法
 * @param {SE2&} init_pose 初始位姿
 * @return {*}
 */
bool Icp2d::AlignG2OP2L_2(SE2& init_pose) {
    int iterations = 10;        // 迭代次数
    double rk_delta = 0.8;
    float max_dis = 0.3;       // 最近邻时的最远距离（平方）
    int min_effect_pts = 20;    // 最小有效点数
    
    SE2 current_pose = init_pose;   // 当前位姿

    using BlockSolverType = g2o::BlockSolver<g2o::BlockSolverTraits<3, 1>>;
    using LinearSolverType = g2o::LinearSolverCholmod<BlockSolverType::PoseMatrixType>;
    auto* solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);

    auto* v = new VertexSE2();      // 新建SE2位姿顶点
    v->setId(0);                    // 设置顶点的id
    v->setEstimate(current_pose);   // 设置顶点的估计值为初始位姿
    optimizer.addVertex(v);         // 将顶点添加到优化器中
    int effective_num = 0;          // 有效点数

    // 遍历源始点云
    for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
        double range = source_scan_->ranges[i]; // 源始点云的距离
        // 判断每个点的距离是否越界
        if (range < source_scan_->range_min || range > source_scan_->range_max) 
            continue;

        // 当前激光点的角度
        double angle = source_scan_->angle_min + i * source_scan_->angle_increment;
        
        auto *edge = new EdgeSE2P2L_2(kdtree_2d, target_cloud_2d, range, angle);

        edge->setVertex(0, v);                  // 设置边的第一个顶点为SE2位姿顶点
        

        // 利用当前点附近的几个有效近邻点，基于SVD奇异值分解，拟合出ax+by+c=0 中的最小直线系数 a,b,c，对应公式（6.11）
        if (edge->isLineFitValid()) {
            effective_num++; // 有效点数
            edge->setInformation(Eigen::Matrix<double, 1, 1>::Identity());// 观测为2维点坐标，因此信息矩阵需设为2x2单位矩阵
            auto rk = new g2o::RobustKernelHuber;   // Huber鲁棒核函数
            rk->setDelta(rk_delta);                 // 设置阈值
            edge->setRobustKernel(rk);              // 为边设置鲁棒核函数
            optimizer.addEdge(edge);                // 将约束边添加到优化器中
        }
    }

    // 判断有效激光点数是否少于最小有效点数阈值
    if (effective_num < min_effect_pts) 
        return false;

    optimizer.setVerbose(false);        // 不输出优化过程
    optimizer.initializeOptimization(); // 初始化优化器
    optimizer.optimize(10);              // g2o内部仅非线性优化求解一次

    // 取出优化后的SE2位姿，更新当前位姿，用于下一次迭代
    current_pose = v->estimate();
    
    init_pose = current_pose;
    LOG(INFO) << "estimated pose: " << current_pose.translation().transpose() << ", theta: " << current_pose.so2().log();
    // LOG(INFO) << "g2o: estimated pose: " << init_pose.translation().transpose() << ", theta: " << init_pose.so2().log();
    return true;
}

/**
 * @description: 构建目标点云的KD树
 * @return {*}
 */
void Icp2d::BuildTargetKdTree() {
    // 判断目标2D点云是否为空
    if (target_scan_ == nullptr) {
        LOG(ERROR) << "target is not set";
        return;
    }

    // 重置PCL形式的目标点云
    target_cloud_2d.reset(new Cloud2d); // pcl::PointCloud<Point2d>
    target_cloud_3d.reset(new Cloud3d);
    // 遍历目标2D点云
    for (size_t i = 0; i < target_scan_->ranges.size(); ++i) {
        // 判断每个点距离是否越界
        if (target_scan_->ranges[i] < target_scan_->range_min || target_scan_->ranges[i] > target_scan_->range_max) 
            continue;

        // 根据最小角度和角分辨率计算每个点的角度
        double real_angle = target_scan_->angle_min + i * target_scan_->angle_increment;

        
        Point2d pt2d; // pcl::PointXY 
        Point3d pt3d; // pcl::PointXYZ
        // 将每个目标点的距离和角度转换为笛卡尔坐标
        pt2d.x = target_scan_->ranges[i] * std::cos(real_angle);
        pt2d.y = target_scan_->ranges[i] * std::sin(real_angle);
        pt3d.x = pt2d.x;
        pt3d.y = pt2d.y;
        pt3d.z = 0.0;
        target_cloud_2d->points.push_back(pt2d); // 将每个目标点添加到PCL形式的目标点云中
        target_cloud_3d->points.push_back(pt3d);
    }

    // 设置目标点云的点数和是否为稠密点云
    target_cloud_2d->width = target_cloud_2d->points.size();
    target_cloud_2d->height = 1;
    target_cloud_2d->is_dense = false;
    target_cloud_3d->width = target_cloud_3d->points.size();
    target_cloud_3d->height = 1;
    target_cloud_3d->is_dense = false;

    pcl::io::savePCDFile("target_cloud2d.pcd", *target_cloud_2d);
    // pcl::io::savePCDFile("target_cloud3d.pcd", *target_cloud_3d);

    kdtree_2d = boost::make_shared<pcl::search::KdTree<Point2d>>();
    kdtree_2d->setInputCloud(target_cloud_2d);
    
    
    // kdtree_3d = boost::make_shared<pcl::search::KdTree<Point3d>>(); // 创建KD树
    // kdtree_3d->setInputCloud(target_cloud_3d);   // 设置输入构建目标点云的KD树

    // clusterExtractor_3d.setClusterTolerance(0.5);  // 设置聚类的距离阈值
    // clusterExtractor_3d.setMinClusterSize(100);     // 设置聚类的最小点数
    // clusterExtractor_3d.setMaxClusterSize(25000);   // 设置聚类的最大点数
    // clusterExtractor_3d.setSearchMethod(kdtree_3d);
    // // 设置点云聚类对象的输入点云数据
    // clusterExtractor_3d.setInputCloud(target_cloud_3d);

    // // 创建一个向量来存储聚类的结果
    // std::vector<pcl::PointIndices> cluster_indices;
    // // 执行点云聚类
    // clusterExtractor_3d.extract(cluster_indices);

    // // 输出聚类结果
    // int clusterNumber = 1;
    // for (const auto& indices : cluster_indices) {
    //     LOG(INFO) << "Cluster " << clusterNumber << " has " << indices.indices.size() << " points.";
    //     clusterNumber++;

    //     // 提取聚类点云
    //     Cloud3d::Ptr clusterCloud(new Cloud3d);
    //     pcl::copyPointCloud(*target_cloud_3d, indices, *clusterCloud);

    //     // 将pcl::PointXYZ类型的点云转换为Eigen::Vector3d类型的点云
    //     std::vector<Vec2d> pts;
    //     pts.reserve(clusterCloud->size());

    //     for (const pcl::PointXYZ& pt : *clusterCloud)
    //         pts.emplace_back(Vec2d(pt.x, pt.y));

    //     // 拟合直线，组装J、H和误差
    //     Vec3d line_coeffs;
    //     // 利用当前点附近的几个有效近邻点，基于SVD奇异值分解，拟合出ax+by+c=0 中的最小直线系数 a,b,c，对应公式（6.11）
    //     if (math::FitLine2D(pts, line_coeffs)) 
    //         LOG(INFO) << "line_coeffs: " << line_coeffs[0] << ", " << line_coeffs[1] << ", " << line_coeffs[2];
        
    // }
}

}  // namespace sad