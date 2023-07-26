//
// Created by xiang on 2021/8/19.
//
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>

#include <pcl/io/pcd_io.h>
#include <pcl/search/kdtree.h>

#include "ch5/bfnn.h"
#include "ch5/gridnn.hpp"
#include "ch5/kdtree.h"
#include "ch5/octo_tree.h"
#include "common/point_cloud_utils.h"
#include "common/point_types.h"
#include "common/sys_utils.h"
#include "common/math_utils.h"

#include "ch5/nanoflann.hpp"

template <typename T>
struct PointCloud_NanoFlann
{
    struct Point
    {
        T x, y, z;
    };

    using coord_t = T;  //!< The type of each coordinate

    std::vector<Point> pts;

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if (dim == 0)
            return pts[idx].x;
        else if (dim == 1)
            return pts[idx].y;
        else
            return pts[idx].z;
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
        return false;
    }
};

DEFINE_string(first_scan_path, "./data/ch5/first.pcd", "第一个点云路径");
DEFINE_string(second_scan_path, "./data/ch5/second.pcd", "第二个点云路径");
DEFINE_double(ANN_alpha, 1.0, "ANN的比例因子");

// 测试暴力近邻匹配算法的正确性与耗时
TEST(CH5_TEST, BFNN) {
    sad::CloudPtr first(new sad::PointCloudType), second(new sad::PointCloudType);
    pcl::io::loadPCDFile(FLAGS_first_scan_path, *first);    // 读取点云
    pcl::io::loadPCDFile(FLAGS_second_scan_path, *second);  

    if (first->empty() || second->empty()) {
        LOG(ERROR) << "cannot load cloud";
        FAIL();
    }

    // voxel grid 至 0.05
    sad::VoxelGrid(first);  // 对点云进行降采样，默认阈值为0.05
    sad::VoxelGrid(second); 

    LOG(INFO) << "points: " << first->size() << ", " << second->size();

    // 评价单线程和多线程版本的暴力匹配
    sad::evaluate_and_call( // 按照固定次数调用指定方法，并衡量其运行时间
        [&first, &second]() { // 将lambda函数作为参数传入evaluate_and_call函数，不需要在外面单独定义函数，很方便，值得学习
            std::vector<std::pair<size_t, size_t>> matches;
            sad::bfnn_cloud(first, second, matches);    // 单线程版本
        },                    // lambda函数 
        "暴力匹配（单线程）",  // 方法名
        5);                  // 调用5次，计算平均运行时间

    sad::evaluate_and_call(
        [&first, &second]() {
            std::vector<std::pair<size_t, size_t>> matches;
            sad::bfnn_cloud_mt(first, second, matches); // 多线程版本
        },
        "暴力匹配（多线程）", 5);   // 调用5次，计算平均运行时间

    SUCCEED();
}

/**
 * 评测最近邻的正确性
 * @param truth 真值
 * @param esti  估计
 */
void EvaluateMatches(const std::vector<std::pair<size_t, size_t>>& truth,
                     const std::vector<std::pair<size_t, size_t>>& esti) {
    int fp = 0;  // false-positive，esti存在但truth中不存在     误检率
    int fn = 0;  // false-negative, truth存在但esti不存在       漏检率

    LOG(INFO) << "truth: " << truth.size() << ", esti: " << esti.size();

    /// 检查某个匹配在另一个容器中存不存在
    auto exist = [](const std::pair<size_t, size_t>& data, const std::vector<std::pair<size_t, size_t>>& vec) -> bool {
        return std::find(vec.begin(), vec.end(), data) != vec.end();
    };

    int effective_esti = 0;
    for (const auto& d : esti) {
        if (d.first != sad::math::kINVALID_ID && d.second != sad::math::kINVALID_ID) {
            effective_esti++;

            if (!exist(d, truth)) {
                fp++;
            }
        }
    }

    for (const auto& d : truth) {
        if (!exist(d, esti)) {
            fn++;
        }
    }

    float precision = 1.0 - float(fp) / effective_esti;
    float recall = 1.0 - float(fn) / truth.size();
    LOG(INFO) << "precision: " << precision << ", recall: " << recall << ", fp: " << fp << ", fn: " << fn;
}

TEST(CH5_TEST, GRID_NN) {
    sad::CloudPtr first(new sad::PointCloudType), second(new sad::PointCloudType);
    pcl::io::loadPCDFile(FLAGS_first_scan_path, *first);
    pcl::io::loadPCDFile(FLAGS_second_scan_path, *second);

    if (first->empty() || second->empty()) {
        LOG(ERROR) << "cannot load cloud";
        FAIL();
    }

    LOG(INFO) << "points: " << first->size() << ", " << second->size();
    // voxel grid 至 0.05
    sad::VoxelGrid(first, 0.5);
    sad::VoxelGrid(second, 0.5);

    LOG(INFO) << "points: " << first->size() << ", " << second->size();

    std::vector<std::pair<size_t, size_t>> truth_matches;
    sad::bfnn_cloud(first, second, truth_matches);

    // 对比不同种类的grid
    sad::GridNN<2>  grid2_0(0.1, sad::GridNN<2>::NearbyType::CENTER), 
                    grid2_4(0.1, sad::GridNN<2>::NearbyType::NEARBY4),
                    grid2_8(0.1, sad::GridNN<2>::NearbyType::NEARBY8);
    sad::GridNN<3>  grid3_0(0.1, sad::GridNN<3>::NearbyType::CENTER),       // 【ClarkWang】
                    grid3_6(0.1, sad::GridNN<3>::NearbyType::NEARBY6),
                    grid3_14(0.1, sad::GridNN<3>::NearbyType::NEARBY14);    // 【ClarkWang】

    grid2_0.SetPointCloud(first);
    grid2_4.SetPointCloud(first);
    grid2_8.SetPointCloud(first);
    
    grid3_0.SetPointCloud(first);       // 【ClarkWang】
    grid3_6.SetPointCloud(first);
    grid3_14.SetPointCloud(first);      // 【ClarkWang】

    // 评价各种版本的Grid NN
    // sorry没有C17的template lambda... 下面必须写的啰嗦一些
    LOG(INFO) << "===================";
    std::vector<std::pair<size_t, size_t>> matches;
    sad::evaluate_and_call(
        [&first, &second, &grid2_0, &matches]() { grid2_0.GetClosestPointForCloud(first, second, matches); },
        "Grid0 2D 单线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid2_0, &matches]() { grid2_0.GetClosestPointForCloudMT(first, second, matches); },
        "Grid0 2D 多线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid2_4, &matches]() { grid2_4.GetClosestPointForCloud(first, second, matches); },
        "Grid4 2D 单线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid2_4, &matches]() { grid2_4.GetClosestPointForCloudMT(first, second, matches); },
        "Grid4 2D 多线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid2_8, &matches]() { grid2_8.GetClosestPointForCloud(first, second, matches); },
        "Grid8 2D 单线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid2_8, &matches]() { grid2_8.GetClosestPointForCloudMT(first, second, matches); },
        "Grid8 2D 多线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid3_0, &matches]() { grid3_0.GetClosestPointForCloud(first, second, matches); },
        "Grid0 3D 单线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid3_0, &matches]() { grid3_0.GetClosestPointForCloudMT(first, second, matches); },
        "Grid0 3D 多线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid3_6, &matches]() { grid3_6.GetClosestPointForCloud(first, second, matches); },
        "Grid6 3D 单线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid3_6, &matches]() { grid3_6.GetClosestPointForCloudMT(first, second, matches); },
        "Grid6 3D 多线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid3_14, &matches]() { grid3_14.GetClosestPointForCloud(first, second, matches); },
        "Grid14 3D 单线程", 10);
    EvaluateMatches(truth_matches, matches);

    LOG(INFO) << "===================";
    sad::evaluate_and_call(
        [&first, &second, &grid3_14, &matches]() { grid3_14.GetClosestPointForCloudMT(first, second, matches); },
        "Grid14 3D 多线程", 10);
    EvaluateMatches(truth_matches, matches);

    SUCCEED();
}

TEST(CH5_TEST, KDTREE_BASICS) {
    sad::CloudPtr cloud(new sad::PointCloudType);   // 构造点云
    sad::PointType p1, p2, p3, p4;                  // 点

    // 初始化四个三维点位置
    // (0, 0, 0)
    p1.x = 0;
    p1.y = 0;
    p1.z = 0;

    // (1, 0, 0)
    p2.x = 1;
    p2.y = 0;
    p2.z = 0;

    // (0, 1, 0)
    p3.x = 0;
    p3.y = 1;
    p3.z = 0;

    // (1, 1, 0)
    p4.x = 1;
    p4.y = 1;
    p4.z = 0;

    // 将四个点加入点云
    cloud->points.push_back(p1);
    cloud->points.push_back(p2);
    cloud->points.push_back(p3);
    cloud->points.push_back(p4);

    sad::KdTree kdtree;         // 构造kd树对象
    kdtree.BuildTree(cloud);    // 利用点云构造kd树
    kdtree.PrintAll();          // 打印kd树

    using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<float, PointCloud_NanoFlann<float>>,
        PointCloud_NanoFlann<float>, 3 /* dim */
        >;

    // nanoflann tree
    PointCloud_NanoFlann<float> cloud_flann;
    cloud_flann.pts.resize(4);
    for (int i = 0; i < cloud->size(); i++)
    {
        cloud_flann.pts[i].x = cloud->points[i].x;
        cloud_flann.pts[i].y = cloud->points[i].y;
        cloud_flann.pts[i].z = cloud->points[i].z;
    }

    my_kd_tree_t kdtree_nanoflann(3, cloud_flann, {10});
    kdtree_nanoflann.buildIndex();

    // 待查询点
    const float query_pt[3] = {0.8, 0, 0};

    // ----------------------------------------------------------------
    // knnSearch():  Perform a search for the N closest points 
    //               执行搜索以查找最近的N个点
    // ----------------------------------------------------------------
    {
        size_t                num_results = 5;              // 最近的5个点
        std::vector<uint32_t> ret_index(num_results);       // 返回的点的索引
        std::vector<float>    out_dist_sqr(num_results);    // 返回的点的距离

        // 调用knnSearch()函数，返回最近的5个点的索引和距离
        num_results = kdtree_nanoflann.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

        // In case of less points in the tree than requested:
        // 如果树中的点少于请求的点，则重新调整返回的点的索引和距离容器的大小
        ret_index.resize(num_results);
        out_dist_sqr.resize(num_results);

        LOG(INFO) << "knnSearch(): num_results=" << num_results;
        for (size_t i = 0; i < num_results; i++)
            LOG(INFO) << "idx[" << i << "]=" << ret_index[i] << " dist[" << i << "]=" << out_dist_sqr[i];
    }

    // // ----------------------------------------------------------------
    // // radiusSearch(): Perform a search for the points within search_radius
    // //                  执行搜索以查找search_radius半径内的点
    // // ----------------------------------------------------------------
    // {
    //     const float search_radius = static_cast<float>(0.1);
    //     std::vector<nanoflann::ResultItem<uint32_t, float>> ret_matches;

    //     // nanoflanSearchParamsameters params;
    //     // params.sorted = false;

    //     // 调用radiusSearch()函数，返回指定搜索半径内的点的索引和距离
    //     const size_t nMatches = kdtree_nanoflann.radiusSearch(&query_pt[0], search_radius, ret_matches);

    //     LOG(INFO) << "radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches";
    //     for (size_t i = 0; i < nMatches; i++)
    //         LOG(INFO) << "idx[" << i << "]=" << ret_matches[i].first << " dist[" << i << "]=" << ret_matches[i].second;
    // }

    SUCCEED();
}

TEST(CH5_TEST, KDTREE_KNN) {
    // 创建两个点云
    sad::CloudPtr first(new sad::PointCloudType), second(new sad::PointCloudType);
    pcl::io::loadPCDFile(FLAGS_first_scan_path, *first);
    pcl::io::loadPCDFile(FLAGS_second_scan_path, *second);

    if (first->empty() || second->empty()) {
        LOG(ERROR) << "cannot load cloud";
        FAIL();
    }

    // voxel grid 至 0.05
    sad::VoxelGrid(first);
    sad::VoxelGrid(second);

    // 比较 bfnn
    std::vector<std::pair<size_t, size_t>> true_matches;
    sad::bfnn_cloud_mt_k(first, second, true_matches);


    /////////////////////////// kdtree sad测试 ////////////////////////////////

    sad::KdTree kdtree;
    sad::evaluate_and_call([&first, &kdtree]() { kdtree.BuildTree(first); }, "Kd Tree build", 1);

    kdtree.SetEnableANN(true, FLAGS_ANN_alpha);

    LOG(INFO) << "Kd tree leaves: " << kdtree.size() << ", points: " << first->size();

    // 对第2个点云执行knn
    std::vector<std::pair<size_t, size_t>> matches;
    sad::evaluate_and_call([&first, &second, &kdtree, &matches]() { kdtree.GetClosestPointMT(second, matches, 5); },
                           "Kd Tree 5NN 多线程", 1);
    EvaluateMatches(true_matches, matches);

    /////////////////////////// kdtree pcl测试 ////////////////////////////////

    LOG(INFO) << "building kdtree pcl";
    // 对比PCL
    pcl::search::KdTree<sad::PointType> kdtree_pcl;
    sad::evaluate_and_call([&first, &kdtree_pcl]() { kdtree_pcl.setInputCloud(first); }, "Kd Tree build", 1);

    LOG(INFO) << "searching pcl";
    matches.clear();
    std::vector<int> search_indices(second->points.size());
    for (int i = 0; i < second->points.size(); i++) {
        search_indices[i] = i;
    }

    std::vector<std::vector<int>> result_index;
    std::vector<std::vector<float>> result_distance;
    sad::evaluate_and_call( [&]() { kdtree_pcl.nearestKSearch(*second, search_indices, 5, result_index, result_distance); 
                                    // 遍历每个点，获取最近的5个点的索引和距离
                                    for (int i = 0; i < second->points.size(); i++) {
                                        // 遍历每个点的最近的5个点
                                        for (int j = 0; j < result_index[i].size(); ++j) {
                                            int m = result_index[i][j];         // 最近的5个点的索引
                                            double d = result_distance[i][j];   // 最近的5个点的距离
                                            matches.push_back({m, i});          // 将最近的5个点的索引和距离存入matches
                                        }
                                    }    
                                },
                            "Kd Tree 5NN in PCL", 1);
    
    EvaluateMatches(true_matches, matches);

    /////////////////////////// kdtree nanoflann测试 ////////////////////////////////

    LOG(INFO) << "building kdtree nanflann 1";
    
    using kdtree_nano = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud_NanoFlann<float>>, 
                                                            PointCloud_NanoFlann<float>, 3>;

    kdtree_nano* my_kdtree_nano;
    // nanoflann tree
    PointCloud_NanoFlann<float> first_cloud_flann;
    first_cloud_flann.pts.resize(first->points.size());
    for (int i = 0; i < first->points.size(); i++)
    {
        first_cloud_flann.pts[i].x = first->points[i].x;
        first_cloud_flann.pts[i].y = first->points[i].y;
        first_cloud_flann.pts[i].z = first->points[i].z;
    }
    sad::evaluate_and_call([&my_kdtree_nano, &first_cloud_flann]() { 
                                my_kdtree_nano = new kdtree_nano(3, first_cloud_flann, nanoflann::KDTreeSingleIndexAdaptorParams(10));
                                my_kdtree_nano->buildIndex();
                            }, "Kd Tree build 1", 1);

    LOG(INFO) << "searching nanflann 1";
    matches.clear();
    int k = 5; 
    std::vector<std::vector<uint32_t>> ret_index_all;
    std::vector<std::vector<float>> out_dist_sqr_all;
    ret_index_all.resize(second->size());
    out_dist_sqr_all.resize(second->size());
    sad::evaluate_and_call([&second, &my_kdtree_nano, &matches, &k, &ret_index_all, &out_dist_sqr_all]() {
            // 索引
            std::vector<int> index(second->size());
            for (int i = 0; i < second->points.size(); ++i) {
                index[i] = i;
            }
            
            std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&second, &my_kdtree_nano, &k, &ret_index_all, &out_dist_sqr_all](int idx) {
                std::vector<uint32_t> ret_index(k);          // 返回的点的索引
                std::vector<float> out_dist_sqr(k);         // 返回的点的距离
                // 获取第二个点云中的每个点，作为当前待查询点
                float query_p[3] = { second->points[idx].x, second->points[idx].y, second->points[idx].z};
                
                // 调用knnSearch()函数，返回最近的5个点的索引和距离
                int num_results = my_kdtree_nano->knnSearch(&query_p[0], k, &ret_index[0], &out_dist_sqr[0]);
                
                ret_index.resize(num_results);
                out_dist_sqr.resize(num_results);

                ret_index_all[idx] = ret_index;
                out_dist_sqr_all[idx] = out_dist_sqr;
                
            });
            // 遍历每个点，获取最近的5个点的索引和距离
            for (int i = 0; i < second->points.size(); i++) {
                // 遍历每个点的最近的5个点
                for (int j = 0; j < ret_index_all[i].size(); ++j) {
                    int m = ret_index_all[i][j];         // 最近的5个点的索引
                    double d = out_dist_sqr_all[i][j];   // 最近的5个点的距离
                    matches.push_back({m, i});          // 将最近的5个点的索引和距离存入matches
                }
            }

        },
        "Kd Tree 5NN in nanflann 1", 1);

    EvaluateMatches(true_matches, matches);

    /////////////////////////// kdtree nanoflann测试 2 ////////////////////////////////

    LOG(INFO) << "building kdtree nanflann 2";

    kdtree_nano* my_kdtree_nano2;
    sad::evaluate_and_call([&my_kdtree_nano2, &first_cloud_flann]() { 
                                my_kdtree_nano2 = new kdtree_nano(3, first_cloud_flann, nanoflann::KDTreeSingleIndexAdaptorParams(10));
                                my_kdtree_nano2->buildIndex();
                            }, "Kd Tree build 2", 1);

    LOG(INFO) << "searching nanflann 2";
    matches.clear();
    // int k = 5; 
    std::vector<std::vector<size_t>> ret_index_all2;
    std::vector<std::vector<float>> out_dist_sqr_all2;
    ret_index_all2.resize(second->size());
    out_dist_sqr_all2.resize(second->size());
    
    sad::evaluate_and_call([&second, &my_kdtree_nano2, &matches, &k, &ret_index_all2, &out_dist_sqr_all2]() {
            // 索引
            std::vector<int> index(second->size());
            for (int i = 0; i < second->points.size(); ++i) {
                index[i] = i;
            }
            
            std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&second, &my_kdtree_nano2, &k, &ret_index_all2, &out_dist_sqr_all2](int idx) {
                
                size_t ret_index[k];          // 返回的点的索引
                float out_dist_sqr[k];         // 返回的点的距离
                // 获取第二个点云中的每个点，作为当前待查询点
                float query_p[3] = { second->points[idx].x, second->points[idx].y, second->points[idx].z};
                nanoflann::KNNResultSet<float> resultSet(k);
                resultSet.init(ret_index, out_dist_sqr);
                my_kdtree_nano2->findNeighbors(resultSet, query_p);

                std::vector<size_t> ret_index_(ret_index, ret_index+k);          // 返回的点的索引
                std::vector<float> out_dist_sqr_(out_dist_sqr, out_dist_sqr+k);  // 返回的点的距离

                ret_index_all2[idx] = ret_index_;
                out_dist_sqr_all2[idx] = out_dist_sqr_;
                
            });
            // 遍历每个点，获取最近的5个点的索引和距离
            for (int i = 0; i < second->points.size(); i++) {
                // 遍历每个点的最近的5个点
                for (int j = 0; j < ret_index_all2[i].size(); ++j) {
                    int m = ret_index_all2[i][j];         // 最近的5个点的索引
                    double d = out_dist_sqr_all2[i][j];   // 最近的5个点的距离
                    matches.push_back({m, i});          // 将最近的5个点的索引和距离存入matches
                }
            }

        },
        "Kd Tree 5NN in nanflann 2", 1);

    EvaluateMatches(true_matches, matches);

    LOG(INFO) << "done.";

    SUCCEED();
}

TEST(CH5_TEST, OCTREE_BASICS) {
    sad::CloudPtr cloud(new sad::PointCloudType);
    sad::PointType p1, p2, p3, p4;
    p1.x = 0;
    p1.y = 0;
    p1.z = 0;

    p2.x = 1;
    p2.y = 0;
    p2.z = 0;

    p3.x = 0;
    p3.y = 1;
    p3.z = 0;

    p4.x = 1;
    p4.y = 1;
    p4.z = 0;

    cloud->points.push_back(p1);
    cloud->points.push_back(p2);
    cloud->points.push_back(p3);
    cloud->points.push_back(p4);

    sad::OctoTree octree;
    octree.BuildTree(cloud);
    octree.SetApproximate(false);
    LOG(INFO) << "Octo tree leaves: " << octree.size() << ", points: " << cloud->size();

    SUCCEED();
}

TEST(CH5_TEST, OCTREE_KNN) {
    sad::CloudPtr first(new sad::PointCloudType), second(new sad::PointCloudType);
    pcl::io::loadPCDFile(FLAGS_first_scan_path, *first);
    pcl::io::loadPCDFile(FLAGS_second_scan_path, *second);

    if (first->empty() || second->empty()) {
        LOG(ERROR) << "cannot load cloud";
        FAIL();
    }

    // voxel grid 至 0.05
    sad::VoxelGrid(first);
    sad::VoxelGrid(second);

    sad::OctoTree octree;
    sad::evaluate_and_call([&first, &octree]() { octree.BuildTree(first); }, "Octo Tree build", 1);

    octree.SetApproximate(true, FLAGS_ANN_alpha);
    LOG(INFO) << "Octo tree leaves: " << octree.size() << ", points: " << first->size();

    /// 测试KNN
    LOG(INFO) << "testing knn";
    std::vector<std::pair<size_t, size_t>> matches;
    sad::evaluate_and_call([&first, &second, &octree, &matches]() { octree.GetClosestPointMT(second, matches, 5); },
                           "Octo Tree 5NN 多线程", 1);

    LOG(INFO) << "comparing with bfnn";
    /// 比较真值
    std::vector<std::pair<size_t, size_t>> true_matches;
    sad::bfnn_cloud_mt_k(first, second, true_matches);
    EvaluateMatches(true_matches, matches);

    LOG(INFO) << "done.";

    SUCCEED();
}

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;

    testing::InitGoogleTest(&argc, argv);
    google::ParseCommandLineFlags(&argc, &argv, true);
    return RUN_ALL_TESTS(); // 执行所有的 test case
}
