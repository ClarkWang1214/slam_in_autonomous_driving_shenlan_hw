//
// Created by xiang on 2021/8/18.
//

#include "ch5/bfnn.h"
#include <execution>

namespace sad {

/**
 * @description: 单个点的暴力最近邻搜索
 * @param {CloudPtr} cloud  点云
 * @param {Vec3f&} point    待搜索点
 * @return {*}
 */
int bfnn_point(CloudPtr cloud, const Vec3f& point) {
    return std::min_element(cloud->points.begin(), cloud->points.end(),     // 遍历点云中的每个点，计算待搜索点到每个点的距离，然后取出最小值
                            [&point](const PointType& pt1, const PointType& pt2) -> bool { // lambda函数，用于比较点pt1到目标point的距离是否比pt2到point的距离小
                                return (pt1.getVector3fMap() - point).squaredNorm() < (pt2.getVector3fMap() - point).squaredNorm();
                            }) - cloud->points.begin();
}

/**
 * @description: 多个点的暴力最近邻搜索
 * @param {CloudPtr} cloud  点云
 * @param {Vec3f&} point    待搜索点
 * @param {int} k           点数
 * @return {*}
 */
std::vector<int> bfnn_point_k(CloudPtr cloud, const Vec3f& point, int k) {
    // 用于存储点云中每个点的索引，以及到目标点的距离
    struct IndexAndDis {
        IndexAndDis() {}
        IndexAndDis(int index, double dis2) : index_(index), dis2_(dis2) {}
        int index_ = 0;
        double dis2_ = 0;
    };

    std::vector<IndexAndDis> index_and_dis(cloud->size());
    for (int i = 0; i < cloud->size(); ++i) {
        index_and_dis[i] = {i, (cloud->points[i].getVector3fMap() - point).squaredNorm()}; // 计算待搜索点到每个点的距离
    }
    // 按照距离排序
    std::sort(index_and_dis.begin(), index_and_dis.end(),
              [](const auto& d1, const auto& d2) { return d1.dis2_ < d2.dis2_; });
    std::vector<int> ret;
    std::transform( index_and_dis.begin(), index_and_dis.begin() + k, 
                    std::back_inserter(ret),                            //从后面插入
                    [](const auto& d1) { return d1.index_; });          // lambda函数，取出已排好序向量的前k个点的索引
    return ret;
}

/**
 * @description: 对两个点云执行暴力匹配（多线程）
 * @param {CloudPtr} cloud1
 * @param {CloudPtr} cloud2
 * @return {*}
 */
void bfnn_cloud_mt(CloudPtr cloud1, CloudPtr cloud2, std::vector<std::pair<size_t, size_t>>& matches) {
    // 先生成索引
    std::vector<size_t> index(cloud2->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });

    // 并行化for_each
    matches.resize(index.size());
    // 使用std::execution，对每个点并发地调用单点最近邻算法。
    std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](auto idx) {
        matches[idx].second = idx;
        matches[idx].first = bfnn_point(cloud1, ToVec3f(cloud2->points[idx]));
    });
}

/**
 * @description: 对两个点云执行暴力匹配（单线程）
 * @param {CloudPtr} cloud1
 * @param {CloudPtr} cloud2
 * @return {*}
 */
void bfnn_cloud(CloudPtr cloud1, CloudPtr cloud2, std::vector<std::pair<size_t, size_t>>& matches) {
    // 单线程版本
    std::vector<size_t> index(cloud2->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });

    matches.resize(index.size());
    std::for_each(std::execution::seq, index.begin(), index.end(), [&](auto idx) {
        matches[idx].second = idx;
        matches[idx].first = bfnn_point(cloud1, ToVec3f(cloud2->points[idx]));
    });
}

/**
 * @description: 对两个点云执行暴力匹配（多线程），k近邻
 * @param {CloudPtr} cloud1
 * @param {CloudPtr} cloud2
 * @param {int} k
 * @return {*}
 */
void bfnn_cloud_mt_k(CloudPtr cloud1, CloudPtr cloud2, std::vector<std::pair<size_t, size_t>>& matches, int k) {
    // 先生成索引
    std::vector<size_t> index(cloud2->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });

    // 并行化for_each
    matches.resize(index.size() * k);
    std::for_each(std::execution::par_unseq, index.begin(), index.end(), [&](auto idx) {
        auto v = bfnn_point_k(cloud1, ToVec3f(cloud2->points[idx]), k);
        for (int i = 0; i < v.size(); ++i) {
            matches[idx * k + i].first = v[i];
            matches[idx * k + i].second = idx;
        }
    });
}

}  // namespace sad
