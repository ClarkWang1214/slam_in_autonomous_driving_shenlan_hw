//
// Created by xiang on 2021/8/25.
//

#ifndef SLAM_IN_AUTO_DRIVING_GRID2D_HPP
#define SLAM_IN_AUTO_DRIVING_GRID2D_HPP

#include "common/eigen_types.h"
#include "common/math_utils.h"
#include "common/point_types.h"

#include <glog/logging.h>
#include <execution>
#include <map>

namespace sad {

/**
 * 栅格法最近邻
 * @tparam dim 模板参数，使用2D或3D栅格
 */
template <int dim> // 2 or 3
class GridNN {
public:
    using KeyType = Eigen::Matrix<int, dim, 1>; // 栅格坐标
    using PtType = Eigen::Matrix<float, dim, 1>;// 点坐标

    enum class NearbyType {
        CENTER,  // 只考虑中心
        // for 2D
        NEARBY4,  // 上下左右
        NEARBY8,  // 上下左右+四角

        // for 3D
        NEARBY6,  // 上下左右前后

        //【新增】上下左右前后，包含8个角
        NEARBY14, 
    };

    /**
     * 构造函数
     * @param resolution 分辨率
     * @param nearby_type 近邻判定方法
     */
    explicit GridNN(float resolution = 0.1, NearbyType nearby_type = NearbyType::NEARBY4)
        : resolution_(resolution), nearby_type_(nearby_type) {

        // 分辨率倒数
        inv_resolution_ = 1.0 / resolution_;    

        // check dim and nearby
        if (dim == 2 && nearby_type_ == NearbyType::NEARBY6) {
            LOG(INFO) << "2D grid does not support nearby6, using nearby4 instead.";
            nearby_type_ = NearbyType::NEARBY4;
        // 【新增】Clark
        } else if (dim == 3 && (nearby_type_ != NearbyType::NEARBY6 && 
                                nearby_type_ != NearbyType::NEARBY14 && 
                                nearby_type_ != NearbyType::CENTER)) {
            LOG(INFO) << "3D grid does not support nearby4/8, using nearby6/14 instead.";
            nearby_type_ = NearbyType::NEARBY6;
        }

        // 生成2D栅格
        GenerateNearbyGrids();
    }

    /// 设置点云，建立栅格
    bool SetPointCloud(CloudPtr cloud);

    /// 获取最近邻
    bool GetClosestPoint(const PointType& pt, PointType& closest_pt, size_t& idx);

    /// 对比两个点云
    bool GetClosestPointForCloud(CloudPtr ref, CloudPtr query, std::vector<std::pair<size_t, size_t>>& matches);
    bool GetClosestPointForCloudMT(CloudPtr ref, CloudPtr query, std::vector<std::pair<size_t, size_t>>& matches);

private:
    /// 根据最近邻的类型，生成附近网格
    void GenerateNearbyGrids();

    /// 空间坐标转到grid
    KeyType Pos2Grid(const PtType& pt);

    float resolution_ = 0.1;       // 分辨率
    float inv_resolution_ = 10.0;  // 分辨率倒数

    NearbyType nearby_type_ = NearbyType::NEARBY4;

    // 哈希表，存放实际的栅格数据，由于在没有数据的地方就不必保留空的栅格，哈希表在这种场景下非常好用
    std::unordered_map<KeyType, std::vector<size_t>, hash_vec<dim>> grids_; 
    CloudPtr cloud_;

    std::vector<KeyType> nearby_grids_;  // 附近的栅格，最近邻
};

/**
 * @description: 设置点云，建立栅格
 * @return {*}
 */
template <int dim>
bool GridNN<dim>::SetPointCloud(CloudPtr cloud) {
    // vector存放点云索引：0, 1, 2, 3, ...，n-1
    std::vector<size_t> index(cloud->size()); 
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });

    // 遍历所有点云，将点云索引插入到栅格中
    std::for_each(index.begin(), index.end(), [&cloud, this](const size_t& idx) {
        auto pt = cloud->points[idx]; // 点云中第idx个点
        auto key = Pos2Grid(ToEigen<float, dim>(pt));   // 计算点所在栅格坐标
        // 判断是否已有该栅格，若没有则新增，若有则将点云索引添加到该栅格中
        if (grids_.find(key) == grids_.end()) 
            grids_.insert({key, {idx}});
        else 
            grids_[key].emplace_back(idx);
    });

    cloud_ = cloud;
    LOG(INFO) << "grids: " << grids_.size();
    return true;
}

/**
 * @description: 点坐标转换为栅格坐标
 * @return {*}
 */
template <int dim>
Eigen::Matrix<int, dim, 1> GridNN<dim>::Pos2Grid(const Eigen::Matrix<float, dim, 1>& pt) {
    // 点所在栅格坐标 = 点坐标 * 分辨率倒数
    return (pt * inv_resolution_).template cast<int>();
}

/**
 * @description: 特化模板类，对2维栅格定义近邻生成方式。3维栅格允许取0、4、8个最近邻
 * @return {*}
 */
template <>
void GridNN<2>::GenerateNearbyGrids() {
    if (nearby_type_ == NearbyType::CENTER) 
        nearby_grids_.emplace_back(KeyType::Zero());
    else if (nearby_type_ == NearbyType::NEARBY4) 
        nearby_grids_ = {Vec2i(0, 0), Vec2i(-1, 0), Vec2i(1, 0), Vec2i(0, 1), Vec2i(0, -1)};
    else if (nearby_type_ == NearbyType::NEARBY8) 
        nearby_grids_ = {Vec2i(0, 0),   
                         Vec2i(-1, 0), Vec2i(1, 0),  Vec2i(0, 1), Vec2i(0, -1),
                         Vec2i(-1, -1), Vec2i(-1, 1), Vec2i(1, -1), Vec2i(1, 1)};
}

/**
 * @description: 特化模板类，对3维栅格定义近邻生成方式。3维栅格允许取0、6、14个最近邻
 * @return {*}
 */
template <>
void GridNN<3>::GenerateNearbyGrids() {
    if (nearby_type_ == NearbyType::CENTER) 
        nearby_grids_.emplace_back(KeyType::Zero());
    else if (nearby_type_ == NearbyType::NEARBY6) 
        nearby_grids_ = {KeyType( 0,  0,  0),  
                         KeyType(-1,  0,  0),   KeyType(1,  0,  0),   // 左右
                         KeyType( 0,  1,  0),   KeyType(0, -1,  0),   // 前后
                         KeyType( 0,  0, -1),   KeyType(0,  0,  1)};  // 上下
    else if (nearby_type_ == NearbyType::NEARBY14) 
        // 【新增】
        nearby_grids_ = {KeyType( 0,  0,  0),  
                         KeyType(-1,  0,  0),   KeyType( 1,  0,  0),   // 左右
                         KeyType( 0,  1,  0),   KeyType( 0, -1,  0),   // 前后
                         KeyType( 0,  0, -1),   KeyType( 0,  0,  1),   // 上下
                         KeyType(-1, -1, -1),   KeyType(-1,  1, -1),   
                         KeyType(-1, -1,  1),   KeyType(-1,  1,  1), 
                         KeyType( 1, -1, -1),   KeyType( 1,  1, -1), 
                         KeyType( 1, -1,  1),   KeyType( 1,  1,  1)};
}

/**
 * @description: 在pt栅格周边寻找最近邻
 * @param {PointType&} pt           待查找点
 * @param {PointType&} closest_pt   返回的最近邻点
 * @param {size_t&} idx
 * @return {*}
 */
template <int dim>
bool GridNN<dim>::GetClosestPoint(const PointType& pt, PointType& closest_pt, size_t& idx) {
    std::vector<size_t> idx_to_check;   // 存储最近邻栅格的索引

    // 【1】计算给定点所在的栅格
    auto key = Pos2Grid(ToEigen<float, dim>(pt));

    // 【2】根据最近邻的定义，查找附近的栅格
    std::for_each(  nearby_grids_.begin(), // 遍历附近栅格，2D栅格为4/8，3D栅格为6/14 
                    nearby_grids_.end(), 
                    [&key, &idx_to_check, this](const KeyType& delta) { // lambda函数，对每个附近栅格进行操作
                        auto dkey = key + delta;
                        auto iter = grids_.find(dkey);
                        if (iter != grids_.end()) {
                            idx_to_check.insert(idx_to_check.end(), iter->second.begin(), iter->second.end());
                        }
                    }
                );

    if (idx_to_check.empty()) return false;

    // brute force nn in cloud_[idx]
    CloudPtr nearby_cloud(new PointCloudType);
    std::vector<size_t> nearby_idx;
    // 遍历附近栅格中的点
    for (auto& idx : idx_to_check) {
        // 加入到nearby_cloud附近点云容器中
        nearby_cloud->points.template emplace_back(cloud_->points[idx]);
        nearby_idx.emplace_back(idx); // 保存附近点云的索引
    }

    // 【3】采样暴力匹配算法计算这些周围栅格中的最近邻，相比较之前对所有点云进行暴力匹配，这里的点云数量更少，所以效率更高
    size_t closest_point_idx = bfnn_point(nearby_cloud, ToVec3f(pt));   // 返回最近邻点的索引
    idx = nearby_idx.at(closest_point_idx); // 获取最近邻点的索引对应的原始点云索引
    closest_pt = cloud_->points[idx];    // 保存最近邻点的坐标

    return true;
}

template <int dim>
bool GridNN<dim>::GetClosestPointForCloud(CloudPtr ref, CloudPtr query,
                                          std::vector<std::pair<size_t, size_t>>& matches) {
    matches.clear();
    std::vector<size_t> index(query->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });
    std::for_each(index.begin(), index.end(), [this, &matches, &query](const size_t& idx) {
        PointType cp;
        size_t cp_idx;
        if (GetClosestPoint(query->points[idx], cp, cp_idx)) {
            matches.emplace_back(cp_idx, idx);
        }
    });

    return true;
}

template <int dim>
bool GridNN<dim>::GetClosestPointForCloudMT(CloudPtr ref, CloudPtr query,
                                            std::vector<std::pair<size_t, size_t>>& matches) {
    matches.clear();
    // 与串行版本基本一样，但matches需要预先生成，匹配失败时填入非法匹配
    std::vector<size_t> index(query->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });
    matches.resize(index.size());

    std::for_each(std::execution::par_unseq, index.begin(), index.end(), [this, &matches, &query](const size_t& idx) {
        PointType cp;
        size_t cp_idx;
        if (GetClosestPoint(query->points[idx], cp, cp_idx)) 
            matches[idx] = {cp_idx, idx};
        else 
            matches[idx] = {math::kINVALID_ID, math::kINVALID_ID};
        
    });

    return true;
}

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_GRID2D_HPP
