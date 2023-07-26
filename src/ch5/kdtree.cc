//
// Created by xiang on 2021/9/22.
//

#include "ch5/kdtree.h"
#include "common/math_utils.h"

#include <glog/logging.h>
#include <execution>
#include <set>

namespace sad {

/**
 * @description: 构建Kd树
 * @param {CloudPtr} &cloud
 * @return {*}
 */
bool KdTree::BuildTree(const CloudPtr &cloud) {
    // 1. 判断点云是否为空
    if (cloud->empty()) return false;
    
    cloud_.clear(); // 清空点云数据
    cloud_.resize(cloud->size());
    // 2. 将PCL的点云类型pcl::PointCloud<pcl::PointXYZI>转换为Eigen的三维向量类型Eigen::Vector3f
    for (size_t i = 0; i < cloud->points.size(); ++i) 
        cloud_[i] = ToVec3f(cloud->points[i]);

    Clear();
    Reset();

    // 3. 生成索引列表
    IndexVec idx(cloud->size());
    for (int i = 0; i < cloud->points.size(); ++i) 
        idx[i] = i;

    // 4. 递归调用，将所有点插入到Kd树中
    Insert(idx, root_.get());
    return true;
}

/**
 * @description: 所有点插入到kd树的递归实现
 * @param {IndexVec} &points
 * @param {KdTreeNode} *node
 * @return {*}
 */
void KdTree::Insert(const IndexVec &points, KdTreeNode *node) {
    // 将当前节点插入到哈希表中
    nodes_.insert({node->id_, node});

    if (points.empty()) return;

    if (points.size() == 1) {
        size_++;
        node->point_idx_ = points[0];
        return;
    }

    IndexVec left, right;
    // 计算Xn在各轴的方差，然后挑选分布最大的那个轴作为分割轴；取平均值作为分割阈值
    if (!FindSplitAxisAndThresh(points, node->axis_index_, node->split_thresh_, left, right)) {
        size_++;                        // Kd树的size加1
        node->point_idx_ = points[0];   // 将节点的索引设为points[0]
        return;
    }

    // 定义一个lambda表达式，用于创建新的节点
    const auto create_if_not_empty  
                = [&node, this](KdTreeNode *&new_node, const IndexVec &index) {
                    if (!index.empty()) {
                        new_node = new KdTreeNode;          // 创建新的节点
                        new_node->id_ = tree_node_id_++;    // 初始化新节点的id，索引加1
                        Insert(index, new_node);            // 递归调用插入新节点
                    }
                };

    create_if_not_empty(node->left_, left);
    create_if_not_empty(node->right_, right);
}

/**
 * @description: 
 * @param {PointType} &pt               待检索点
 * @param {vector<int>} &closest_idx    返回的最近邻点集合
 * @param {int} k                       前k个最近邻
 * @return {*}
 */
bool KdTree::GetClosestPoint(const PointType &pt, std::vector<int> &closest_idx, int k) {
    // 检查k是否大于点云的size
    if (k > size_) {
        LOG(ERROR) << "cannot set k larger than cloud size: " << k << ", " << size_;
        return false;
    }
    k_ = k; // 将k赋值给成员变量k_

    std::priority_queue<NodeAndDistance> knn_result;    // 用于存储最近邻点的优先队列
    Knn(ToVec3f(pt), root_.get(), knn_result);          // 调用Knn函数，寻找待检索点的最近邻

    // 排序并返回结果
    closest_idx.resize(knn_result.size());
    for (int i = closest_idx.size() - 1; i >= 0; --i) {
        // 倒序插入
        closest_idx[i] = knn_result.top().node_->point_idx_;
        knn_result.pop();
    }
    return true;
}
/**
 * @description: 并行，为点云寻找最近邻
 * @param {CloudPtr} &cloud                             待检索点云
 * @param {vector<std::pair<size_t, size_t>>} &matches  返回的匹配结果
 * @param {int} k                                       前k个最近邻
 * @return {*}
 */
bool KdTree::GetClosestPointMT(const CloudPtr &cloud, std::vector<std::pair<size_t, size_t>> &matches, int k) {
    // 调整返回的匹配结果的大小
    matches.resize(cloud->size() * k);

    // 构建索引
    std::vector<int> index(cloud->size());
    for (int i = 0; i < cloud->points.size(); ++i) 
        index[i] = i; // 待检索点云中每个点的索引

    // lambda函数中的参数idx是索引列表中的每个索引，并行遍历该索引列表
    std::for_each(std::execution::par_unseq, index.begin(), index.end(), [this, &cloud, &matches, &k](int idx) {
        // 用于存储最近邻点的索引
        std::vector<int> closest_idx;   
        // 获取待检索点云中第idx个点，在kd树种前k个最近邻点的索引集合
        GetClosestPoint(cloud->points[idx], closest_idx, k); 
        // 将最近邻点的索引和待检索点的索引存入matches中
        for (int i = 0; i < k; ++i) {
            matches[idx * k + i].second = idx;                  // 待检索点的索引作为匹配向量的second值
            if (i < closest_idx.size()) 
                matches[idx * k + i].first = closest_idx[i];    // 最近邻点的索引作为匹配向量的first值
            else 
                matches[idx * k + i].first = math::kINVALID_ID; // 如果最近邻点的索引不存在，则将其设为非法值
        }
    });

    return true;
}

/**
 * @description: Kd树的最近邻搜索，递归实现
 * @param {Vec3f} &pt
 * @param {KdTreeNode} *node
 * @param {priority_queue<NodeAndDistance>} &knn_result
 * @return {*}
 */
void KdTree::Knn(const Vec3f &pt, KdTreeNode *node, std::priority_queue<NodeAndDistance> &knn_result) const {
    // 判断当前节点是否为叶子节点
    if (node->IsLeaf()) {
        // 如果是叶子，检查该节点与待检索点的距离，看是否小于
        ComputeDisForLeaf(pt, node, knn_result);
        return;
    }

    // 看pt落在左还是右，优先搜索pt所在的子树
    // 然后再看另一侧子树是否需要搜索
    
    // 如果当前节点不是叶子节点，
    KdTreeNode *this_side, *that_side;
    // 根据待检索点的坐标与当前节点的分割阈值进行比较，判断待检索点落在哪一侧
    if (pt[node->axis_index_] < node->split_thresh_) {
        // 若小于阈值，则待检索点落在左侧
        this_side = node->left_;
        that_side = node->right_;
    } else {
        // 否则，待检索点落在右侧
        this_side = node->right_;
        that_side = node->left_;
    }

    // 递归调用Knn函数，搜索待检索点的最近邻
    Knn(pt, this_side, knn_result);

    // 检查是否需要展开当前节点的另一侧
    if (NeedExpand(pt, node, knn_result))   // 注意这里是跟当前节点node自己进行比较
        Knn(pt, that_side, knn_result); // 继续递归调用另一侧的近邻搜索算法
}

/**
 * @description: 展开条件判断
 * @param {Vec3f} &pt
 * @param {KdTreeNode} *node
 * @param {priority_queue<NodeAndDistance>} &knn_result
 * @return {*}
 */
bool KdTree::NeedExpand(const Vec3f &pt, KdTreeNode *node, std::priority_queue<NodeAndDistance> &knn_result) const {
    // 判断此时优先队列中的结果是否小于k个
    if (knn_result.size() < k_) 
        // 如果小于k个，则需要继续展开
        return true;

    // 判断是否使用ANN 近似最近邻
    if (approximate_) {
        float d = pt[node->axis_index_] - node->split_thresh_; 
        // 【待检索点与当前节点的分割面的距离】 是否小于 【k近邻优先队列中最大距离的alpha倍】
        if ((d * d) < knn_result.top().distance2_ * alpha_) 
            // 如果小于，则需要继续展开
            return true;
        else 
            // 否则，不需要继续展开
            return false;
    } else {
        // 检测切面距离，看是否有比现在更小的
        float d = pt[node->axis_index_] - node->split_thresh_;
        if ((d * d) < knn_result.top().distance2_) 
            return true;
        else 
            return false;
    }
}

/**
 * @description: 计算叶子节点与待检索点的距离
 * @param {Vec3f} &pt                                       待检索点
 * @param {KdTreeNode} *node                                叶子节点
 * @param {priority_queue<NodeAndDistance>} &knn_result     用于存储最近邻点的优先队列，插入时自动排序，默认降序，队头为最大值
 * @return {*}
 */
void KdTree::ComputeDisForLeaf(const Vec3f &pt, KdTreeNode *node, std::priority_queue<NodeAndDistance> &knn_result) const {
    // 比较与结果队列的差异，如果优于最远距离，则插入

    // 计算待检索点与叶子节点的距离
    float dis2 = Dis2(pt, cloud_[node->point_idx_]);
    if (knn_result.size() < k_) 
        // results中结果 不足 k个，直接插入
        knn_result.push({node, dis2});
    else {
        // results等于k，比较current与max_dis_iter之间的差异
        // 判断待检索点与叶子节点的距离是否小于K近邻优先队列中最大距离
        if (dis2 < knn_result.top().distance2_) {
            // 如果小于，则将该叶子节点插入到优先队列中，然后弹出优先队列中的距离最大的匹配点
            knn_result.push({node, dis2});  // 在优先队列最后插入一个元素
            knn_result.pop();               // 弹出优先队列中的第一个元素（距离最大）
        }
    }
}

/**
 * @description: 计算三个轴上的散布情况
 * @param {IndexVec} &point_idx
 * @param {int} &axis
 * @param {float} &th
 * @param {IndexVec} &left
 * @param {IndexVec} &right
 * @return {*}
 */
bool KdTree::FindSplitAxisAndThresh(const IndexVec &point_idx, int &axis, float &th, IndexVec &left, IndexVec &right) {
    // 计算三个轴上的散布情况，我们使用math_utils.h里的函数
    Vec3f var;
    Vec3f mean;
    // 计算每个轴上的均值和方差
    math::ComputeMeanAndCovDiag(point_idx, mean, var, [this](int idx) { return cloud_[idx]; });
    int max_i, max_j;
    var.maxCoeff(&max_i, &max_j);
    axis = max_i;
    th = mean[axis];

    for (const auto &idx : point_idx) {
        if (cloud_[idx][axis] < th) 
            // 中位数可能向左取整
            left.emplace_back(idx);
        else 
            right.emplace_back(idx);
    }

    // 边界情况检查：输入的points等于同一个值，上面的判定是>=号，所以都进了右侧
    // 这种情况不需要继续展开，直接将当前节点设为叶子就行
    if (point_idx.size() > 1 && (left.empty() || right.empty())) 
        return false;

    return true;
}

/**
 * @description: 重置Kd树
 * @return {*}
 */
void KdTree::Reset() {
    tree_node_id_ = 0;
    root_.reset(new KdTreeNode());
    root_->id_ = tree_node_id_++;
    size_ = 0;
}

/**
 * @description: 清空Kd树
 * @return {*}
 */
void KdTree::Clear() {
    // 遍历所有节点
    for (const auto &np : nodes_) 
        // 判断是否是根节点
        if (np.second != root_.get()) 
            delete np.second; // 删除节点

    nodes_.clear();     // 清空节点
    root_ = nullptr;    // 根节点置空
    size_ = 0;
    tree_node_id_ = 0;
}

/**
 * @description: 打印所有节点
 * @return {*}
 */
void KdTree::PrintAll() {
    // 遍历所有节点
    for (const auto &np : nodes_) {
        auto node = np.second;
        // 判断当前节点左右子节点是否都为空
        if (node->left_ == nullptr && node->right_ == nullptr) 
            // 若为空，则为叶子节点
            LOG(INFO) << "leaf node: " << node->id_ << ", idx: " << node->point_idx_;
        else 
            // 否则，就取当前节点的轴索引和分割阈值
            LOG(INFO) << "node: " << node->id_ << ", axis: " << node->axis_index_ << ", th: " << node->split_thresh_;
    }
}

}  // namespace sad
