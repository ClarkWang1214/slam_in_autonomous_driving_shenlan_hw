//
// Created by xiang on 22-12-20.
//

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>

#include "common/eigen_types.h"

#include "common/point_cloud_utils.h"
#include "keyframe.h"

#include "ch7/ndt_3d.h"
#include "common/math_utils.h"
#include <Eigen/SVD>
#include <execution>

DEFINE_string(map_path, "./data/ch9/", "导出数据的目录");
DEFINE_double(voxel_size, 0.1, "导出地图分辨率");

using namespace sad;

Ndt3d::Options options_;

using KeyType = Eigen::Matrix<int, 3, 1>;  // 体素的索引

// 【新增】将局部地图区块内的NDT体素保存到指定txt文件中
void SaveNDTVoxelToFile(const std::string &filePath, std::unordered_map<KeyType, Ndt3d::VoxelData, hash_vec<3>> grids) {
    // std::ofstream fout_ndt_map(filePath, std::ios::binary);
    std::ofstream fout_ndt_map(filePath);
    
    if (fout_ndt_map.is_open()) {
        // 遍历所有体素栅格
        for (const auto& voxel : grids) {
            // 将矩阵的数据拷贝到行向量中
            Eigen::VectorXd rowVector(9);
            Mat3d mat = voxel.second.info_;
            Eigen::Map<Eigen::VectorXd>(rowVector.data(), rowVector.size()) = Eigen::Map<Eigen::VectorXd>(mat.data(), mat.size());
            // 将体素索引、均值和协方差矩阵写入文件的一行
            fout_ndt_map << voxel.first.transpose() << " " << voxel.second.mu_.transpose() << " " << rowVector.transpose() << std::endl;
        }
    } 
    else 
        LOG(INFO) << "Failed to open the file.";
}

/**
 * @description: 构建NDT体素
 * @param target_ 目标点云
 * @return {*}
 */
void BuildNdtMapVoxels(Vec2i submap_key, CloudPtr target_, std::ofstream fout_ndt_map, std::unordered_map<KeyType, Ndt3d::VoxelData, hash_vec<3>>& grids_, double voxel_size) {
    assert(target_ != nullptr); // 目标点云指针不能为空
    assert(target_->empty() == false);  // 目标点云不能为空
    grids_.clear(); // 清空体素栅格

    // 【新增】对目标点云进行体素滤波，用NDT的体素分辨率作为滤波分辨率
    sad::VoxelGrid(target_, voxel_size);

    double inv_voxel_size = 1.0 / voxel_size;

    /// 分配体素索引
    std::vector<size_t> index(target_->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });

    // 生成体素栅格
    std::for_each(index.begin(), index.end(), [&](const size_t& idx) {
        auto pt = ToVec3d(target_->points[idx]);
        // 对目标点云中的每个点，计算它所在的体素栅格ID
        auto key = (pt * inv_voxel_size).cast<int>(); 

        // 查看该栅格是否已存在
        if (grids_.find(key) == grids_.end()) 
            grids_.insert({key, {idx}}); // 若不存在，则插入该栅格
        else 
            grids_[key].idx_.emplace_back(idx); // 若存在，则将该点的索引插入到该体素栅格中
    });

    // 并发遍历所有体素栅格
    std::for_each(std::execution::par_unseq, 
                grids_.begin(), grids_.end(), 
                [&](auto& v) {
                    // 判断体素中的点数是否大于阈值3个
                    if (v.second.idx_.size() > 3) {
                        // 要求至少有３个点，才会计算每个体素中的均值和协方差
                        math::ComputeMeanAndCov(v.second.idx_, 
                                                v.second.mu_, v.second.sigma_,
                                                [&](const size_t& idx) { 
                                                    return ToVec3d(target_->points[idx]); 
                                                });
                        // SVD 检查最大与最小奇异值，限制最小奇异值
                        Eigen::JacobiSVD svd(v.second.sigma_, Eigen::ComputeFullU | Eigen::ComputeFullV);
                        Vec3d lambda = svd.singularValues();
                        if (lambda[1] < lambda[0] * 1e-3) 
                            lambda[1] = lambda[0] * 1e-3;

                        if (lambda[2] < lambda[0] * 1e-3) 
                            lambda[2] = lambda[0] * 1e-3;

                        Mat3d inv_lambda = Vec3d(1.0 / lambda[0], 1.0 / lambda[1], 1.0 / lambda[2]).asDiagonal();
                        v.second.info_ = svd.matrixV() * inv_lambda * svd.matrixU().transpose(); // 信息矩阵
                    }
                });

    // 遍历所有体素栅格
    for (auto iter = grids_.begin(); iter != grids_.end();) {
        if (iter->second.idx_.size() > 3) 
            iter++;
        else 
            iter = grids_.erase(iter); // 删除点数不够3个的体素栅格
    }

    if(grids_.size() != 0){
        fout_ndt_map << submap_key[0] << " " << submap_key[1] << std::endl;
        // 将每个地图区块（xy方向上：100x100 m^2）内的所有NDT体素块保存逐行保存在txt文件中，后续修改为二进制保存减少文件内存大小
        SaveNDTVoxelToFile("./data/ch9/ndt_map_data_" + std::to_string(int(voxel_size)) + "/" + std::to_string(submap_key[0]) + "_" + std::to_string(submap_key[1]) + ".txt", grids_);
    }
}

/**
 * @description: 地图的导出
 * 
 *  最后，需要将整个点云地图进行导出。将所有点云放入一个地图固然可以方便查看，但在实时定位系统中这样做不是很好。
 *  在多数应用中，希望控制实时点云的载入规模，比如，只加载自身周围200米内的点云，其它范围的点云则视情况卸载掉，
 *  这样可以控制实时系统的计算量。下面重新组织本章建立的点云，按照100米的变成进行分块，并设计分块加载和卸载接口。
 *  
 *  点云的切分，实际上根据点的坐标计算所在的网格，然后将它投到对应的网格中去。
 * 
 * 
 * @param {int} argc
 * @param {char**} argv
 * @return {*}
 */
int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);

    using namespace sad;

    // 加载所有关键帧点云
    std::map<IdType, KFPtr> keyframes;
    if (!LoadKeyFrames("./data/ch9/keyframes.txt", keyframes)) {
        LOG(ERROR) << "failed to load keyframes";
        return 0;
    }

    std::map<Vec2i, CloudPtr, less_vec<2>> pts_map_data;  // 以网格ID为索引的地图区块点云数据
    
    float resolution = FLAGS_voxel_size; // 地图分辨率，默认为0.1米，即10cm

    // pcl::VoxelGrid<PointType> voxel_grid_filter;
    // voxel_grid_filter.setLeafSize(resolution, resolution, resolution);

    // 逻辑和dump map差不多，但每个点查找它的网格ID，没有的话会创建
    // 遍历所有关键帧点云
    for (auto& kfp : keyframes) {
        auto kf = kfp.second; // 取出关键帧

        // 加载以当前关键帧kf索引id命名的点云
        kf->LoadScan("./data/ch9/");

        // 存储世界系下点云
        CloudPtr cloud_trans(new PointCloudType);
        // 将关键帧点云转换到世界坐标系下，使用第二轮优化后的位姿
        pcl::transformPointCloud(*kf->cloud_, *cloud_trans, kf->opti_pose_2_.matrix());

        // 体素滤波
        // CloudPtr kf_cloud_voxeled(new PointCloudType);
        // voxel_grid_filter.setInputCloud(cloud_trans);
        // voxel_grid_filter.filter(*kf_cloud_voxeled);
        sad::VoxelGrid(cloud_trans, resolution);

        LOG(INFO) << "building kf " << kf->id_ << " in " << keyframes.size();

        // add to grid
        // 遍历当前关键帧滤波后点云中的每个点
        for (const auto& pt : cloud_trans->points) { // kf_cloud_voxeled->points
            // 计算每个点所在的网格ID
            int gx = floor((pt.x - 50.0) / 100); // 以100米为一个网格，计算x方向上的网格ID
            int gy = floor((pt.y - 50.0) / 100); // 以100米为一个网格，计算y方向上的网格ID
            Vec2i key(gx, gy);

            // auto key_ndtvoxel = (pt * )

            // 在pts_map_data中查找是否存在该网格ID
            auto iter = pts_map_data.find(key);
            if (iter == pts_map_data.end()) {
                // 若不存在，则创建一个新的点云
                // create point cloud
                CloudPtr cloud(new PointCloudType);
                cloud->points.emplace_back(pt); // 将该点加入到点云中
                cloud->is_dense = false;
                cloud->height = 1;
                pts_map_data.emplace(key, cloud); // 然后将该点云加入到pts_map_data中
            } else 
                // 若存在，则将该点加入到该网格的点云中
                iter->second->points.emplace_back(pt); 
        }
    }

    // 存储点云和索引文件
    LOG(INFO) << "saving maps, grids: " << pts_map_data.size();
    std::system("mkdir -p ./data/ch9/pts_map_data/");
    std::system("rm -rf ./data/ch9/pts_map_data/*");  // 清理一下文件夹
    std::ofstream fout_pts_map("./data/ch9/pts_map_data/pts_map_index.txt");

    // 【新增】存储NDT体素均值与协方差，以及索引文件
    std::system("mkdir -p ./data/ch9/ndt_map_data_10/");
    std::system("mkdir -p ./data/ch9/ndt_map_data_5/");
    std::system("mkdir -p ./data/ch9/ndt_map_data_4/");
    std::system("mkdir -p ./data/ch9/ndt_map_data_3/");
    std::system("mkdir -p ./data/ch9/ndt_map_data_1/");
    
    // 清理一下文件夹
    std::system("rm -rf ./data/ch9/ndt_map_data_10/*");  
    std::system("rm -rf ./data/ch9/ndt_map_data_5/*");  
    std::system("rm -rf ./data/ch9/ndt_map_data_4/*");
    std::system("rm -rf ./data/ch9/ndt_map_data_3/*");  
    std::system("rm -rf ./data/ch9/ndt_map_data_1/*"); 

    std::ofstream fout_ndt_map_10("./data/ch9/ndt_map_data_10/ndt_map_index.txt");
    std::ofstream fout_ndt_map_5("./data/ch9/ndt_map_data_5/ndt_map_index.txt");
    std::ofstream fout_ndt_map_4("./data/ch9/ndt_map_data_4/ndt_map_index.txt");
    std::ofstream fout_ndt_map_3("./data/ch9/ndt_map_data_3/ndt_map_index.txt");
    std::ofstream fout_ndt_map_1("./data/ch9/ndt_map_data_1/ndt_map_index.txt");

    using NdtmapVoxel = std::unordered_map<KeyType, Ndt3d::VoxelData, hash_vec<3>>;

    //【新增】grids中保存了当前地图区块下所有体素格子的均值和协方差，需要将体素所在格子的三维空间坐标
    NdtmapVoxel grids_10, grids_5, grids_4, grids_3, grids_1;

    std::map<Vec2i, NdtmapVoxel, less_vec<2>> ndt_map_data;  // 以网格ID为索引的地图数据

    // 遍历所有地图区块
    for (auto& dp : pts_map_data) {
        // 将索引文件以文本格式存储在map_index.txt中，方便快速读取到地图区块的位置
        fout_pts_map << dp.first[0] << " " << dp.first[1] << std::endl;
        
        // 获取该地图区块的点云数量
        dp.second->width = dp.second->size();

        // 函数内部通过swap交换指针，将滤波后的点云赋值回dp.second
        sad::VoxelGrid(dp.second, 0.1); // 体素滤波分辨率为0.1米

        // 将该地图区块的点云保存到指定pcd文件中
        sad::SaveCloudToFile("./data/ch9/pts_map_data/" + std::to_string(dp.first[0]) + "_" + std::to_string(dp.first[1]) + ".pcd", *dp.second);
        
        // // 容器清空
        // grids_10.clear();
        // grids_5.clear();
        // grids_4.clear();
        // grids_3.clear();
        // grids_1.clear();

        
        // 构建五种不同分辨率的NDT体素
        BuildNdtMapVoxels(dp.first, dp.second, fout_ndt_map_10, grids_10, 10.0);   // 体素分辨率为10米
        BuildNdtMapVoxels(dp.first, dp.second, fout_ndt_map_5, grids_5, 5.0);      // 体素分辨率为5米 
        BuildNdtMapVoxels(dp.first, dp.second, fout_ndt_map_4, grids_4, 4.0);      // 体素分辨率为4米
        BuildNdtMapVoxels(dp.first, dp.second, fout_ndt_map_3, grids_3, 3.0);      // 体素分辨率为3米
        BuildNdtMapVoxels(dp.first, dp.second, fout_ndt_map_1, grids_1, 1.0);      // 体素分辨率为1米

        // if(grids_10.size() != 0){
        //     fout_ndt_map_10 << dp.first[0] << " " << dp.first[1] << std::endl;
        //     // 将每个地图区块（xy方向上：100x100 m^2）内的所有NDT体素块（分辨率默认为 10x10x10 m^3）保存逐行保存在txt文件中，后续修改为二进制保存减少文件内存大小
        //     SaveNDTVoxelToFile("./data/ch9/ndt_map_data_10/" + std::to_string(dp.first[0]) + "_" + std::to_string(dp.first[1]) + ".txt", grids_10);
        // }

        // if(grids_5.size() != 0){
        //     fout_ndt_map_5 << dp.first[0] << " " << dp.first[1] << std::endl;
        //     // 将每个地图区块（xy方向上：100x100 m^2）内的所有NDT体素块（分辨率默认为 5x5x5 m^3）保存逐行保存在txt文件中，后续修改为二进制保存减少文件内存大小
        //     SaveNDTVoxelToFile("./data/ch9/ndt_map_data_5/" + std::to_string(dp.first[0]) + "_" + std::to_string(dp.first[1]) + ".txt", grids_5);
        // }

        // if(grids_4.size() != 0) {
        //     fout_ndt_map_4 << dp.first[0] << " " << dp.first[1] << std::endl;
        //     // 将每个地图区块（xy方向上：100x100 m^2）内的所有NDT体素块（分辨率默认为 4x4x4 m^3）保存逐行保存在txt文件中，后续修改为二进制保存减少文件内存大小
        //     SaveNDTVoxelToFile("./data/ch9/ndt_map_data_4/" + std::to_string(dp.first[0]) + "_" + std::to_string(dp.first[1]) + ".txt", grids_4);
        // }

        // if(grids_3.size() != 0) {
        //     fout_ndt_map_3 << dp.first[0] << " " << dp.first[1] << std::endl;
        //     // 将每个地图区块（xy方向上：100x100 m^2）内的所有NDT体素块（分辨率默认为 3x3x3 m^3）保存逐行保存在txt文件中，后续修改为二进制保存减少文件内存大小
        //     SaveNDTVoxelToFile("./data/ch9/ndt_map_data_3/" + std::to_string(dp.first[0]) + "_" + std::to_string(dp.first[1]) + ".txt", grids_3);
        // }

        // if(grids_1.size() != 0) {
        //     fout_ndt_map_1 << dp.first[0] << " " << dp.first[1] << std::endl;
        //     // 将每个地图区块（xy方向上：100x100 m^2）内的所有NDT体素块（分辨率默认为 1x1x1 m^3）保存逐行保存在txt文件中，后续修改为二进制保存减少文件内存大小
        //     SaveNDTVoxelToFile("./data/ch9/ndt_map_data_1/" + std::to_string(dp.first[0]) + "_" + std::to_string(dp.first[1]) + ".txt", grids_1);
        // }
    }
    fout_pts_map.close();
    fout_ndt_map_10.close();
    fout_ndt_map_5.close();
    fout_ndt_map_4.close();
    fout_ndt_map_3.close();
    fout_ndt_map_1.close();

    LOG(INFO) << "done.";
    return 0;
}
