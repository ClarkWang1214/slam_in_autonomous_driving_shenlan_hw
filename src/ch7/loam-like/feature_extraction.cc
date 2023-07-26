//
// Created by xiang on 2022/7/21.
//
#include "ch7/loam-like/feature_extraction.h"
#include <glog/logging.h>


// 地面点云提取相关
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <pcl/ModelCoefficients.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <opencv2/opencv.hpp>

#define FLT_MAX          3.402823466e+38F        // max value
#define M_PI       3.14159265358979323846

namespace sad {

/**
 * @description: 读取原始的Velodyne packets，转换为点云
 * @param {FullCloudPtr} pc_in
 * @param {CloudPtr} pc_out_edge
 * @param {CloudPtr} pc_out_surf
 * @return {*}
 */
void FeatureExtraction::Extract(FullCloudPtr pc_in, CloudPtr pc_out_edge, CloudPtr pc_out_surf, CloudPtr groundCloud_legoloam) {
    int num_scans = 16;     // 16线
    std::vector<CloudPtr> scans_in_each_line;  // 分线数的点云
    CloudPtr full_cloud(new PointCloudType),            // 全点云
             range_cloud(new PointCloudType),           // 距离图像  range image projection
             current_ground(new PointCloudType);        // 当前帧的地面点云【新增】  PCL RANSAC方法
            //  groundCloud_legoloam(new PointCloudType);  // 地面点云【新增】 LegoLOAM方法
    // 初始化16条激光线束的点云
    for (int i = 0; i < num_scans; i++) 
        scans_in_each_line.emplace_back(new PointCloudType); // 16个激光线束点云对象

    // 遍历所有激光点，按照分配到对应的激光线束中
    int k = 0;
    unsigned int rowIdn, columnIdn, index;
    float verticalAngle, horizonAngle, range;
    float azimuth_resolution_deg = 0.3; // 方位角分辨率
    float sensorMinimumRange = 1.0;
    int Horizon_SCAN = int(360 / azimuth_resolution_deg);  // 水平为360度，按分辨率切分即可，360/0.3=1200

    range_cloud->points.resize(num_scans*Horizon_SCAN);
    // std::fill(full_cloud->points.begin(), full_cloud->points.end(), nanPoint);

    // 我们生成一个HSV图像以更好地显示图像
    cv::Mat groundMat(num_scans, Horizon_SCAN, CV_8S, cv::Scalar::all(0));
    for (auto &pt : pc_in->points) {
        // 点云中已经携带了每个点的线束信息，无需计算，直接读取即可。
        assert(pt.ring >= 0 && pt.ring < num_scans);
        // if (rowIdn < 0 || rowIdn >= num_scans) // 0~15
        //         continue;

        // LOG(INFO) << "pt.ring: " << unsigned(pt.ring); // pt.ring 是uint8_t类型，范围是0~255，打印时需要转换为unsigned类型
        PointType p;
        p.x = pt.x; // x,y,z坐标
        p.y = pt.y;
        p.z = pt.z;
        p.intensity = pt.intensity; // intensity存储的是强度
        // 将点按照线束分组，第ring行激光线束的点
        scans_in_each_line[pt.ring]->points.emplace_back(p);
        full_cloud->points.emplace_back(p);  // 用于PCL RANSAC方法

        /// 【LegoLOAM】
        rowIdn = unsigned(pt.ring); // 如果当前点的ring值[0~15]已知，行索引就是ring值

        horizonAngle = atan2(pt.x, pt.y) * 180 / M_PI;
        // LOG(INFO) << "horizonAngle: " << horizonAngle; 

        // 列索引
        columnIdn = -round((horizonAngle-90.0)/azimuth_resolution_deg) + Horizon_SCAN/2;
        // LOG(INFO) << "columnIdn: " << columnIdn; 
        if (columnIdn >= Horizon_SCAN)
            columnIdn -= Horizon_SCAN;

        if (columnIdn < 0 || columnIdn >= Horizon_SCAN)
                continue;

        
        // 计算距离
        range = sqrt(pt.x * pt.x + pt.y * pt.y + pt.z * pt.z);
        if (range < sensorMinimumRange)
                continue;


        // LOG(INFO) << "range: " << range; 
        // 计算出当前点在距离图像中的索引
        index = columnIdn  + rowIdn * Horizon_SCAN;
        // LOG(INFO) << "index: " << index; 
        range_cloud->points[index] = p; // 用于LegoLOAM
        k++;
    }

    // // TODO: 地面点云分割，PCL RANSAC方法 
    // //创建分割时所需要的模型系数对象coefficients及存储内点的点索引集合对象inliers。
    // pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    // pcl::PointIndices::Ptr inliers(new pcl::PointIndices());
    // // 创建分割对象
    // pcl::SACSegmentation<PointType> seg;    // 可选择配置，设置模型系数需要优化
    // seg.setOptimizeCoefficients(true);          // 设置对估计的模型参数进行优化处理
    // seg.setModelType(pcl::SACMODEL_PLANE);      // 设置分割的模型类型为平面
    // seg.setMethodType(pcl::SAC_RANSAC);         // 设置用ransac随机参数估计方法
    // seg.setMaxIterations(100);
    // seg.setDistanceThreshold(0.15);             // 平面厚度距离阈值，表示距离平面多少米之内的点都算平面内点inlier。
    // seg.setInputCloud(full_cloud);                   // 输入点云
    // seg.segment(*inliers, *coefficients);       // 实现分割，并存储分割结果到点集合inliers及存储平面模型系数coefficients
                        
    // pcl::ExtractIndices<PointType> extract;
    // extract.setInputCloud(full_cloud);
    // extract.setIndices(inliers);
    // extract.setNegative(false); // false表示提取地面点
    // extract.filter(*current_ground);
    // // extract.setNegative(true);  // true表示剔除地面点
    // // extract.filter(*current_no_ground);

    // LOG(INFO) << "ground RANSAC: " << current_ground->size();
    // pcl::io::savePCDFileBinaryCompressed("./data/ch7/ground_RANSAC.pcd", *current_ground);


    // TODO: 地面点云分割，legoloam 
    // legoloam: 通过相邻两线的点云来计算水平夹角，水平夹角在特定范围（10度）内为地面点云。

    size_t lowerInd, upperInd;
    float diffX, diffY, diffZ, angle;
    
    float sensorMountAngle = 0.0;
    // groundMat
    // -1, no valid info to check if ground of not  无效点
    //  0, initial value, after validation, means not ground  表示不是地面点
    //  1, ground 地面点标记
    for (size_t j = 0; j < Horizon_SCAN; ++j){ // 一根线束上共1800个点
        for (size_t i = 0; i < 7; ++i){ // 对于VLP16线来说，下面7根为地面线束，[0-6]

            lowerInd = j + ( i )*Horizon_SCAN;  // 第 i   根线束的第j个点的索引
            upperInd = j + (i+1)*Horizon_SCAN;  // 第 i+1 根线束的第j个点的索引

            if (range_cloud->points[lowerInd].intensity == -1 ||
                range_cloud->points[upperInd].intensity == -1){
                // no info to check, invalid points
                groundMat.at<int8_t>(i,j) = -1; // 无效点
                continue;
            }
                
            // 上下两根激光线束扫到的垂直相邻点的坐标差
            diffX = range_cloud->points[upperInd].x - range_cloud->points[lowerInd].x;
            diffY = range_cloud->points[upperInd].y - range_cloud->points[lowerInd].y;
            diffZ = range_cloud->points[upperInd].z - range_cloud->points[lowerInd].z;

            // 计算相邻两根线束的点云来计算水平夹角，水平夹角在特定范围内，则认为是地面点
            angle = atan2(diffZ, sqrt(diffX*diffX + diffY*diffY) ) * 180 / M_PI;
            if (abs(angle - sensorMountAngle) <= 10){
                // 那么，这两根相邻线束上相同索引处扫到点，就是地面点，标记为1
                groundMat.at<int8_t>(i,j) = 1;      // 第i根线束的第j个点是地面点
                groundMat.at<int8_t>(i+1,j) = 1;    // 第i+1根线束的第j个点是地面点
            }
        }
    }

    // extract ground cloud (groundMat == 1) 
    // mark entry that doesn't need to label (ground and invalid point) for segmentation
    // note that ground remove is from 0~N_SCAN-1, need rangeMat for mark label matrix for the 16th scan
    // 遍历16线激光线束的底下的7根线束，将这7根线束上标记为1的地面点提取出来，放到groundCloud中
    for (size_t i = 0; i <= 7; ++i){    // 7根地面线束
        for (size_t j = 0; j < Horizon_SCAN; ++j){  // 1800个点
            // 判断第i根线束上的第j个点是否是地面点
            if (groundMat.at<int8_t>(i,j) == 1)
                // 如果是地面点，则放到groundCloud中
                groundCloud_legoloam->push_back(range_cloud->points[j + i*Horizon_SCAN]);
        }
    }
    // pcl::io::savePCDFileBinaryCompressed("./data/ch7/full_cloud.pcd", *full_cloud);
    
    // pcl::io::savePCDFileBinaryCompressed("./data/ch7/ground_LegoLOAM.pcd", *groundCloud_legoloam);

    LOG(INFO) << "ground LegoLOAM: " << groundCloud_legoloam->size();
    
    // 处理曲率
    for (int i = 0; i < num_scans; i++) {
        LOG(INFO) << "scans_in_each_line[i]->points.size(): " << scans_in_each_line[i]->points.size();
        // 判断每条线束中的点数是否大于131
        if (scans_in_each_line[i]->points.size() < 131)  
            continue;   // 小于131的线束，直接跳过

        std::vector<IdAndValue> cloud_curvature;  // 存储每条线对应的曲率

        // 当前激光线束上总的点数
        int total_points = scans_in_each_line[i]->points.size() - 10;

        // 依次计算每个点的曲率，两头留一定余量
        for (int j = 5; j < (int)scans_in_each_line[i]->points.size() - 5; j++) {
            // 分别计算当前点与前后各5个点之间的距离值
            double diffX = scans_in_each_line[i]->points[j - 5].x + scans_in_each_line[i]->points[j - 4].x +
                           scans_in_each_line[i]->points[j - 3].x + scans_in_each_line[i]->points[j - 2].x +
                           scans_in_each_line[i]->points[j - 1].x - 10 * scans_in_each_line[i]->points[j].x +
                           scans_in_each_line[i]->points[j + 1].x + scans_in_each_line[i]->points[j + 2].x +
                           scans_in_each_line[i]->points[j + 3].x + scans_in_each_line[i]->points[j + 4].x +
                           scans_in_each_line[i]->points[j + 5].x;
            double diffY = scans_in_each_line[i]->points[j - 5].y + scans_in_each_line[i]->points[j - 4].y +
                           scans_in_each_line[i]->points[j - 3].y + scans_in_each_line[i]->points[j - 2].y +
                           scans_in_each_line[i]->points[j - 1].y - 10 * scans_in_each_line[i]->points[j].y +
                           scans_in_each_line[i]->points[j + 1].y + scans_in_each_line[i]->points[j + 2].y +
                           scans_in_each_line[i]->points[j + 3].y + scans_in_each_line[i]->points[j + 4].y +
                           scans_in_each_line[i]->points[j + 5].y;
            double diffZ = scans_in_each_line[i]->points[j - 5].z + scans_in_each_line[i]->points[j - 4].z +
                           scans_in_each_line[i]->points[j - 3].z + scans_in_each_line[i]->points[j - 2].z +
                           scans_in_each_line[i]->points[j - 1].z - 10 * scans_in_each_line[i]->points[j].z +
                           scans_in_each_line[i]->points[j + 1].z + scans_in_each_line[i]->points[j + 2].z +
                           scans_in_each_line[i]->points[j + 3].z + scans_in_each_line[i]->points[j + 4].z +
                           scans_in_each_line[i]->points[j + 5].z;
            // 第j个点，对应的曲率值，构造结构体对象
            IdAndValue distance(j, diffX * diffX + diffY * diffY + diffZ * diffZ);
            cloud_curvature.push_back(distance);    // 存入点云曲率容器中
        }

        // 对每个区间选取特征，把360度分为6个区间
        for (int j = 0; j < 6; j++) {
            int sector_length = (int)(total_points / 6);
            int sector_start = sector_length * j;
            int sector_end = sector_length * (j + 1) - 1;
            if (j == 5) 
                sector_end = total_points - 1;
            // 根据索引值，截取曲率容器中的一段区间
            std::vector<IdAndValue> sub_cloud_curvature(cloud_curvature.begin() + sector_start, cloud_curvature.begin() + sector_end);

            /// 从每个区间中提取角点和平面点
            ExtractFromSector(scans_in_each_line[i], sub_cloud_curvature, pc_out_edge, pc_out_surf);
        }
    }
}

/**
 * @description: 对单独一段区域提取角点和面点
 * @param {CloudPtr} &pc_in
 * @param {vector<IdAndValue>} &cloud_curvature
 * @param {CloudPtr} &pc_out_edge
 * @param {CloudPtr} &pc_out_surf
 * @return {*}
 */
void FeatureExtraction::ExtractFromSector(const CloudPtr &pc_in, std::vector<IdAndValue> &cloud_curvature, CloudPtr &pc_out_edge, CloudPtr &pc_out_surf) {
    // 按曲率排序
    std::sort(cloud_curvature.begin(), cloud_curvature.end(),
              [](const IdAndValue &a, const IdAndValue &b) { return a.value_ < b.value_; });

    int largest_picked_num = 0;
    int point_info_count = 0;

    /// 按照曲率最大的开始搜，选取曲率最大的角点
    std::vector<int> picked_points;  // 标记被选中的角点，角点附近的点都不会被选取
    for (int i = cloud_curvature.size() - 1; i >= 0; i--) {
        int ind = cloud_curvature[i].id_;
        if (std::find(picked_points.begin(), picked_points.end(), ind) == picked_points.end()) {
            if (cloud_curvature[i].value_ <= 0.1) 
                break;

            largest_picked_num++;
            picked_points.push_back(ind);

            if (largest_picked_num <= 20) {
                pc_out_edge->push_back(pc_in->points[ind]);
                point_info_count++;
            } else 
                break;

            for (int k = 1; k <= 5; k++) {
                double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k - 1].x;
                double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k - 1].y;
                double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k - 1].z;
                if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05) 
                    break;
                picked_points.push_back(ind + k);
            }
            for (int k = -1; k >= -5; k--) {
                double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k + 1].x;
                double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k + 1].y;
                double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k + 1].z;
                if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05) 
                    break;
                picked_points.push_back(ind + k);
            }
        }
    }

    /// 选取曲率较小的平面点
    for (int i = 0; i <= (int)cloud_curvature.size() - 1; i++) {
        int ind = cloud_curvature[i].id_;
        if (std::find(picked_points.begin(), picked_points.end(), ind) == picked_points.end()) 
            pc_out_surf->push_back(pc_in->points[ind]);
    }
}

}  // namespace sad