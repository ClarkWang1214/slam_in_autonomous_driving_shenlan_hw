//
// Created by xiang on 2022/3/23.
//

#include "ch6/occupancy_map.h"
#include "common/eigen_types.h"
#include "common/math_utils.h"

#include <glog/logging.h>
#include <execution>

namespace sad {

/**
 * @description: 占据栅格地图构造函数，用于初始化占据栅格图像
 * @return {*}
 */
OccupancyMap::OccupancyMap() {
    // 预先计算固定大小的模板区域
    BuildModel();
    // 初始化占据栅格图像，大小为1000x1000，每个像素的值为127
    occupancy_grid_ = cv::Mat(image_size_, image_size_, CV_8U, 127);
}

/**
 * @description: 预先计算固定大小的模板区域，利用这个模板来更新栅格信息
 * @return {*}
 */
void OccupancyMap::BuildModel() {
    // 默认模板大小 400x400
    for (int x = -model_size_; x <= model_size_; x++) {
        for (int y = -model_size_; y <= model_size_; y++) {
            Model2DPoint pt;
            pt.dx_ = x;
            pt.dy_ = y;
            // 预先计算每个模板点的距离和夹角
            pt.range_ = sqrt(x * x + y * y) * inv_resolution_;  // 距离
            pt.angle_ = std::atan2(y, x);   // 角度
            pt.angle_ = pt.angle_ > M_PI ? pt.angle_ - 2 * M_PI : pt.angle_;  // limit in 2pi 限制在2pi内
            model_.push_back(pt);
        }
    }
}

/**
 * @description: 计算激光扫描在该角度下的距离值
 * @param {double} angle    角度
 * @param {Ptr} scan        当前帧的激光点云
 * @return {*}
 */
double OccupancyMap::FindRangeInAngle(double angle, Scan2d::Ptr scan) {
    // 将任意角度控制在[-PI, PI]以内
    math::KeepAngleInPI(angle);
    // 判断该角度是否在激光扫描数据的角度范围内
    if (angle < scan->angle_min || angle > scan->angle_max) 
        return 0.0;

    // 计算该角度在激光扫描数据中的索引，即用该角度减去激光扫描数据的最小角度，再除以激光扫描数据的角分辨率
    int angle_index = int((angle - scan->angle_min) / scan->angle_increment);
    // 判断该索引是否越界
    if (angle_index < 0 || angle_index >= scan->ranges.size()) 
        return 0.0; // 若越界，则返回距离值为0

    // 该角度索引值加一
    int angle_index_p = angle_index + 1;
    double real_angle = angle;  // 将该模板点角度设为真实角度初值

    // take range
    double range = 0;
    // 判断该索引的下一个是否越界
    if (angle_index_p >= scan->ranges.size()) 
        range = scan->ranges[angle_index]; // 若越界，说明该索引是最后一个，则距离值为对应的距离值
    else { // 若不是，说明该索引不是最后一个，则需要进行插值
        // 
        double s = ((angle - scan->angle_min) / scan->angle_increment) - angle_index;
        double range1 = scan->ranges[angle_index];      // 当前索引对应的距离值
        double range2 = scan->ranges[angle_index_p];    // 下一个索引对应的距离值
        double real_angle1 = scan->angle_min + scan->angle_increment * angle_index;     // 当前索引对应的真实角度
        double real_angle2 = scan->angle_min + scan->angle_increment * angle_index_p;   // 下一个索引对应的真实角度

        if (range2 < scan->range_min || range2 > scan->range_max) { // 判断下一个索引对应的距离值是否越界
            range = range1;             // 若越界，则距离值为当前索引对应的距离值
            real_angle = real_angle1;   // 将该模板点角度设为当前索引对应的真实角度
        } else if (range1 < scan->range_min || range1 > scan->range_max) {  // 判断当前索引对应的距离值是否越界
            range = range2;             // 若越界，则距离值为下一个索引对应的距离值
            real_angle = real_angle2;   // 将该模板点角度设为下一个索引对应的真实角度
        } else if (std::fabs(range1 - range2) > 0.3) {  // 判断当前索引对应的距离值与下一个索引对应的距离值是否相差较大
            range = s > 0.5 ? range2 : range1;
            real_angle = s > 0.5 ? real_angle2 : real_angle1;
        } else 
            range = range1 * (1 - s) + range2 * s;
    }
    return range;
}

/**
 * @description: 
 * 
 *  激光扫描数据：    局部坐标系，世界坐标系，子地图坐标系，像素坐标系
 *                      scan -> world -> submap -> pixel
 * 
 * @param {shared_ptr<Frame>} frame
 * @param {GridMethod} method
 * @return {*}
 */
void OccupancyMap::AddLidarFrame(std::shared_ptr<Frame> frame, GridMethod method) {
    auto& scan = frame->scan_; // 获取当前帧的激光点云
    // pose_：Tws， submap -> world 子地图与世界系之间的位姿
    // pose_.inverse()：Tsw  world -> submap
    // frame->pose_：Twc，current_frame -> world  每个扫描数据的世界坐标与世界坐标系之间的位姿
    // Tsc = Tws.inv() * Twc = Tsw * Twc  当前激光扫描与子地图之间的位姿关系
    SE2 pose_in_submap = pose_.inverse() * frame->pose_;    // 公式（6.23）Twc = Tws * Tsc
    float theta = pose_in_submap.so2().log(); // 当前帧scan在子地图中位姿对应的角度
    has_outside_pts_ = false;

    // 先计算末端点所在的网格
    std::set<Vec2i, less_vec<2>> endpoints;
    // 遍历激光点云中的每个点
    for (size_t i = 0; i < scan->ranges.size(); ++i) {
        // 跳过无效点
        if (scan->ranges[i] < scan->range_min || scan->ranges[i] > scan->range_max) 
            continue;
        // 计算激光点的真实夹角
        double real_angle = scan->angle_min + i * scan->angle_increment;
        // 极坐标转笛卡尔坐标
        double x = scan->ranges[i] * std::cos(real_angle);
        double y = scan->ranges[i] * std::sin(real_angle);
        // Pwc * pc，激光扫描数据转换到世界系下，再转换为子地图坐标系，最后转换为像素坐标系
        endpoints.emplace(World2Image(frame->pose_ * Vec2d(x, y))); 
    }

    // 根据栅格化方法，更新占据栅格地图
    if (method == GridMethod::MODEL_POINTS) { // 模板化算法
        // 并行遍历所有模板点，生成白色点
        std::for_each(std::execution::par_unseq, 
                      model_.begin(), model_.end(), 
                      [&](const Model2DPoint& pt) {
                        // 将激光扫描数据的世界系坐标 P_w 转换为子地图坐标系下坐标 P_s
                        Vec2i pos_in_image = World2Image(frame->pose_.translation());

                        // 子地图系下的点，加上当前模板点的偏移量，得到当前模板点在子地图系下的坐标
                        Vec2i pw = pos_in_image + Vec2i(pt.dx_, pt.dy_);  

                        // 判断该模板点预先计算的距离值是否小于最小距离阈值
                        if (pt.range_ < closest_th_) {
                            // 若小于，则小距离内认为无物体
                            SetPoint(pw, false);
                            return;
                        }

                        // 模板点预先计算的角度值 减去 当前激光scan在子地图中位姿对应的角度
                        double angle = pt.angle_ - theta;  // Tsc
                        // 计算激光扫描在该角度下的距离值
                        double range = FindRangeInAngle(angle, scan); 
                        // 判断该角度下的距离值是否越界
                        if (range < scan->range_min || range > scan->range_max) {
                            /// 某方向无测量值时，认为无效
                            /// 但离机器比较近时，涂白，表示有物体
                            if (pt.range_ < endpoint_close_th_) 
                                SetPoint(pw, false); 
                            return;
                        }

                        // 将当前模板点的距离值 与 激光在该角度下的距离值进行比较
                        if (range > pt.range_ && endpoints.find(pw) == endpoints.end()) 
                            /// 末端点与车体连线上的点，涂白
                            SetPoint(pw, false);
                    });
    } else {
        Vec2i start = World2Image(frame->pose_.translation());
        std::for_each(  std::execution::par_unseq, endpoints.begin(), endpoints.end(),
                        [this, &start](const auto& pt) { 
                            BresenhamFilling(start, pt); // 调用bresenham算法，对每条激光扫描线计算直线的栅格化
                        }
        );
    }

    // 遍历所有末端点，涂成黑色，表示障碍物
    std::for_each(  endpoints.begin(), endpoints.end(), 
                    [this](const auto& pt) { 
                        SetPoint(pt, true); 
                    }
    );
}

/**
 * @description: 在某个点处填入占据或者非占据信息
 * @param {Vec2i&} pt
 * @param {bool} occupy
 * @return {*}
 */
void OccupancyMap::SetPoint(const Vec2i& pt, bool occupy) {
    int x = pt[0], y = pt[1];
    // 判断该点是否越界
    if (x < 0 || y < 0 || x >= occupancy_grid_.cols || y >= occupancy_grid_.rows) {
        // 外部点是否被占据
        if (occupy) 
            // 若是，说明在栅格化过程中有一个占据栅格落在了外部
            has_outside_pts_ = true; 
        return;
    }

    /// 这里设置了一个上下限，栅格值范围 [117, 137]
    uchar value = occupancy_grid_.at<uchar>(y, x); // 获取该点的占据栅格值
    // 判断该点是否被占据
    if (occupy) {
        // 被占据的情况下，占据值是否大于117
        if (value > 117)  // 下限
            // 若是，占据值减一，多次看到一个栅格是障碍物，占据概率会不断下降
            occupancy_grid_.ptr<uchar>(y)[x] -= 1;
    } else {
        // 未被占据的情况下，占据值是否小于137
        if (value < 137)  // 上限
            // 若是，占据值加一，多次观测到栅格为空，占据概率也会随之上升
            occupancy_grid_.ptr<uchar>(y)[x] += 1;
    }
}

/**
 * @description: 根据占据栅格值，绘制黑白灰三色的占据栅格地图
 * 
 *  原则上，栅格地图是以概率形式进行计算，每个栅格存储被占据的概率，值为0-1之间的浮点数。
 *  当每次被观测到占据或者非占据时，按照概率原理更新自身的占据概率估计。
 *  但这种方式需要引入logit概念，在实际计算中会引入指数和对数的计算，反而会比较复杂。
 * 
 *  在工程上，可以直接使用图像灰度来描述一个格子是否被占据。由于8位图像使用0-255之间的整数来描述灰度。
 *  所以，如果一个栅格数值为127，可以认为该格子状态为“未知”。每次观测到障碍物时，让该格子的数值减一；反之则加一。
 *  
 *  同时限制每个格子的最大，最小计数。相当于统计了每个栅格被观测到障碍物的次数，既能形成更新的效果，也节省了计算时间。这里主要采用这种方式实现概率栅格计算。
 * 
 * @return {*}
 */
cv::Mat OccupancyMap::GetOccupancyGridBlackWhite() const {
    cv::Mat image(image_size_, image_size_, CV_8UC3);
    // 遍历占据栅格地图
    for (int x = 0; x < occupancy_grid_.cols; ++x) {
        for (int y = 0; y < occupancy_grid_.rows; ++y) {
            // 判断占据栅格值是否大于、等于、小于127
            if (occupancy_grid_.at<uchar>(y, x) == 127) 
                // 值等于127，为灰色，未定义区域
                image.at<cv::Vec3b>(y, x) = cv::Vec3b(127, 127, 127);   
            else if (occupancy_grid_.at<uchar>(y, x) < 127) 
                // 值小于127，为黑色，障碍物非可通行区域
                image.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 0, 0);         
            else if (occupancy_grid_.at<uchar>(y, x) > 127) 
                // 值大于127，为白色，无障碍可通行区域
                image.at<cv::Vec3b>(y, x) = cv::Vec3b(255, 255, 255);   
        }
    }
    return image;
}

/**
 * @description: 直接栅格化算法，对每条激光扫描线计算直线的栅格化
 * @param {Vec2i&} p1
 * @param {Vec2i&} p2
 * @return {*}
 */
void OccupancyMap::BresenhamFilling(const Vec2i& p1, const Vec2i& p2) {
    // 计算p1到p2坐标增长的方向
    int dx = p2.x() - p1.x(); 
    int dy = p2.y() - p1.y();
    int ux = dx > 0 ? 1 : -1;
    int uy = dy > 0 ? 1 : -1;

    // 取绝对值
    dx = abs(dx);
    dy = abs(dy);

    // p1坐标
    int x = p1.x();
    int y = p1.y();

    // 比较|dx|和|dy|的大小，取增量大的作为主要增长方向
    if (dx > dy) {
        // 以x轴为增量
        int e = -dx; // 初始误差
        for (int i = 0; i < dx; ++i) {
            x += ux;        // 每当x自增1
            e += 2 * dy;    // 误差增加2dy
            if (e >= 0) {   // 判断误差是否大于零
                y += uy;    // 若是，y自增1
                e -= 2 * dx;// 误差减去2dx
            }
            // 重复上述步骤直到x,y坐标到达p2
            if (Vec2i(x, y) != p2) 
                SetPoint(Vec2i(x, y), false); // 
        }
    } else {
        int e = -dy;
        for (int i = 0; i < dy; ++i) {
            y += uy;
            e += 2 * dx;
            if (e >= 0) {
                x += ux;
                e -= 2 * dy;
            }
            if (Vec2i(x, y) != p2) {
                SetPoint(Vec2i(x, y), false);
            }
        }
    }
}

}  // namespace sad