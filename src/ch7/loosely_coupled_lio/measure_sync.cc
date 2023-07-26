//
// Created by xiang on 23-2-2.
//
#include "measure_sync.h"

namespace sad {

/**
 * @description: 激光雷达数据与IMU数据同步
 *              将一帧激光以及这一帧点云扫描scan期间的IMU数据一起打包到MeasureGroup中
 * @return {*}
 */
bool MessageSync::Sync() {
    // 如果激光雷达数据或者IMU数据为空，则返回false
    if (lidar_buffer_.empty() || imu_buffer_.empty()) 
        return false;

    // 判断激光雷达数据是否已经推送到队列中
    if (!lidar_pushed_) {
        // 若还没有，
        // 则从雷达数据队列中取出第一个数据，放入measures_中
        measures_.lidar_ = lidar_buffer_.front();               
        measures_.lidar_begin_time_ = time_buffer_.front();     // 记录激光雷达数据的开始时间

        // 记录当前帧雷达数据的结束时间，该帧雷达数据开始时间 + 最后一个点的时间，再单位转为秒
        lidar_end_time_ = measures_.lidar_begin_time_ + measures_.lidar_->points.back().time / double(1000);
        measures_.lidar_end_time_ = lidar_end_time_;

        // 标记激光雷达数据已经推送到队列中
        lidar_pushed_ = true;
    }

    // 判断IMU队列中的最后一个IMU数据的时间戳 是否早于 雷达数据结束的时间，说明IMU数据不够用
    if (last_timestamp_imu_ < lidar_end_time_) 
        return false; // 若是，则返回，同步失败

    // 获取IMU数据缓存队列中的第一个IMU数据的时间戳
    double imu_time = imu_buffer_.front()->timestamp_;

    // 清空此时激光雷达数据到来时的IMU测量数据队列
    measures_.imu_.clear();

    // 遍历IMU数据缓存队列，将时间戳早于激光雷达数据结束时间的IMU数据放入measures_中
    while ((!imu_buffer_.empty()) && (imu_time < lidar_end_time_)) {
        // 获取此时IMU数据缓存队列中的第一个IMU数据的时间戳
        imu_time = imu_buffer_.front()->timestamp_;

        // 判断时间戳是否晚于激光雷达最后一个点的时间戳
        if (imu_time > lidar_end_time_) 
            break; // 若是，说明此时的IMU数据已经超过了激光雷达数据的时间范围，退出while循环

        // 将此时IMU数据缓存队列中的第一个IMU数据放入measures_中
        measures_.imu_.push_back(imu_buffer_.front());

        // 弹出队列中的第一个IMU数据
        imu_buffer_.pop_front();
    }

    // 弹出队列中的第一个激光雷达数据、时间戳
    lidar_buffer_.pop_front();
    time_buffer_.pop_front();

    // 激光雷达数据已被弹出，标记此时队列中没有数据
    lidar_pushed_ = false;

    // 判断是否有回调函数
    if (callback_) 
        callback_(measures_);   // 将激光雷达到来期间的IMU数据传入回调函数处理

    return true;
}

void MessageSync::Init(const std::string& yaml) { conv_->LoadFromYAML(yaml); }

}  // namespace sad