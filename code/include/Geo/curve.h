#pragma once
#include <Eigen/Dense>

class curve
{
private:
    /* data */
public:
    curve(/* args */);
    ~curve();
    virtual Eigen::Vector3d get_start_point() = 0;
    virtual Eigen::Vector3d get_end_point() = 0;
};