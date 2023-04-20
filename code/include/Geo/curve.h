#pragma once
// #include "../../../3rd/Eigen3/include/eigen3/Eigen/Dense"
// #include "../../../3rd/eigen/Eigen/Dense"
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