#pragma once
#include "curve.h"

class line : curve
{
private:
    Eigen::Vector3d startPoint;
    Eigen::Vector3d endPoint;

public:
    line(Eigen::Vector3d s, Eigen::Vector3d e): startPoint(s), endPoint(e) {}
    ~line();
    virtual Eigen::Vector3d get_start_point() { return startPoint; }
    virtual Eigen::Vector3d get_end_point() { return endPoint; }
};

