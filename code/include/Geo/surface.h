
#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
class surface
{
private:
public:
    surface();
    ~surface();
    virtual Eigen::Vector3d get_normal() = 0;
};