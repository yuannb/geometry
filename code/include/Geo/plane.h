#pragma once
// #include "../../../3rd/Eigen3/include/eigen3/Eigen/Dense"
// #include "../../../3rd/eigen/Eigen/Dense"
#include <Eigen/Dense>
// #include <Eigen/Dense>
#include "surface.h"

class plane : public surface
{
public:
    Eigen::Vector4d equation;
    // Eigen::Vector3d udir;
    // Eigen::Vector3d vdir;
    // Eigen::Vector3d zdir;
    // Eigen::Vector3d origin;
    void set_eqution(Eigen::Vector4d vec) { equation = vec; }
    Eigen::Vector3d virtual get_normal() { return Eigen::Vector3d(equation[0], equation[1], equation[2]); }

public:
    plane(Eigen::Vector4d equ);
    ~plane() { };
};
