#pragma once
// #include "../../../3rd/Eigen3/include/eigen3/Eigen/Dense"
// #include "../../../3rd/eigen/Eigen/Dense"
#include <Eigen/Dense>
// #include <Eigen/Dense>
#include "surface.h"

class plane : surface
{
private:
    Eigen::Vector4d equation;
    // Eigen::Vector3d udir;
    // Eigen::Vector3d vdir;
    // Eigen::Vector3d zdir;
    // Eigen::Vector3d origin;

public:
    plane(Eigen::Vector4d equ);
    ~plane() { };
};
