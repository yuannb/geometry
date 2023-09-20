#pragma once
#include <Eigen/Dense>
#include "surface.h"

class plane : public old_surface
{
public:
    Eigen::Vector4d equation;
    void set_eqution(Eigen::Vector4d vec) { equation = vec; }
    Eigen::Vector3d virtual get_normal() { return Eigen::Vector3d(equation[0], equation[1], equation[2]); }

public:
    plane(Eigen::Vector4d equ);
    ~plane() { };
};
