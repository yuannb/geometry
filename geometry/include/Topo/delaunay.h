#pragma once
#include "geodef.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "memory"
// #include "nurbsCurve.h"

struct triangel
{
    std::vector<std::shared_ptr<Vertex>> vtxarry;
    std::vector<std::vector<int>> face;
};

triangel  delaunay(std::shared_ptr<Solid> s);

triangel  discret(std::shared_ptr<Solid> s);

// std::vector<Eigen::Vector3d> discret1(Bezier<3, 1, false> &curve);
