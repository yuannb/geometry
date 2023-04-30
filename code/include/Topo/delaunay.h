#pragma once
#include "geodef.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "memory"

struct triangel
{
    std::vector<std::shared_ptr<Vertex>> vtxarry;
    std::vector<std::vector<int>> face;
};

triangel  delaunay(std::shared_ptr<Solid> s);

triangel  discret(std::shared_ptr<Solid> s);