#pragma once
#include "geodef.h"
#include <vector>
#include <unordered_map>
#include <algorithm>

struct triangel
{
    std::vector<Vertex*> vtxarry;
    std::vector<std::vector<int>> face;
};

triangel  delaunay(Solid *s);