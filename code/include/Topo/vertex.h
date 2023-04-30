#pragma once

#include <Eigen/Dense>
#include <memory>
        #include <fstream>
#include <algorithm>
#include <iosfwd>

class vertex : public std::enable_shared_from_this<vertex>
{
public:
    vertex()
    {
        vedge.reset();
        nextv = nullptr;
        prevv.reset();
    };
    vertex(double x, double y, double z) : vcoord(x, y, z) {        
        vedge.reset();
        nextv = nullptr;
        prevv.reset();}
    bool RemoveListFromSolid(std::shared_ptr<Solid> s);

    ~vertex();
public:
    Id          vertexno; // vertex identifier
    std::weak_ptr<HalfEdge> vedge;   // pointer to halfedge
    std::shared_ptr<Vertex> nextv;   // pointer to next vertex
    std::weak_ptr<Vertex> prevv;   // pointer to previous vertex
    Eigen::Vector3d vcoord; // vertex coordinates
};
