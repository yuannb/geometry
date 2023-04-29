#pragma once
#include "globalSymbol.h"
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
    };
    vertex(double x, double y, double z) : vcoord(x, y, z) {}
    // vertex(std::shared_ptr<Solid> s);
    //have to call it before deconstruct
    bool RemoveListFromSolid(std::shared_ptr<Solid> s);

    ~vertex() { }
public:
    Id          vertexno; // vertex identifier
    std::weak_ptr<HalfEdge> vedge;   // pointer to halfedge
    // HalfEdge    *vedge;   // pointer to halfedge
    std::shared_ptr<Vertex> nextv;   // pointer to next vertex
    // Vertex      *nextv;   // pointer to next vertex
    std::weak_ptr<Vertex> prevv;   // pointer to previous vertex
    // Vertex      *prevv;   // pointer to previous vertex
    Eigen::Vector3d vcoord; // vertex coordinates
};
