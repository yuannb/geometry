#pragma once
#include "globalSymbol.h"
// #include    "../../../3rd/Eigen3/include/eigen3/Eigen/Core"
// #include "../../../3rd/eigen/Eigen/Dense"
#include <Eigen/Dense>
// #include <Eigen/Dense>

class vertex
{
public:
    vertex() = default;
    vertex(Solid *s);
    //have to call it before deconstruct
    bool RemoveListFromSolid(Solid *s);

    ~vertex() { }
public:
    Id          vertexno; // vertex identifier
    HalfEdge    *vedge;   // pointer to halfedge
    // vector      vcoord;   // vertex coordinates
    Vertex      *nextv;   // pointer to next vertex
    Vertex      *prevv;   // pointer to previous vertex
    Eigen::Vector3d vcoord;
};
