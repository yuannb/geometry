#pragma once
#include "globalSymbol.h"

class halfedge
{
public:
    halfedge() = default;
    ~halfedge() { }

public:
    Edge        *edg; // pointer to parent edge
    Vertex      *vtx;  // pointer to starting vertex
    Loop        *wloop;// back pointer to loop
    HalfEdge    *nxt;  // pointer to next halfedge
    HalfEdge    *prv;  // pointer to previous halfedge
};
