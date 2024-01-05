#pragma once

#include <memory>
#include <fstream>
#include <algorithm>
#include <iosfwd>
#include "params.h"

class halfedge : public std::enable_shared_from_this<halfedge>
{
public:
    halfedge(){
    };
    ~halfedge() { 
    }
    std::shared_ptr<HalfEdge> mate();

public:
    std::shared_ptr<Edge> edg; // pointer to parent edge
    std::shared_ptr<Vertex> vtx;  // pointer to starting vertex
    std::weak_ptr<Loop> wloop;// back pointer to loop
    std::shared_ptr<HalfEdge> nxt;  // pointer to next halfedge
    std::weak_ptr<HalfEdge> prv; //pointer to previous halfedge
};
