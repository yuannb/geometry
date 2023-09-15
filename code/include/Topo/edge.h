#pragma once

#include <memory>
#include "params.h"

class old_curve;

// class Solid;

class edge : public std::enable_shared_from_this<edge>
{
public: 
    edge();
    bool RemoveListFromSolid(std::shared_ptr<Solid> s);

    ~edge();
public:
    std::weak_ptr<HalfEdge> he1;  // pointer to right halfedge
    std::weak_ptr<HalfEdge> he2;  // pointer to left halfedge
    std::shared_ptr<Edge>   nexte;// pointer to next edge
    std::weak_ptr<Edge>     preve;// pointer to previous edge
    std::shared_ptr<old_curve>    cur;
};
