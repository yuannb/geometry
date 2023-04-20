#pragma once
#include "globalSymbol.h"

class curve;

// class Solid;

class edge
{
public: 
    edge(Solid *s);

    //have to call it before deconstruct
    bool RemoveListFromSolid(Solid *s);
    
    ~edge();
public:
    HalfEdge    *he1;  // pointer to right halfedge
    HalfEdge    *he2;  // pointer to left halfedge
    Edge        *nexte;// pointer to next edge
    Edge        *preve;// pointer to previous edge
    curve       *cur;
};
