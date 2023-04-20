#pragma once
#include "globalSymbol.h"

class loop
{
public:
    loop(Face *f);
    //have to call it before deconstruct
    bool RemoveListFromFace(Face *f);

    ~loop() { }
public:
    HalfEdge    *ledg; // ptr to ring of halfedges
    Face        *lface; // back pointer to face
    Loop        *nextl; // pointer to next loop
    Loop        *prevl; // pointer to previous loop
};
