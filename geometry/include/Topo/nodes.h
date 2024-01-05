#pragma once
#include    "solid.h"
#include    "face.h"
#include    "loop.h"
#include    "halfedge.h"
#include    "halfedge.h"
#include    "vertex.h"
#include    "edge.h"


union nodes
{
    Solid       s;
    Face        f;
    Loop        l;
    HalfEdge    h;
    Vertex      v;
    Edge        e;
};
