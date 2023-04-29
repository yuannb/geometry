#include "halfedge.h"
#include "edge.h"

std::shared_ptr<HalfEdge> halfedge::mate()
{
    std::shared_ptr<HalfEdge> thisptr = shared_from_this();
    std::shared_ptr<HalfEdge> he1 = edg->he1.lock();
    if (thisptr != he1)
    {
        return he1;
    }
    return edg->he2.lock();
}