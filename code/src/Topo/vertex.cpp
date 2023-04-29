#include "vertex.h"
#include "solid.h"

// vertex::vertex(std::shared_ptr<Solid> s)
// {
//     nextv = s->svertes;
//     prevv.reset();
//     std::shared_ptr<Vertex> thisv = shared_from_this();
//     if (s->svertes)
//         s->svertes->prevv = thisv;
//     s->svertes = thisv;
//     this->vedge.reset();
// }

bool vertex::RemoveListFromSolid(std::shared_ptr<Solid> s)
{
    if (this == s->svertes.get())
    {
        s->svertes = s->svertes->nextv;
        if (s->svertes)
            s->svertes->prevv.reset();
    }
    else
    {
        prevv.lock()->nextv = nextv;
        if (nextv)
            nextv->prevv = prevv;
    }
    return true;
}