#include "vertex.h"
#include "solid.h"

vertex::vertex(Solid *s)
{
    nextv = s->svertes;
    prevv = nullptr;
    if (s->svertes)
        s->svertes->prevv = this;
    s->svertes = this;
    this->vedge = nullptr;
}

bool vertex::RemoveListFromSolid(Solid *s)
{
    if (this == s->svertes)
    {
        s->svertes = s->svertes->nextv;
        if (s->svertes)
            s->svertes->prevv = nullptr;
    }
    else
    {
        prevv->nextv = nextv;
        if (nextv)
            nextv->prevv = prevv;
    }
    return true;
}