#include "vertex.h"
#include "solid.h"
#include "params.h"
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
    nextv = nullptr;
    prevv.reset();
    return true;
}

vertex::~vertex() 
{ 
    std::shared_ptr<Vertex> snexv = nextv;
    std::shared_ptr<Vertex> presnext = nextv ? nextv->prevv.lock() : nullptr;
    while (snexv)
    {
        presnext = snexv;
        snexv = snexv->nextv;
    }
    presnext = presnext ? presnext->prevv.lock() : nullptr;
    while (presnext)
    {
        presnext->nextv = nullptr;
        presnext = presnext->prevv.lock();
    }
}