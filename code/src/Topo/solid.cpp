#include "solid.h"
#include "params.h"
#include "face.h"
#include "vertex.h"

solid::solid()
{
    nexts = firsts;
    prevs = nullptr;
    if (firsts)
        firsts->prevs = this;
    firsts = this;
    this->sedges = nullptr;
    this->sfaces = nullptr;
    this->svertes = nullptr;
}

bool solid::RemoveListFromSolid()
{
    if (this == firsts)
    {
        firsts = firsts->nexts;
        if (firsts)
            firsts->prevs = nullptr;
    }
    else
    {
        prevs->nexts = nexts;
        if (nexts)
            nexts->prevs = prevs;
    }
    
    return true;
}

Face *solid::getface(Id faceno)
{
    Face *f = sfaces;
    while (f)
    {
        if (f->faceno == faceno)
            return f;
        f = f->nextf;
    }
    return nullptr;
}

Vertex *solid::getvertex(Id vertexno)
{
    Vertex *v = svertes;
    while (v)
    {
        if (v->vertexno == vertexno)
            return v;
        v = v->nextv;
    }
    return nullptr;
}