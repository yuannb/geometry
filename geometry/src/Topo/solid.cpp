#include "solid.h"
#include "params.h"
#include "face.h"
#include "vertex.h"
#include "edge.h"

solid::solid()
{
    sfaces = nullptr;
    sedges = nullptr;
    svertes = nullptr;
    nexts = nullptr;
    prevs.reset();
}
solid::~solid()
{
    std::shared_ptr<Solid> snext = nexts;
    std::shared_ptr<Solid> presnext = snext ? snext->prevs.lock() : nullptr;
    while (snext)
    {
        presnext = snext;
        snext = snext->nexts;
    }
    presnext = presnext ? presnext->prevs.lock() : nullptr;
    while (presnext)
    {
        presnext->nexts = nullptr;
        presnext = presnext->prevs.lock();
    }
}

bool solid::addlist()
{
    nexts = firsts.lock();
    prevs.reset();
    if (nexts)
        nexts->prevs = shared_from_this();
    firsts = shared_from_this();
    this->sedges = nullptr;
    this->sfaces = nullptr;
    this->svertes = nullptr;
    return true;
};

bool solid::RemoveListFromSolid()
{
    std::shared_ptr<Solid> first = firsts.lock();
    if (shared_from_this() == first)
    {
        firsts = first->nexts;
        if (first->nexts)
            first->nexts->prevs.reset();
    }
    else
    {
        prevs.lock()->nexts = nexts;
        if (nexts)
            nexts->prevs = prevs;
    }
    prevs.reset();
    nexts = nullptr;
    
    return true;
}

std::shared_ptr<Face> solid::getface(Id faceno)
{
    std::shared_ptr<Face> f = sfaces;
    while (f)
    {
        if (f->faceno == faceno)
            return f;
        f = f->nextf;
    }
    return nullptr;
}

std::shared_ptr<Vertex> solid::getvertex(Id vertexno)
{
    std::shared_ptr<Vertex> v = svertes;
    while (v)
    {
        if (v->vertexno == vertexno)
            return v;
        v = v->nextv;
    }
    return nullptr;
}

void solid::addlist(std::shared_ptr<Edge> edg)
{
    edg->nexte = sedges;
    edg->preve.reset();
    if (sedges)
        sedges->preve = edg;
    sedges = edg;
}
void solid::addlist(std::shared_ptr<Face> fac)
{
    fac->nextf = sfaces;
    fac->prevf.reset();
    if (sfaces)
        sfaces->prevf = fac;
    sfaces = fac;
    fac->fsolid = shared_from_this();
}
void solid::addlist(std::shared_ptr<Vertex> vtx)
{
    vtx->nextv = svertes;
    vtx->prevv.reset();
    if (svertes)
        svertes->prevv = vtx;
    svertes = vtx;
}