#include "edge.h"
#include "solid.h"
#include <fstream>
#include <algorithm>
#include <iosfwd>
#include "params.h"
edge::edge()
{
    he1.reset();
    he2.reset();
    nexte = nullptr;
    preve.reset();
    cur = nullptr;
}
bool edge::RemoveListFromSolid(std::shared_ptr<Solid> s)
{
    if (this == s->sedges.get())
    {
        s->sedges = s->sedges->nexte;
        if (s->sedges)
            s->sedges->preve.reset();
    }
    else
    {
        preve.lock()->nexte = nexte;
        if (this->nexte)
            this->nexte->preve = preve;
    }
    nexte = nullptr;
    preve.reset();
    return true;
}
edge::~edge()
{
    std::shared_ptr<Edge> edg = nexte;
    std::shared_ptr<Edge> pre = nexte ? nexte->preve.lock() : nullptr;
    while (edg)
    {
        pre = edg;
        edg = edg->nexte;
    }
    
    pre = pre ? pre->preve.lock() : nullptr;
    while (pre)
    {
        pre->nexte = nullptr;
        pre = pre->preve.lock();
    } 
}