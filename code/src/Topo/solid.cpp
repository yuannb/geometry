#include "solid.h"
#include "params.h"

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