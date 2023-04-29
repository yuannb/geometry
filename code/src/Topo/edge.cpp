#include "edge.h"
#include "solid.h"
#include <fstream>
#include <algorithm>
#include <iosfwd>
edge::edge()
{
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

    return true;
}
edge::~edge()
{
}