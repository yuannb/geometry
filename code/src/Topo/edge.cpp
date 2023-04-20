#include "edge.h"
#include "solid.h"

edge::edge(Solid *s)
{
    he1 = nullptr;
    he2 = nullptr;
    nexte = s->sedges;
    preve = nullptr;

    if (s->sedges)
        s->sedges->preve = this;
    s->sedges = this;
}
bool edge::RemoveListFromSolid(Solid *s)
{
    if (this == s->sedges)
    {
        s->sedges = s->sedges->nexte;
        if (s->sedges)
            s->sedges->preve = nullptr;
    }
    else
    {
        this->preve->nexte = this->nexte;
        if (this->nexte)
            this->nexte->preve = this->preve;
    }

    return true;
}
edge::~edge()
{
}