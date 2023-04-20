#include "face.h"
#include "solid.h"

face::face(Solid *s)
{
    floops = nullptr;
    flout = nullptr;

    nextf = s->sfaces;
    prevf = nullptr;
    if (s->sfaces)
        s->sfaces->prevf = this;
    s->sfaces = this;
    fsolid = s;
}

bool face::RemoveListFromSolid(Solid *s)
{
        // break;
    if (this == s->sfaces)
    {
        s->sfaces = s->sfaces->nextf;
        if (s->sfaces)
            s->sfaces->prevf = nullptr;
    }
    else
    {
        prevf->nextf = nextf;
        if (nextf)
            nextf->prevf = prevf;
    }
    return true;
}