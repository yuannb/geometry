#include "loop.h"
#include "face.h"
loop::loop(Face *f)
{
    nextl = f->floops;
    prevl = nullptr;
    if (f->floops)
        f->floops->prevl = this;
    f->floops = this;
    lface = f;
}


bool loop::RemoveListFromFace(Face *f)
{
    if (this == f->flout)
        f->flout = nullptr;
    if (this == f->floops)
    {
        f->floops = this->nextl;
        if (f->floops)
            f->floops->prevl = nullptr;
    }
    else
    {
        prevl->nextl = nextl;
        if (nextl)
            nextl->prevl = prevl;
    }
    return true;
}