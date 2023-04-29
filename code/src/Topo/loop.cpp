#include "loop.h"
#include "face.h"
loop::loop()
{
}


bool loop::RemoveListFromFace(std::shared_ptr<Face> f)
{
    if (this == f->flout.get())
        f->flout = nullptr;
    if (this == f->floops.get())
    {
        f->floops = this->nextl;
        if (f->floops)
            f->floops->prevl.reset();
    }
    else
    {
        prevl.lock()->nextl = nextl;
        if (nextl)
            nextl->prevl = prevl;
    }
    return true;
}