#include "loop.h"
#include "face.h"
#include "params.h"
loop::loop()
{
    ledg = nullptr; 
    lface.reset(); 
    nextl = nullptr; 
    prevl.reset(); 
}

loop::~loop() 
{ 
    std::shared_ptr<Loop> cl = nextl;
    std::shared_ptr<Loop> prel = nextl ? nextl->prevl.lock() : nullptr;
    while (cl)
    {
        prel = cl;
        cl = cl->nextl;
    }
    
    prel = prel ? prel->prevl.lock() : nullptr;
    while (prel)
    {
        prel->nextl = nullptr;
        prel = prel->prevl.lock();
    } 
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
    nextl = nullptr;
    prevl.reset();
    return true;
}