#include "face.h"
#include "solid.h"
#include "surface.h"
#include "loop.h"
#include "params.h"

face::face()
{
    floops = nullptr;
    flout = nullptr;
    surf = nullptr;
    nextf = nullptr;
    prevf.reset();
    fsolid.reset();
}

face::~face() 
{
        std::shared_ptr<Face> fac = nextf;
        std::shared_ptr<Face> prefac = fac ? fac->prevf.lock() : nullptr;
        while (fac)
        {
            prefac = fac;
            fac = fac->nextf;
        }
        
        prefac = prefac ? prefac->prevf.lock() : nullptr;
        while (prefac)
        {
            prefac->nextf = nullptr;
            prefac = prefac->prevf.lock();
        } 

};

bool face::RemoveListFromSolid(std::shared_ptr<Solid> s)
{
        // break;
    if (shared_from_this() == s->sfaces)
    {
        s->sfaces = s->sfaces->nextf;
        if (s->sfaces)
            s->sfaces->prevf.reset();
    }
    else
    {
        prevf.lock()->nextf = nextf;
        if (nextf)
            nextf->prevf = prevf;
    }
    nextf = nullptr;
    prevf.reset();
    return true;
}

Eigen::Vector3d face::get_normal()
{
    return surf->get_normal();
}

void face::addlist(std::shared_ptr<Loop> l)
{
    l->nextl = floops;
    l->prevl.reset();
    if (floops)
        floops->prevl = l;
    floops = l;
    l->lface = shared_from_this();
}