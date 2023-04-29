#include "face.h"
#include "solid.h"
#include "surface.h"
#include "loop.h"

face::face()
{
    floops = nullptr;
    flout = nullptr;
    surf = nullptr;
    nextf = nullptr;
    prevf.reset();
}

// bool face::addtofacelist(std::shared_ptr<Solid> s)
// {
//     floops = nullptr;
//     flout = nullptr;

//     nextf = s->sfaces;
//     prevf.reset();
//     std::shared_ptr<Face> thisface = shared_from_this();
//     if (s->sfaces)
//         s->sfaces->prevf = thisface;
//     s->sfaces = thisface;
//     fsolid = s;
// }

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