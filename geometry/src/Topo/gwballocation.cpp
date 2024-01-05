#include "gwballocation.h"

void addlist(std::shared_ptr<Solid> s)
{
    if (s == firsts.lock())
    {
        firsts = s->nexts;
        if (s->nexts)
            s->nexts->prevs.reset();
    }
    else
    {
        s->prevs.lock()->nexts = s->nexts;
        if (s->nexts)
            s->nexts->prevs = s->prevs;
    }
}