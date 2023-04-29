#include "gwballocation.h"

// unsigned nodesize[] = 
// {
//     sizeof(Solid), sizeof(Face), sizeof(Loop),
//     sizeof(HalfEdge), sizeof(Edge), sizeof(Vertex), 0
// };

// void addlist(int what, Node *which, Node *where)
// {
//     switch (what)
//     {
//     case SOLID:
//         which->s.nexts = firsts;
//         which->s.prevs = (Solid*) NIL;
//         if (firsts)
//             firsts->prevs = (Solid *)which;
//         firsts = (Solid*)which;
//         break;
//     case FACE:
//         which->f.nextf = where->s.sfaces;
//         which->f.prevf = (Face*)NIL;
//         if (where->s.sfaces)
//             where->s.sfaces->prevf = (Face*)which;
//         where->s.sfaces = (Face*)which;
//         which->f.fsolid = (Solid*)where;
//         break;
//     case LOOP:
//         which->l.nextl = where->f.floops;
//         which->l.prevl = (Loop*)NIL;
//         if (where->f.floops)
//             where->f.floops->prevl = (Loop*)which;
//         where->f.floops = (Loop*)which;
//         which->l.lface = (Face*)where;
//         break;
//     case EDGE:
//         which->e.nexte = where->s.sedges;
//         which->e.preve = (Edge*)NIL;
//         if (where->s.sedges)
//             where->s.sedges->preve = (Edge*)which;
//         where->s.sedges = (Edge*)which;
//         break;
//     case VERTEX:
//         which->v.nextv = where->s.svertes;
//         which->v.prevv = (Vertex*)NIL;
//         if (where->s.svertes)
//             where->s.svertes->prevv = (Vertex*)which;
//         where->s.svertes = (Vertex*)which;
//         break;
//     default:
//         break;
//     }
// }

// void removelist(int what, Node *which, Node *where)
// {
//     switch (what)
//     {
//     case SOLID:
//         if (&which->s == firsts)
//             firsts = firsts->nexts;
//             if (firsts)
//                 firsts->prevs = (Solid*)NIL;
//         else
//         {
//             which->s.prevs->nexts = which->s.nexts;
//             if (which->s.nexts)
//                 which->s.nexts->prevs = which->s.prevs;
//         }
//         break;
//     case FACE:
//         if (&which->f == where->s.sfaces)
//         {
//             where->s.sfaces = where->s.sfaces->nextf;
//             if (where->s.sfaces)
//                 where->s.sfaces->prevf = (Face*)NIL;
//         }
//         else
//         {
//             which->f.prevf->nextf = which->f.nextf;
//             if (which->f.nextf)
//                 which->f.nextf->prevf = which->f.prevf;
//         }
//         break;
//     case EDGE:
//         if (&which->e == where->s.sedges)
//         {
//             where->s.sedges = where->s.sedges->nexte;
//             if (where->s.sedges)
//                 where->s.sedges->preve = (Edge*)NIL;
//         }
//         else
//         {
//             which->e.preve->nexte = which->e.nexte;
//             if (which->e.nexte)
//                 which->e.nexte->preve = which->e.preve;
//         }
//         break;
//     case LOOP:
//         if (&which->l == where->f.flout)
//             where->f.flout = (Loop*)NIL;
//         if (&which->l == where->f.floops)
//         {
//             where->f.floops = where->f.floops->nextl;
//             if (where->f.floops)
//                 where->f.floops->prevl = nullptr;
//         }
//         else
//         {
//             which->l.prevl->nextl = which->l.nextl;
//             if (which->l.nextl)
//                 which->l.nextl->prevl = which->l.prevl;
//         }
//         break;
//     case VERTEX:
//         if (&which->v == where->s.svertes)
//         {
//             where->s.svertes = where->s.svertes->nextv;
//             if (where->s.svertes)
//                 where->s.svertes->prevv == (Vertex*)NIL;
//         }
//         else
//         {
//             which->v.prevv->nextv = which->v.nextv;
//             if (which->v.nextv)
//                 which->v.nextv->prevv = which->v.prevv;
//         }
//         break;
//     default:
//         break;
//     }
// }

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