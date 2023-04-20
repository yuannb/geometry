#include "query.h"

void listsolid(Solid *s)
{
    Face *f;
    Loop *l;
    HalfEdge *he;
    f = s->sfaces;
    while (f)
    {
        std::cout << "face : " << f->faceno << " " << std::endl;
        l = f->floops;
        while (l)
        {
            std::cout << "loop : "<< std::endl;
            he = l->ledg;
            do
            {
                std::cout << he->vtx->vertexno << " :("
                << he->vtx->vcoord[0] << ", "
                << he->vtx->vcoord[1] << ", "
                << he->vtx->vcoord[2] << "); ";
            } while ((he = he->nxt) != l->ledg);
            std::cout << std::endl;
            l = l->nextl;
        }
        f = f->nextf;
    }   
}

void listneighbors(Vertex *v)
{
    HalfEdge    *adj;
    adj = v->vedge;
    if (adj)
    {
        do
        {
            std::cout << adj->nxt->vtx->vertexno << ", ";
        } while ((adj = mate(adj)->nxt) != v->vedge);
    }
    else
    {
        std::cout << "no adjacent vertices";
    }
    std::cout << std::endl;
}