#include "query.h"

void listsolid(std::shared_ptr<Solid> s)
{
    std::shared_ptr<Face> f = nullptr;
    std::shared_ptr<Loop> l = nullptr;
    std::shared_ptr<HalfEdge> he = nullptr;
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

void listneighbors(std::shared_ptr<Vertex> v)
{
    std::shared_ptr<HalfEdge> adj = v->vedge.lock();
    std::shared_ptr<HalfEdge> first = adj;
    if (adj)
    {
        do
        {
            std::cout << adj->nxt->vtx->vertexno << ", ";
        } while ((adj = adj->mate()->nxt) != first);
    }
    else
    {
        std::cout << "no adjacent vertices";
    }
    std::cout << std::endl;
}