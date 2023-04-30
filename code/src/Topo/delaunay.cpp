#include "delaunay.h"

triangel  delaunay(std::shared_ptr<Solid> s)
{
    std::shared_ptr<Vertex> sv = s->svertes;
    triangel tri;
    while (sv)
    {
        tri.vtxarry.push_back(sv);
        sv = sv->nextv; 
    }
    std::shared_ptr<Face> sf = s->sfaces;
    while (sf)
    {
        std::shared_ptr<Loop> l = sf->flout;
        std::shared_ptr<HalfEdge> he = l->ledg;
        std::vector<int> vec;
        do
        {
            auto it = std::find(tri.vtxarry.begin(), tri.vtxarry.end(), he->vtx);
            int num = it - tri.vtxarry.begin() + 1;
            vec.push_back(num);
        } while ((he = he->nxt) != l->ledg);
        tri.face.push_back(vec);
        sf = sf->nextf;
    }
    return tri;
}

triangel  discret(std::shared_ptr<Solid> s)
{
    std::shared_ptr<Vertex> v = s->svertes;
    int num = 1;
    triangel tri;
    while (v)
    {
        v->vertexno = num++;
        tri.vtxarry.push_back(v);
        v = v->nextv;
    }
    std::shared_ptr<Face> f = s->sfaces;
    while (f)
    {
        std::shared_ptr<Loop> l = f->floops;
        std::vector<int> facno;
        while (l)
        {
            std::shared_ptr<HalfEdge> he = l->ledg;
            do
            {
                facno.push_back(he->vtx->vertexno);

            }while((he = he->nxt) != l->ledg);
            l = l->nextl;
        }
        tri.face.push_back(facno);
        f = f->nextf;       
    }
    return tri;
}
