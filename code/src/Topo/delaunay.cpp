#include "delaunay.h"

triangel  delaunay(Solid *s)
{
    Vertex *sv = s->svertes;
    triangel tri;
    while (sv)
    {
        tri.vtxarry.push_back(sv);
        sv = sv->nextv; 
    }
    Face *sf = s->sfaces;
    while (sf)
    {
        Loop *l = sf->flout;
        HalfEdge *he = l->ledg;
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