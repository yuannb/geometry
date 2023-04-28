#include "gwballocation.h"

unsigned nodesize[] = 
{
    sizeof(Solid), sizeof(Face), sizeof(Loop),
    sizeof(HalfEdge), sizeof(Edge), sizeof(Vertex), 0
};

void addlist(int what, Node *which, Node *where)
{
    switch (what)
    {
    case SOLID:
        which->s.nexts = firsts;
        which->s.prevs = (Solid*) NIL;
        if (firsts)
            firsts->prevs = (Solid *)which;
        firsts = (Solid*)which;
        break;
    case FACE:
        which->f.nextf = where->s.sfaces;
        which->f.prevf = (Face*)NIL;
        if (where->s.sfaces)
            where->s.sfaces->prevf = (Face*)which;
        where->s.sfaces = (Face*)which;
        which->f.fsolid = (Solid*)where;
        break;
    case LOOP:
        which->l.nextl = where->f.floops;
        which->l.prevl = (Loop*)NIL;
        if (where->f.floops)
            where->f.floops->prevl = (Loop*)which;
        where->f.floops = (Loop*)which;
        which->l.lface = (Face*)where;
        break;
    case EDGE:
        which->e.nexte = where->s.sedges;
        which->e.preve = (Edge*)NIL;
        if (where->s.sedges)
            where->s.sedges->preve = (Edge*)which;
        where->s.sedges = (Edge*)which;
        break;
    case VERTEX:
        which->v.nextv = where->s.svertes;
        which->v.prevv = (Vertex*)NIL;
        if (where->s.svertes)
            where->s.svertes->prevv = (Vertex*)which;
        where->s.svertes = (Vertex*)which;
        break;
    default:
        break;
    }
}

// Node *gnew(int what, Node *where)
// {
//     Node *node;
//     node = (Node*) std::malloc(nodesize[what]);
//     switch (what)
//     {
//     case SOLID:
//         addlist(SOLID, node, NIL);
//         node->s.sfaces = (Face*)NIL;
//         node->s.sedges = (Edge*)NIL;
//         node->s.svertes = (Vertex*)NIL;
//         break;
//     case FACE:
//         addlist(FACE, node, where);
//         node->f.floops = (Loop*)NIL;
//         node->f.flout = (Loop*)NIL;
//         break;
//     case LOOP:
//         addlist(LOOP, node, where);
//         break;
//     case EDGE:
//         addlist(EDGE, node, where);
//         break;;
//     case VERTEX:
//         addlist(VERTEX, node, where);
//         // node->v.vedge = (Vertex*)NIL; //this should be a bug
//         break;
//     default:
//         break;
//     }
//     return node;
// }

void removelist(int what, Node *which, Node *where)
{
    switch (what)
    {
    case SOLID:
        if (&which->s == firsts)
            firsts = firsts->nexts;
            if (firsts)
                firsts->prevs = (Solid*)NIL;
        else
        {
            which->s.prevs->nexts = which->s.nexts;
            if (which->s.nexts)
                which->s.nexts->prevs = which->s.prevs;
        }
        break;
    case FACE:
        if (&which->f == where->s.sfaces)
        {
            where->s.sfaces = where->s.sfaces->nextf;
            if (where->s.sfaces)
                where->s.sfaces->prevf = (Face*)NIL;
        }
        else
        {
            which->f.prevf->nextf = which->f.nextf;
            if (which->f.nextf)
                which->f.nextf->prevf = which->f.prevf;
        }
        break;
    case EDGE:
        if (&which->e == where->s.sedges)
        {
            where->s.sedges = where->s.sedges->nexte;
            if (where->s.sedges)
                where->s.sedges->preve = (Edge*)NIL;
        }
        else
        {
            which->e.preve->nexte = which->e.nexte;
            if (which->e.nexte)
                which->e.nexte->preve = which->e.preve;
        }
        break;
    case LOOP:
        if (&which->l == where->f.flout)
            where->f.flout = (Loop*)NIL;
        if (&which->l == where->f.floops)
        {
            where->f.floops = where->f.floops->nextl;
            if (where->f.floops)
                where->f.floops->prevl = nullptr;
        }
        else
        {
            which->l.prevl->nextl = which->l.nextl;
            if (which->l.nextl)
                which->l.nextl->prevl = which->l.prevl;
        }
        break;
    case VERTEX:
        if (&which->v == where->s.svertes)
        {
            where->s.svertes = where->s.svertes->nextv;
            if (where->s.svertes)
                where->s.svertes->prevv == (Vertex*)NIL;
        }
        else
        {
            which->v.prevv->nextv = which->v.nextv;
            if (which->v.nextv)
                which->v.nextv->prevv = which->v.prevv;
        }
        break;
    // case HALFEDGE:
    //     if (&which->h == which->h.wloop->ledg)
    //     {
    //         if (which->h.nxt == &which->h)
    //             which->h.wloop->ledg = (HalfEdge*)NIL;
    //         else
    //         {
    //             // which->h.nxt = which->h.prv;
    //             // which->h.prv = which->h.nxt;
    //             which->h.wloop->ledg = which->h.prv;
    //         }
    //     }
    //     // else
    //     // {
    //     //         which->h.nxt = which->h.prv;
    //     //         which->h.prv = which->h.nxt;
    //     // }
    //     break;
    default:
        break;
    }
}

// void Free(Solid *dsolid)
// {
//     // free vertex
//     Vertex *vtx = dsolid->svertes;
//     while (vtx)
//     {
//         Vertex *nxtvex = vtx->nextv;
//         free(vtx);
//         vtx = nxtvex;
//     }

//     //free face
//     Face *face = dsolid->sfaces;
//     while (face)
//     {
//         Face *nxtface = face->nextf;
//         Loop *dl = face->floops;
//         while (dl)
//         {
//             Loop *nxtdl = dl->nextl;
//             HalfEdge *he = dl->ledg;
//             if (he == NIL)
//                 continue;
//             if (he->nxt = he)
//             {
//                 free(he);
//                 continue;
//             }
//             else
//             {
//                 he->prv->nxt = NIL;
//             }
//             while (he)
//             {
//                 HalfEdge *nxthe = he->nxt;
//                 free(he);
//                 he = nxthe;
//             }
            
//             free(dl);
//         }
//         free(face);
//         face = nxtface;
//     }
    
//     //free edge
//     Edge *edg = dsolid->sedges;
//     while (edg)
//     {
//         Edge *nxtedg = edg->nexte;
//         free(edg);
//         edg = nxtedg;
//     }
// }

// void Free(Face *face)
// {
//     Loop *dl = face->floops;
//     while (dl)
//     {
//         Loop *nxtdl = dl->nextl;
//         HalfEdge *he = dl->ledg;
//         if (he == NIL)
//             continue;
//         if (he->nxt = he)
//         {
//             free(he);
//             continue;
//         }
//         else
//         {
//             he->prv->nxt = NIL;
//         }
//         while (he)
//         {
//             HalfEdge *nxthe = he->nxt;
//             free(he);
//             he = nxthe;
//         }
        
//         free(dl);
//     }
//     free(face);
// }

// void Free(Loop *dl)
// {
//     HalfEdge *he = dl->ledg;
//     if (he == NIL)
//     {
//         free(dl);
//         return;
//     }

//     if (he->nxt = he)
//     {
//         free(he);
//     }
//     else
//     {
//         he->prv->nxt = NIL;
//     }
//     while (he)
//     {
//         HalfEdge *nxthe = he->nxt;
//         free(he);
//         he = nxthe;
//     }
    
//     free(dl);
// }

// void Free(Edge *edg)
// {
//     if (edg)
//         free(edg);
// }

// void Free(Vertex *vtx)
// {
//     if (vtx)
//         free(vtx);
// }

// void Free(HalfEdge *he)
// {
//     if (he)
//         free(he);
// }

// HalfEdge *findOtherHalfEdge(HalfEdge *he)
// {
//     HalfEdge *otherhe = mate(he)->nxt;
//     if (otherhe == NIL || otherhe == he)
//         return NIL;
//     return otherhe;
// }



// // TODO : free storage
// void gdel(int what, Node *which, Node *where)
// {
//     // Node *node;
//     switch (what)
//     {
//     case SOLID:
//         removelist(what, which, where);
//         Free(&which->s);
//         break;
//     case FACE:
//         removelist(what, which, where);
//         // Loop *l = which->f.floops;
//         // while (l)
//         // {
//         //     HalfEdge *he = l->ledg;
//         //     if (he)
//         //     {
//         //         HalfEdge *otherhe = findOtherHalfEdge(he);
//         //         he->vtx->vedge = otherhe;
//         //         do
//         //         {
//         //             if (he->edg->he1 == he)
//         //                 he->edg->he1 = (HalfEdge*)NIL;
//         //             else
//         //                 he->edg->he2 = (HalfEdge*)NIL;
//         //         } while ((he = he->nxt) != l->ledg);
//         //     }
//         //     l = l->nextl;
//         // }
//         Free(&which->f);
//         break;
//     case LOOP:
//         removelist(what, which, where);
//         // HalfEdge *he = which->l.ledg;
//         // if (he)
//         // {
//         //     HalfEdge *otherhe = findOtherHalfEdge(he);
//         //     he->vtx->vedge = otherhe;
//         //     do
//         //     {
//         //         if (he->edg->he1 == he)
//         //             he->edg->he1 = (HalfEdge*)NIL;
//         //         else
//         //             he->edg->he2 = (HalfEdge*)NIL;
//         //     } while ((he = he->nxt) != which->l.ledg);
//         // }
//         Free(&which->l);
//         break;
//     case EDGE:
//         removelist(what, which, where);
//         // which->e.he1->edg = (Edge*)NIL;
//         // which->e.he2->edg = (Edge*)NIL;
//         Free(&which->e);
//         break;
//     case VERTEX:
//         removelist(what, which, where);
//         // HalfEdge *he =  which->v.vedge;
//         // if (he)
//         // {
//         //     do
//         //     {
//         //         he->vtx = NIL;
//         //     } while ((he = mate(he)->nxt) != which->v.vedge);
            
//         // }
//         Free(&which->v);
//         break;
//     case HALFEDGE:
//         removelist(HALFEDGE, which, NIL);
//         // HalfEdge *he = findOtherHalfEdge(&which->h);
//         // which->h.vtx->vedge = he;
//         Free(&which->h);
//         break;
//     }
// }