#include "euler.h"
// #include <cstdlib>
Id maxv = 0;
Id maxf = 0;

HalfEdge *addhe(Edge *e, Vertex *v, HalfEdge *where, int sign)
{
    HalfEdge *he;

    if (where->edg == NIL)
    {
        he = where;
    }
    else
    {
        // he = (HalfEdge*) gnew(HALFEDGE, NIL);
        he = new HalfEdge();
        where->prv->nxt = he;
        he->prv = where->prv;
        where->prv = he;
        he->nxt = where;
    }
    he->edg = e;
    he->vtx = v;
    he->wloop = where->wloop;
    if (sign == PLUS)
        e->he1 = he;
    else e->he2 = he;
    return he;
}

HalfEdge *delhe(HalfEdge *he)
{
    if (he->edg == NIL)
    {
        // gdel(HALFEDGE, (Node*)he, NIL);
        delete he;
        return NIL;
    }
    else if (he->nxt == he)
    {
        he->edg = NIL;
        return he;
    }
    else
    {
        he->prv->nxt = he->nxt;
        he->nxt->prv = he->prv;
        HalfEdge *hepre = he->prv;
        // gdel(HALFEDGE, (Node*)he, NIL);
        delete he;
        return hepre;
    }
}

Solid *getsolid(Id sn)
{
    Solid *s;
    for (s = firsts; s != NIL; s = s->nexts)
    {
        if (s->solidno == sn)
            return s;
    }
    return NIL;
}

Face    *fface(Solid *s, Id fn)
{
    for (Face *f = s->sfaces; f != NIL; f = f->nextf)
    {
        if (f->faceno == fn)
            return f;
    }
    return NIL;
}

HalfEdge *fhe(Face *f, Id vn1, Id vn2)
{
    for (Loop *l = f->floops; l != NIL; l = l->nextl)
    {
        HalfEdge *h = l->ledg;
        do
        {
            if (h->vtx->vertexno == vn1 &&
                h->nxt->vtx->vertexno == vn2)
                    return h;
        } while ((h = h->nxt) != l->ledg);
    }
    return NIL;
}

HalfEdge *fhe(Face *f, Id vn1)
{
    for (Loop *l = f->floops; l != NIL; l = l->nextl)
    {
        HalfEdge *h = l->ledg;
        do
        {
            if (h->vtx->vertexno == vn1)
                    return h;
        } while ((h = h->nxt) != l->ledg);
    }
    return NIL;
}

Vertex *getvertex(Solid *s, Id vno)
{
    Vertex *v = s->svertes;
    while (v)
    {
        if (v->vertexno == vno)
            return v;
        v = v->nextv;
    }
    return NIL;
}

Solid *mvfs(Id s, Id f, Id v, double x, double y, double z)
{
    Solid *newsolid;
    Face *newface;
    Vertex *newvertex;
    HalfEdge *newhe;
    Loop *newloop;
    // newsolid = (Solid*)gnew(SOLID, NIL);
    newsolid = new Solid();
    // newface = (Face*)gnew(FACE, (Node*)newsolid);
    newface = new face(newsolid);
    newloop = new Loop(newface);
    // newloop = (Loop*)gnew(LOOP, (Node*)newface);
    newvertex = new Vertex(newsolid);
    // newvertex = (Vertex*)gnew(VERTEX, (Node*)newsolid);
    // newhe = (HalfEdge*)gnew(HALFEDGE, NIL);
    newhe = new HalfEdge();

    newsolid->solidno = s;
    newface->faceno = f;
    newface->flout = newloop;
    newloop->ledg = newhe;
    newhe->wloop = newloop;
    newhe->nxt = newhe;
    newhe->prv = newhe;
    newhe->vtx = newvertex;
    newhe->edg = NIL;
    newvertex->vertexno = v;
    newvertex->vcoord[0] = x;
    newvertex->vcoord[1] = y;
    newvertex->vcoord[2] = z;
    return newsolid;
}

void lmev(HalfEdge *he1, HalfEdge *he2, Id v, double x, double y, double z)
{
    HalfEdge *he;
    Vertex  *newvertex;
    Edge    *newedge;

    // newedge = (Edge*)gnew(EDGE, (Node*)he1->wloop->lface->fsolid);
    newedge = new Edge(he1->wloop->lface->fsolid);
    // newvertex = (Vertex*)gnew(VERTEX, (Node*)he1->wloop->lface->fsolid);
    newvertex = new Vertex(he1->wloop->lface->fsolid);
    newvertex->vcoord[0] = x;
    newvertex->vcoord[1] = y;
    newvertex->vcoord[2] = z;
    newvertex->vertexno = v;

    he = he1;
    while (he != he2)
    {
        he->vtx = newvertex;
        he = mate(he)->nxt;
    }
    
    // the code of book may be wrong (P184)
    addhe(newedge, he2->vtx, he1, MINUS);
    addhe(newedge, newvertex, he2, PLUS);
    newvertex->vedge = he2->prv;
    he2->vtx->vedge = he2;
}

Face *lmef(HalfEdge *he1, HalfEdge *he2, Id f)
{
    Face *newface;
    Loop *newloop;
    HalfEdge *he, *nhe1, *nhe2, *temp;
    Edge    *newedge;
    newedge = new Edge(he1->wloop->lface->fsolid);
    newface = new Face(he1->wloop->lface->fsolid);
    newloop = new Loop(newface);
    newface->faceno = f;
    newface->flout = newloop;

    he = he1;
    while (he != he2)
    {
        he->wloop = newloop;
        he = he->nxt;
    }
    nhe1 = addhe(newedge, he2->vtx, he1, MINUS);
    nhe2 = addhe(newedge, he1->vtx, he2, PLUS);

    nhe1->prv->nxt = nhe2;
    nhe2->prv->nxt = nhe1;
    temp = nhe1->prv;
    nhe1->prv = nhe2->prv;
    nhe2->prv = temp;

    newloop->ledg = nhe1;
    he2->wloop->ledg = nhe2;
    return newface;
}

void lkef(HalfEdge *he1, HalfEdge *he2)
{
    Edge *edg = he1->edg;
    if (edg != he2->edg)
    {
        std::cerr << "lkef : he1->edg != he2.edg" << std::endl;
        return;
    }
    if (he1->wloop == he2->wloop)
    {
        std::cerr << "lkef : he1->wloop == he2->wloop" << std::endl;      
        return;
    }
    if (he1->wloop->lface == he2->wloop->lface)
    {
        std::cerr << "lkef : he1->wloop->lface == he2->wloop->lface" << std::endl;
        return;
    }

    HalfEdge *nxthe1 = he1->nxt;
    halfedge *nxthe2 = he2->nxt;

    //TODO: loop only has one halfedge

    //we assume loop has two halfedge at least
    HalfEdge *he = he2->nxt;
    he1->wloop->ledg = nxthe1;
    do
    {
        he->wloop = he1->wloop;
    } while ((he = he->nxt) != he2);
    
    he2->vtx->vedge = he1->nxt; 
    he1->vtx->vedge = he2->nxt;
    
    if (he2 == he2->vtx->vedge)
    {
        he2->vtx->vedge = he1->nxt;
    }

    he2->prv->nxt = he1->nxt;
    he1->nxt->prv = he2->prv;
    he1->nxt = he2;
    he2->prv = he1;

    Solid *s = he1->wloop->lface->fsolid;
    Loop *dl = he2->wloop;
    Face *df = he2->wloop->lface;
    delhe(he1);
    delhe(he2);
    edg->RemoveListFromSolid(s);
    delete edg;
    dl->RemoveListFromFace(df);
    delete dl;

    // TODO : consider inner loop
    // gdel(FACE, (Node*)df, (Node*)s);
    df->RemoveListFromSolid(s);
    delete df;
}

int kef(Id s, Id f, Id v1, Id v2)
{
    Solid *oldsolid;
    Face *oldface;
    Edge *e;
    HalfEdge *h1, *h2;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "kef : solid " << s << " not found" << std::endl;
        return ERROR; 
    }
    if ((oldface = fface(oldsolid, f)) == NIL)
    {
        std::cerr << "kef : face " << f << " not found in solid " << s << std::endl;
        return ERROR; 
    }
    Loop *l = oldface->floops;
    bool flag = false;
    while(l)
    {
        h1 = oldface->floops->ledg;
        do
        {
            if (h1->vtx->vertexno == v1 && h1->nxt->vtx->vertexno == v2)
            {
                h2 = mate(h1);
                flag = true;
                break;
            }
            else
            {
                if (h1->vtx->vertexno == v2 && h1->nxt->vtx->vertexno == v1)
                {
                   h2 = mate(h1);
                   flag = true;
                   break; 
                }
            }
            if (flag == true)
                break;
        } while ((h1 = h1->nxt) != oldface->floops->ledg);
        l = l->nextl;
    }
    if (flag == false)
    {
        std::cerr << "kef : edge " << v1 << "-" << v2 << " not found in face " << f << std::endl;
        return ERROR; 
    }
    lkef(h1, h2);
    return SUCCESS;
}

void lkemr(HalfEdge *h1, HalfEdge *h2)
{
    HalfEdge *h3, *h4;
    Loop *nl;
    Loop *ol;
    Edge *killedge;

    ol = h1->wloop;
    nl = new Loop(ol->lface);
    killedge = h1->edg;

    h3 = h1->nxt;
    h1->nxt = h2->nxt;
    h2->nxt->prv = h1;
    h3->prv = h2;
    h4 = h2;
    do
    {
        h4->wloop = nl;
    } while ((h4 = h4->nxt) != h2);
    ol->ledg = h3 = delhe(h1);
    nl->ledg = h4 = delhe(h2);

    // this may have a bug(P187)
    h3 = h3->nxt;
    h4 = h4->nxt;

    h3->vtx->vedge = (h3->edg) ? h3 : (HalfEdge*)NIL;
    h4->vtx->vedge = (h4->edg) ? h4 : (HalfEdge*)NIL;

    killedge->RemoveListFromSolid(ol->lface->fsolid);
    delete killedge;
}

int kemr(Id s, Id f, Id v1, Id v2)
{
    Solid *oldsolid;
    Face *oldface;
    Edge *e;
    HalfEdge *h1, *h2;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "kemr : solid " << s << " not found" << std::endl;
        return ERROR; 
    }
    if ((oldface = fface(oldsolid, f)) == NIL)
    {
        std::cerr << "kemr : face " << f << " not found in solid " << s << std::endl;
        return ERROR; 
    }
    Loop *l = oldface->floops;
    bool flag = false;
    while(l)
    {
        h1 = oldface->floops->ledg;
        do
        {
            if (h1->vtx->vertexno == v1 && h1->nxt->vtx->vertexno == v2)
            {
                h2 = mate(h1);
                flag = true;
                break;
            }
            else
            {
                if (h1->vtx->vertexno == v2 && h1->nxt->vtx->vertexno == v1)
                {
                   h2 = mate(h1);
                   flag = true;
                   break; 
                }
            }
            if (flag == true)
                break;
        } while ((h1 = h1->nxt) != oldface->floops->ledg);
        l = l->nextl;
    }
    if (flag == false)
    {
        std::cerr << "kemr : edge " << v1 << "-" << v2 << " not found in face " << f << std::endl;
        return ERROR; 
    }
    lkemr(h1, h2);
    return SUCCESS;
}

void lmekr(HalfEdge *h1, HalfEdge *h2)
{
    if (h1->wloop == h2->wloop)
    {
        std::cerr << "lmekr : h1->wloop == h2->wloop" << std::endl;
        return;
    }
    if (h1->wloop->lface != h2->wloop->lface)
    {
        std::cerr << "lmekr : h1->wloop->lface != h2->wloop->lface" << std::endl;
        return;
    }

    HalfEdge *he = h2;
    Loop *dl = h2->wloop;
    do
    {
        he->wloop = h1->wloop;
    } while ((he = he->nxt) != h2);
    
    Solid *s = h1->wloop->lface->fsolid;
    Edge *newedge = new Edge(s);

    // how to ensure h1->wloop is out loop? there code may have bug
    if (h1->wloop->lface->flout == NIL)
        h1->wloop->lface->flout = h1->wloop;
    
    HalfEdge *temp = h1->prv;
    HalfEdge *preh1 = addhe(newedge, h2->vtx, h1, PLUS);
    HalfEdge *preh2 = addhe(newedge, h1->vtx, h2, MINUS);

    preh1->prv = preh2->prv;
    preh2->prv->nxt = preh1;
    preh2->prv = temp;
    temp->nxt = preh2;

    dl->RemoveListFromFace(h1->wloop->lface);
    delete dl;
}

int mekr(Id s, Id f, Id v1, Id v2, Id v3, Id v4)
{
    Solid *oldsolid;
    Face *oldface1;
    HalfEdge *he1, *he2;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid, f)) == NIL)
    {
        std::cerr << "mev : face " << f << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he1 = fhe(oldface1, v1, v2)) == NIL)
    {
        std::cerr << "mev : edge " << v1 << "-" << v2 << " not found in face " << f << std::endl;
        return ERROR;          
    }
    if ((he2 = fhe(oldface1, v3, v4)) == NIL)
    {
        std::cerr << "mev : edge " << v3 << "-" << v4 << " not found in face " << f << std::endl;
        return ERROR;          
    }
    lmekr(he1, he2);
    return SUCCESS;
}

int smekr(Id s, Id f, Id v1, Id v3)
{
    Solid *oldsolid;
    Face *oldface1;
    HalfEdge *he1, *he2;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid, f)) == NIL)
    {
        std::cerr << "mev : face " << f << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he1 = fhe(oldface1, v1)) == NIL)
    {
        std::cerr << "mev : edge " << v1 << "-" << " not found in face " << f << std::endl;
        return ERROR;          
    }
    if ((he2 = fhe(oldface1, v3)) == NIL)
    {
        std::cerr << "mev : edge " << v3 << "-" << " not found in face " << f << std::endl;
        return ERROR;          
    }
    lmekr(he1, he2);
    return SUCCESS;
}

void lkfmrh(Face *fac1, Face *fac2)
{
    // assume fac2 has just one out loop
    Loop *l2 = fac2->floops;
    addlist(LOOP, (Node*)l2, (Node*)fac1);
    
    fac2->RemoveListFromSolid(fac2->fsolid);
    delete fac2;
}

int kfmrh(Id s, Id f1, Id f2)
{
    Solid *oldsolid;
    Face *oldface1, *oldface2;
    HalfEdge *he1, *he2;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid, f1)) == NIL)
    {
        std::cerr << "mev : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((oldface2 = fface(oldsolid, f2)) == NIL)
    {
        std::cerr << "mev : face " << f2 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    lkfmrh(oldface1, oldface2);
    return ERROR;
}

void lmfkrh(Loop *l, Id f)
{
    // assume l is a inner loop
    // Face *newface = (Face*)gnew(FACE, (Node*)l->lface->fsolid);
    Face *newface = new Face(l->lface->fsolid);
    newface->faceno = f;
    newface->floops = l;
    newface->flout = l;
    if (l == l->lface->floops)
    {
        l->lface->floops = l->nextl;
        l->nextl->prevl = NIL;
        l->nextl = NIL;
    }
    else
    {
        l->prevl->nextl = l->nextl;
        l->nextl->prevl = l->prevl;
        l->prevl = NIL;
        l->nextl = NIL;
    }
    l->lface = newface;
}

int mfkrh(Id s, Id f1, Id v1, Id v2, Id f2)
{
    Solid *oldsolid;
    Face *oldface1;
    HalfEdge *he;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mfkrh : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid, f1)) == NIL)
    {
        std::cerr << "mfkrh : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he = fhe(oldface1, v1, v2)) == NIL)
    {
        std::cerr << "mfkrh : " << v1 << "-" << v2 << " not found in face " << f1 << std::endl;
        return ERROR;
    }
    lmfkrh(he->wloop, f2);
    return ERROR;
}

void lringmv(Loop *l, Face *toface, int inout)
{
    if (inout != 0)
        toface->flout = l;
    
    // ensure l not in toface
    Loop *fl = toface->floops;
    while (fl)
    {
        if (fl == l)
            return;
        fl = fl->nextl;
    }

    toface->floops->prevl = l;
    l->nextl = toface->floops;
    l->lface = toface;
}

int ringmv(Solid *s, Id f1, Id f2, Id v1, Id v2, int inout)
{
    Face *oldface1, *oldface2;
    HalfEdge *he;
    if ((oldface1 = fface(s, f1)) == NIL)
    {
        std::cerr << "ringmv : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((oldface2 = fface(s, f2)) == NIL)
    {
        std::cerr << "ringmv : face " << f2 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he = fhe(oldface1, v1, v2)) == NIL)
    {
        std::cerr << "mev : " << v1 << "-" << v2 << " not found in face " << f1 << std::endl;
        return ERROR;
    }
    lringmv(he->wloop, oldface2, inout);
    return SUCCESS;
}

int mev(Id s, Id f1, Id f2, Id v1, Id v2, Id v3, Id v4, double x, double y, double z)
{
    Solid *oldsolid;
    Face *oldface1, *oldface2;
    HalfEdge *he1, *he2;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid, f1)) == NIL)
    {
        std::cerr << "mev : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((oldface2 = fface(oldsolid, f2)) == NIL)
    {
        std::cerr << "mev : face " << f2 << " not found in solid " << s << std::endl;
        return ERROR;       
    }
    if ((he1 = fhe(oldface1, v1, v2)) == NIL)
    {
        std::cerr << "mev : edge " << v1 << "-" << v2 << " not found in face " << f1 << std::endl;
        return ERROR;          
    }
    if ((he2 = fhe(oldface2, v1, v3)) == NIL)
    {
        std::cerr << "mev : edge " << v1 << "-" << v3 << " not found in face " << f2 << std::endl;
        return ERROR;          
    }
    lmev(he1, he2, v4, x, y, z);
    return SUCCESS;
}

int smev(Id s, Id f1, Id v1, Id v4, double x, double y, double z)
{
    Solid *oldsolid;
    Face *oldface1, *oldface2;
    HalfEdge *he1, *he2;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid, f1)) == NIL)
    {
        std::cerr << "mev : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he1 = fhe(oldface1, v1)) == NIL)
    {
        std::cerr << "mev : edge " << v1 << "-" << "v2" << " not found in face " << f1 << std::endl;
        return ERROR;          
    }
    lmev(he1, he1, v4, x, y, z);
    return SUCCESS;
}

int smef(Id s, Id f1, Id v1, Id v3, Id f2)
{
    Solid *oldsolid;
    Face *oldface1, *oldface2;
    HalfEdge *he1, *he2;
    oldsolid = getsolid(s);
    oldface1 = fface(oldsolid, f1);
    Loop *lp = oldface1->floops;
    while (lp)
    {
        HalfEdge *he = lp->ledg;
        do
        {
            if (he->vtx->vertexno == v1)
                he1 = he;
            if (he->vtx->vertexno == v3)
                he2 = he;
        } while ((he = he->nxt) != (lp->ledg));
        lp = lp->nextl;
    }
    lmef(he1, he2, f2);
    return SUCCESS;
}

void lkvfs(Solid *s)
{
    // gdel(SOLID, (Node*)s, NIL);
    s->RemoveListFromSolid();
    delete s;
}

void kvfs(Id s)
{
    Solid *oldsolid;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return;
    }
    lkvfs(oldsolid);
}

void lkev(HalfEdge *he1, HalfEdge *he2)
{
    if (he1->vtx == he2->vtx)
    {
        std::cerr << "he1->vtx == he2->vtx" << std::endl;
        return;
    }
    Solid *s = he1->wloop->lface->fsolid;
    Edge *de = he1->edg;
    Vertex *dv = he1->vtx;
    HalfEdge *he = he1;
    do
    {
        he->vtx = he2->vtx;
    } while ((he = mate(he)->nxt) != he1);
    he2->vtx->vedge = he2->nxt;
    delhe(he1);
    delhe(he2);

    de->RemoveListFromSolid(s);
    delete de;
    dv->RemoveListFromSolid(s);
    delete dv;
}

int kev(Id s, Id f, Id v1, Id v2)
{   
    Solid *oldsolid;
    Face *oldface;
    Edge *e;
    HalfEdge *h1, *h2;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "kev : solid " << s << " not found" << std::endl;
        return ERROR; 
    }
    if ((oldface = fface(oldsolid, f)) == NIL)
    {
        std::cerr << "kev : face " << f << " not found in solid " << s << std::endl;
        return ERROR; 
    }
    Loop *l = oldface->floops;
    bool flag = false;
    while(l)
    {
        h1 = oldface->floops->ledg;
        do
        {
            if (h1->vtx->vertexno == v1 && h1->nxt->vtx->vertexno == v2)
            {
                h2 = mate(h1);
                flag = true;
                break;
            }
            else
            {
                if (h1->vtx->vertexno == v2 && h1->nxt->vtx->vertexno == v1)
                {
                   h2 = mate(h1);
                   flag = true;
                   break; 
                }
            }
            if (flag == true)
                break;
        } while ((h1 = h1->nxt) != oldface->floops->ledg);
        l = l->nextl;
    }
    if (flag == false)
    {
        std::cerr << "kev : edge " << v1 << "-" << v2 << " not found in face " << f << std::endl;
        return ERROR; 
    }
    lkev(h1, h2);
    return SUCCESS;
}

void getmaxnames(Id sn)
{
    Solid *s;
    Vertex *v;
    Face *f;

    s = getsolid(sn);
    for (v = s->svertes, maxv = 0; v != NIL; v = v->nextv)
        if (v->vertexno > maxv) maxv = v->vertexno;
    for (f = s->sfaces, maxf = 0; f != NIL; f = f->nextf)
        if (f->faceno > maxf) maxf = f->faceno;
}

void merge(Solid *s1, Solid *s2)
{
    while (s2->sedges)
    {
        Edge *e = s2->sedges;
        removelist(EDGE, (Node*)s2->sedges, (Node*)s2);
        addlist(EDGE, (Node*)e, (Node*)s1);
    }
    while (s2->sfaces)
    {
        Face *f = s2->sfaces;
        removelist(FACE, (Node*)s2->sfaces, (Node*)s2);
        addlist(FACE, (Node*)f, (Node*)s1);
    }
    while (s2->svertes)
    {
        Vertex *v = s2->svertes;
        removelist(VERTEX, (Node*)s2->svertes, (Node*)s2);
        addlist(VERTEX, (Node*)v, (Node*)s1);
    }
    s2->RemoveListFromSolid();
    delete s2;
}

double distancetwovector(Eigen::Vector3d v1, Eigen::Vector3d v2)
{
    double dis = 0;
    dis += (v1[0] - v2[0]) * (v1[0] - v2[0]);
    dis += (v1[1] - v2[1]) * (v1[1] - v2[1]);
    dis += (v1[2] - v2[2]) * (v1[2] - v2[2]);
    return dis;
}

bool match(HalfEdge *h1, HalfEdge *h2)
{
    if (distancetwovector(h1->vtx->vcoord, h2->vtx->vcoord) < 1e-4)
        return true;
    return false;
}

void loopglue(Face *f)
{
    HalfEdge *h1, *h2, *h1next;
    h1 = f->floops->ledg;
    h2 = f->floops->nextl->ledg;
    while (!match(h1, h2))
    {
        h2 = h2->nxt;
    }
    if (h1->wloop == f->flout)
    {
        lmekr(h1, h2);
    }
    else
    {
        lmekr(h2, h1);
    }
    lkev(h1->prv, h2->prv);
    while (h1->nxt != h2)
    {
        h1next = h1->nxt;
        lmef(h1->nxt, h1->prv, -1);
        lkev(h1->nxt, mate(h1->nxt));
        lkef(mate(h1), h1);
        h1 = h1next;
    }
    lkef(mate(h1), h1);
}

void glue(Solid *s1, Solid *s2, Face *f1, Face *f2)
{
    merge(s1, s2);
    lkfmrh(f1, f2);
    loopglue(f1);
}


