#include "euler.h"
// #include <cstdlib>
Id maxv = 0;
Id maxf = 0;

std::shared_ptr<HalfEdge> addhe(std::shared_ptr<Edge> e, std::shared_ptr<Vertex> v, std::shared_ptr<HalfEdge> where, int sign)
{
    std::shared_ptr<HalfEdge> he = nullptr;

    if (where->edg == nullptr)
    {
        he = where;
    }
    else
    {
        he = std::make_shared<HalfEdge>();
        where->prv.lock()->nxt = he;
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

std::shared_ptr<HalfEdge> delhe(std::shared_ptr<HalfEdge> he)
{
    if (he->edg == nullptr)
    {
        // delete he;
        return nullptr;
    }
    else if (he->nxt == he)
    {
        he->edg = nullptr;
        return he;
    }
    else
    {
        std::shared_ptr<HalfEdge> heprv = he->prv.lock();
        heprv->nxt = he->nxt;
        he->nxt->prv = he->prv;
        // HalfEdge *hepre = he->prv;
        // gdel(HALFEDGE, (Node*)he, NIL);
        // delete he;
        return heprv;
    }
}

std::shared_ptr<Solid> getsolid(Id sn)
{
    std::shared_ptr<Solid> s;
    for (s = firsts.lock(); s != nullptr; s = s->nexts)
    {
        if (s->solidno == sn)
            return s;
    }
    return nullptr;
}

std::shared_ptr<Face>  fface(Solid *s, Id fn)
{
    for (std::shared_ptr<Face> f = s->sfaces; f != nullptr; f = f->nextf)
    {
        if (f->faceno == fn)
            return f;
    }
    return nullptr;
}

std::shared_ptr<HalfEdge> fhe(Face *f, Id vn1, Id vn2)
{
    for (Loop *l = f->floops.get(); l != nullptr; l = l->nextl.get())
    {
        std::shared_ptr<HalfEdge> h = l->ledg;
        do
        {
            if (h->vtx->vertexno == vn1 &&
                h->nxt->vtx->vertexno == vn2)
                    return h;
        } while ((h = h->nxt) != l->ledg);
    }
    return nullptr;
}

std::shared_ptr<HalfEdge> fhe(Face *f, Id vn1)
{
    for (Loop *l = f->floops.get(); l != nullptr; l = l->nextl.get())
    {
        std::shared_ptr<HalfEdge> h = l->ledg;
        do
        {
            if (h->vtx->vertexno == vn1)
                    return h;
        } while ((h = h->nxt) != l->ledg);
    }
    return nullptr;
}

std::shared_ptr<Vertex> getvertex(Solid *s, Id vno)
{
    std::shared_ptr<Vertex> v = s->svertes;
    while (v)
    {
        if (v->vertexno == vno)
            return v;
        v = v->nextv;
    }
    return nullptr;
}

std::shared_ptr<Solid> mvfs(Id s, Id f, Id v, double x, double y, double z)
{
    // Solid *newsolid;
    std::shared_ptr<Solid> newsolid = std::make_shared<Solid> ();
    newsolid->addlist();
    std::shared_ptr<Face> newface = std::make_shared<Face> ();
    newsolid->addlist(newface);
    // Face *newface;
    std::shared_ptr<Vertex> newvertex = std::make_shared<Vertex> ();
    newsolid->addlist(newvertex);
    // Vertex *newvertex;
    std::shared_ptr<HalfEdge> newhe = std::make_shared<HalfEdge> ();
    // HalfEdge *newhe;
    std::shared_ptr<Loop> newloop = std::make_shared<Loop> ();
    newface->addlist(newloop);
    // Loop *newloop;
    // newsolid = new Solid();
    // newface = new face(newsolid);
    // newloop = new Loop(newface);
    // newvertex = new Vertex(newsolid);
    // newhe = new HalfEdge();

    newsolid->solidno = s;
    newface->faceno = f;
    newface->flout = newloop;
    newloop->ledg = newhe;
    newhe->wloop = newloop;
    newhe->nxt = newhe;
    newhe->prv = newhe;
    newhe->vtx = newvertex;
    newhe->edg = nullptr;
    newvertex->vertexno = v;
    newvertex->vcoord[0] = x;
    newvertex->vcoord[1] = y;
    newvertex->vcoord[2] = z;
    return newsolid;
}

void lmev(std::shared_ptr<HalfEdge> he1, std::shared_ptr<HalfEdge> he2, Id v, double x, double y, double z)
{
    // HalfEdge *he;
    // Vertex  *newvertex;
    // Edge    *newedge;
    std::shared_ptr<Solid> s = he1->wloop.lock()->lface.lock()->fsolid.lock();
    std::shared_ptr<HalfEdge> he = nullptr;
    std::shared_ptr<Vertex> newvertex = std::make_shared<Vertex>();
    s->addlist(newvertex);
    std::shared_ptr<Edge> newedge = std::make_shared<Edge>();
    s->addlist(newedge);

    newvertex->vcoord[0] = x;
    newvertex->vcoord[1] = y;
    newvertex->vcoord[2] = z;
    newvertex->vertexno = v;

    he = he1;
    while (he != he2)
    {
        he->vtx = newvertex;
        he = he->mate()->nxt;
    }
    
    // the code of book may be wrong (P184)
    addhe(newedge, he2->vtx, he1, MINUS);
    addhe(newedge, newvertex, he2, PLUS);
    newvertex->vedge = he2->prv;
    he2->vtx->vedge = he2;
}

std::shared_ptr<Face> lmef(std::shared_ptr<HalfEdge> he1, std::shared_ptr<HalfEdge> he2, Id f)
{
    std::shared_ptr<HalfEdge> he = nullptr, nhe1 = nullptr, nhe2 = nullptr, temp = nullptr;
    std::shared_ptr<Solid> s = he1->wloop.lock()->lface.lock()->fsolid.lock();
    std::shared_ptr<Edge> newedge = std::make_shared<Edge> ();
    s->addlist(newedge);
    std::shared_ptr<Face> newface = std::make_shared<Face> ();
    s->addlist(newface);
    std::shared_ptr<Loop> newloop = std::make_shared<Loop> ();
    newface->addlist(newloop);

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

    nhe1->prv.lock()->nxt = nhe2;
    nhe2->prv.lock()->nxt = nhe1;
    temp = nhe1->prv.lock();
    nhe1->prv = nhe2->prv;
    nhe2->prv = temp;

    newloop->ledg = nhe1;
    he2->wloop.lock()->ledg = nhe2;
    return newface;
}

void lkef(std::shared_ptr<HalfEdge> he1, std::shared_ptr<HalfEdge> he2)
{
    std::shared_ptr<Edge> edg = he1->edg;
    if (edg != he2->edg)
    {
        std::cerr << "lkef : he1->edg != he2.edg" << std::endl;
        return;
    }
    std::shared_ptr<Loop> he1wloop = he1->wloop.lock();
    std::shared_ptr<Loop> he2wloop = he2->wloop.lock();
    if (he1wloop == he2wloop)
    {
        std::cerr << "lkef : he1->wloop == he2->wloop" << std::endl;      
        return;
    }
    if (he1wloop->lface.lock() == he2wloop->lface.lock())
    {
        std::cerr << "lkef : he1->wloop->lface == he2->wloop->lface" << std::endl;
        return;
    }

    std::shared_ptr<HalfEdge> nxthe1 = he1->nxt;
    std::shared_ptr<HalfEdge> nxthe2 = he2->nxt;

    //TODO: loop only has one halfedge

    //we assume loop has two halfedge at least
    std::shared_ptr<HalfEdge> he = he2->nxt;
    he1wloop->ledg = nxthe1;
    do
    {
        he->wloop = he1->wloop;
    } while ((he = he->nxt) != he2);
    
    he2->vtx->vedge = he1->nxt; 
    he1->vtx->vedge = he2->nxt;
    
    if (he2 == he2->vtx->vedge.lock())
    {
        he2->vtx->vedge = he1->nxt;
    }

    he2->prv.lock()->nxt = he1->nxt;
    he1->nxt->prv = he2->prv;
    he1->nxt = he2;
    he2->prv = he1;

    std::shared_ptr<Solid> s = he1wloop->lface.lock()->fsolid.lock();
    // Loop *dl = he2->wloop;
    std::shared_ptr<Face> df = he2wloop->lface.lock();
    delhe(he1);
    delhe(he2);
    edg->RemoveListFromSolid(s);
    // delete edg;
    he2wloop->RemoveListFromFace(df);
    // delete dl;

    // TODO : consider inner loop
    // gdel(FACE, (Node*)df, (Node*)s);
    df->RemoveListFromSolid(s);
    // delete df;
}

int kef(Id s, Id f, Id v1, Id v2)
{
    // Solid *oldsolid;
    // Face *oldface;
    // Edge *e;
    // HalfEdge *h1, *h2;
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface = nullptr;
    std::shared_ptr<Edge> e = nullptr;
    std::shared_ptr<HalfEdge> h1 = nullptr, h2 = nullptr;
    if ((oldsolid = getsolid(s)) == nullptr)
    {
        std::cerr << "kef : solid " << s << " not found" << std::endl;
        return ERROR; 
    }
    if ((oldface = fface(oldsolid.get(), f)) == NIL)
    {
        std::cerr << "kef : face " << f << " not found in solid " << s << std::endl;
        return ERROR; 
    }
    std::shared_ptr<Loop> l = oldface->floops;
    bool flag = false;
    while(l)
    {
        h1 = oldface->floops->ledg;
        do
        {
            if (h1->vtx->vertexno == v1 && h1->nxt->vtx->vertexno == v2)
            {
                h2 = h1->mate();
                flag = true;
                break;
            }
            else
            {
                if (h1->vtx->vertexno == v2 && h1->nxt->vtx->vertexno == v1)
                {
                   h2 = h1->mate();
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

void lkemr(std::shared_ptr<HalfEdge> h1, std::shared_ptr<HalfEdge> h2)
{
    // HalfEdge *h3, *h4;
    // Loop *nl;
    // Loop *ol;
    // Edge *killedge;

    std::shared_ptr<HalfEdge> h3 = nullptr, h4 = nullptr;
    std::shared_ptr<Loop> ol = nullptr, nl = nullptr;
    std::shared_ptr<Edge> killedge = nullptr;

    ol = h1->wloop.lock();
    std::shared_ptr<Face> fac = ol->lface.lock();
    nl = std::make_shared<Loop> ();
    fac->addlist(nl);
    killedge = h1->edg;

    h3 = h1->nxt;
    h1->nxt = h2->nxt;
    h2->nxt->prv = h1;
    h2->nxt = h3;
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

    if (h3->edg)
        h3->vtx->vedge = h3;
    else
        h3->vtx->vedge.reset();
    // h3->vtx->vedge.lock() = (h3->edg) ? h3 : nullptr;
    if (h4->edg)
        h4->vtx->vedge = h4;
    else
        h4->vtx->vedge.reset();
    // h4->vtx->vedge.lock() = (h4->edg) ? h4 : nullptr;

    killedge->RemoveListFromSolid(fac->fsolid.lock());
    // delete killedge;
}

int kemr(Id s, Id f, Id v1, Id v2)
{
    // Solid *oldsolid;
    // Face *oldface;
    // Edge *e;
    // HalfEdge *h1, *h2;
    
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface = nullptr;
    std::shared_ptr<Edge> e = nullptr;
    std::shared_ptr<HalfEdge> h1 = nullptr, h2 = nullptr;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "kemr : solid " << s << " not found" << std::endl;
        return ERROR; 
    }
    if ((oldface = fface(oldsolid.get(), f)) == NIL)
    {
        std::cerr << "kemr : face " << f << " not found in solid " << s << std::endl;
        return ERROR; 
    }
    std::shared_ptr<Loop> l = oldface->floops;
    // Loop *l = oldface->floops;
    bool flag = false;
    while(l)
    {
        h1 = oldface->floops->ledg;
        do
        {
            if (h1->vtx->vertexno == v1 && h1->nxt->vtx->vertexno == v2)
            {
                h2 = h1->mate();
                flag = true;
                break;
            }
            else
            {
                if (h1->vtx->vertexno == v2 && h1->nxt->vtx->vertexno == v1)
                {
                   h2 = h1->mate();
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

void lmekr(std::shared_ptr<HalfEdge> h1, std::shared_ptr<HalfEdge> h2)
{
    std::shared_ptr<Loop> h1wloop = h1->wloop.lock();
    std::shared_ptr<Loop> h2wloop = h2->wloop.lock();
    if (h1wloop == h2wloop)
    {
        std::cerr << "lmekr : h1->wloop == h2->wloop" << std::endl;
        return;
    }
    std::shared_ptr<Face> h1face = h1wloop->lface.lock();
    std::shared_ptr<Face> h2face = h2wloop->lface.lock();
    if (h1face != h2face)
    {
        std::cerr << "lmekr : h1->wloop->lface != h2->wloop->lface" << std::endl;
        return;
    }

    std::shared_ptr<HalfEdge> he = h2;
    std::shared_ptr<Loop> dl = h2wloop;
    do
    {
        he->wloop = h1->wloop;
    } while ((he = he->nxt) != h2);
    
    std::shared_ptr<Solid> s = h1face->fsolid.lock();
    std::shared_ptr<Edge> newedge = std::make_shared<Edge>();
    s->addlist(newedge);

    // how to ensure h1->wloop is out loop? there code may have bug
    if (h1face->flout == nullptr)
        h1face->flout = h1wloop;
    
    std::shared_ptr<HalfEdge> temp = h1->prv.lock();
    std::shared_ptr<HalfEdge> preh1 = addhe(newedge, h2->vtx, h1, PLUS);
    std::shared_ptr<HalfEdge> preh2 = addhe(newedge, h1->vtx, h2, MINUS);

    preh1->prv = preh2->prv;
    preh2->prv.lock()->nxt = preh1;
    preh2->prv = temp;
    temp->nxt = preh2;

    h2wloop->RemoveListFromFace(h1face);
    // delete dl;
}

int mekr(Id s, Id f, Id v1, Id v2, Id v3, Id v4)
{
    // Solid *oldsolid;
    // Face *oldface1;
    // HalfEdge *he1, *he2;
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface1 = nullptr;
    std::shared_ptr<HalfEdge> he1 = nullptr, he2 = nullptr;
    if ((oldsolid = getsolid(s)) == nullptr)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid.get(), f)) == nullptr)
    {
        std::cerr << "mev : face " << f << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he1 = fhe(oldface1.get(), v1, v2)) == NIL)
    {
        std::cerr << "mev : edge " << v1 << "-" << v2 << " not found in face " << f << std::endl;
        return ERROR;          
    }
    if ((he2 = fhe(oldface1.get(), v3, v4)) == NIL)
    {
        std::cerr << "mev : edge " << v3 << "-" << v4 << " not found in face " << f << std::endl;
        return ERROR;          
    }
    lmekr(he1, he2);
    return SUCCESS;
}

int smekr(Id s, Id f, Id v1, Id v3)
{
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface1 = nullptr;
    std::shared_ptr<HalfEdge> he1 = nullptr, he2 = nullptr;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid.get(), f)) == NIL)
    {
        std::cerr << "mev : face " << f << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he1 = fhe(oldface1.get(), v1)) == NIL)
    {
        std::cerr << "mev : edge " << v1 << "-" << " not found in face " << f << std::endl;
        return ERROR;          
    }
    if ((he2 = fhe(oldface1.get(), v3)) == NIL)
    {
        std::cerr << "mev : edge " << v3 << "-" << " not found in face " << f << std::endl;
        return ERROR;          
    }
    lmekr(he1, he2);
    return SUCCESS;
}

void lkfmrh(std::shared_ptr<Face> fac1, std::shared_ptr<Face> fac2)
{
    // assume fac2 has just one out loop
    std::shared_ptr<Loop> l2 = fac2->floops;
    fac1->addlist(l2);
    // addlist(LOOP, (Node*)l2, (Node*)fac1);
    
    fac2->RemoveListFromSolid(fac2->fsolid.lock());
    // delete fac2;
}

int kfmrh(Id s, Id f1, Id f2)
{
    // Solid *oldsolid;
    // Face *oldface1, *oldface2;
    // HalfEdge *he1, *he2;
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface1 = nullptr, oldface2 = nullptr;
    std::shared_ptr<HalfEdge> he1 = nullptr, he2 = nullptr;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid.get(), f1)) == NIL)
    {
        std::cerr << "mev : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((oldface2 = fface(oldsolid.get(), f2)) == NIL)
    {
        std::cerr << "mev : face " << f2 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    lkfmrh(oldface1, oldface2);
    return ERROR;
}

std::shared_ptr<Face> lmfkrh(std::shared_ptr<Loop> l, Id f)
{
    // assume l is a inner loop
    std::shared_ptr<Face> oldface = l->lface.lock();
    std::shared_ptr<Face> newface = std::make_shared<Face> ();
    oldface->fsolid.lock()->addlist(newface);
    // Face *newface = new Face(l->lface->fsolid);
    newface->faceno = f;
    l->RemoveListFromFace(oldface);
    // removelist(LOOP, (Node*)l, (Node*)l->lface);
    // addlist(LOOP, (Node*)l, (Node*)newface);
    newface->addlist(l);
    newface->flout = l;
    return newface;
}

int mfkrh(Id s, Id f1, Id v1, Id v2, Id f2)
{
    // Solid *oldsolid;
    // Face *oldface1;
    // HalfEdge *he;
    
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface1 = nullptr;
    std::shared_ptr<HalfEdge> he = nullptr;
    if ((oldsolid = getsolid(s)) == nullptr)
    {
        std::cerr << "mfkrh : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid.get(), f1)) == nullptr)
    {
        std::cerr << "mfkrh : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he = fhe(oldface1.get(), v1, v2)) == nullptr)
    {
        std::cerr << "mfkrh : " << v1 << "-" << v2 << " not found in face " << f1 << std::endl;
        return ERROR;
    }
    lmfkrh(he->wloop.lock(), f2);
    return ERROR;
}

void lringmv(std::shared_ptr<Loop> l, std::shared_ptr<Face> toface, int inout)
{
    if (inout != 0)
        toface->flout = l;
    
    // ensure l not in toface
    std::shared_ptr<Loop> fl = toface->floops;
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

int ringmv(std::shared_ptr<Solid> s, Id f1, Id f2, Id v1, Id v2, int inout)
{
    // Face *oldface1, *oldface2;
    // HalfEdge *he;

    // std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface1 = nullptr, oldface2 = nullptr;
    std::shared_ptr<HalfEdge> he = nullptr;
    if ((oldface1 = fface(s.get(), f1)) == NIL)
    {
        std::cerr << "ringmv : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((oldface2 = fface(s.get(), f2)) == NIL)
    {
        std::cerr << "ringmv : face " << f2 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he = fhe(oldface1.get(), v1, v2)) == NIL)
    {
        std::cerr << "mev : " << v1 << "-" << v2 << " not found in face " << f1 << std::endl;
        return ERROR;
    }
    lringmv(he->wloop.lock(), oldface2, inout);
    return SUCCESS;
}

int mev(Id s, Id f1, Id f2, Id v1, Id v2, Id v3, Id v4, double x, double y, double z)
{
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface1 = nullptr, oldface2 = nullptr;
    std::shared_ptr<HalfEdge> he1 = nullptr, he2 = nullptr;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid.get(), f1)) == NIL)
    {
        std::cerr << "mev : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((oldface2 = fface(oldsolid.get(), f2)) == NIL)
    {
        std::cerr << "mev : face " << f2 << " not found in solid " << s << std::endl;
        return ERROR;       
    }
    if ((he1 = fhe(oldface1.get(), v1, v2)) == NIL)
    {
        std::cerr << "mev : edge " << v1 << "-" << v2 << " not found in face " << f1 << std::endl;
        return ERROR;          
    }
    if ((he2 = fhe(oldface2.get(), v1, v3)) == NIL)
    {
        std::cerr << "mev : edge " << v1 << "-" << v3 << " not found in face " << f2 << std::endl;
        return ERROR;          
    }
    lmev(he1, he2, v4, x, y, z);
    return SUCCESS;
}

int smev(Id s, Id f1, Id v1, Id v4, double x, double y, double z)
{
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface1 = nullptr, oldface2 = nullptr;
    std::shared_ptr<HalfEdge> he1 = nullptr, he2 = nullptr;
    if ((oldsolid = getsolid(s)) == nullptr)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return ERROR;
    }
    if ((oldface1 = fface(oldsolid.get(), f1)) == nullptr)
    {
        std::cerr << "mev : face " << f1 << " not found in solid " << s << std::endl;
        return ERROR;
    }
    if ((he1 = fhe(oldface1.get(), v1)) == nullptr)
    {
        std::cerr << "mev : edge " << v1 << "-" << "v2" << " not found in face " << f1 << std::endl;
        return ERROR;          
    }
    lmev(he1, he1, v4, x, y, z);
    return SUCCESS;
}

int smef(Id s, Id f1, Id v1, Id v3, Id f2)
{
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface1 = nullptr, oldface2 = nullptr;
    std::shared_ptr<HalfEdge> he1 = nullptr, he2 = nullptr;
    oldsolid = getsolid(s);
    oldface1 = fface(oldsolid.get(), f1);
    std::shared_ptr<Loop> lp = oldface1->floops;
    while (lp)
    {
        std::shared_ptr<HalfEdge> he = lp->ledg;
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

void lkvfs(std::shared_ptr<Solid> s)
{
    // gdel(SOLID, (Node*)s, NIL);
    s->svertes->RemoveListFromSolid(s);
    s->sfaces->RemoveListFromSolid(s);
    s->RemoveListFromSolid();
    // delete s;
}

void kvfs(Id s)
{
    std::shared_ptr<Solid> oldsolid;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "mev : solid " << s << " not found" << std::endl;    
        return;
    }
    lkvfs(oldsolid);
}

void lkev(std::shared_ptr<HalfEdge> he1, std::shared_ptr<HalfEdge> he2)
{
    if (he1->vtx == he2->vtx)
    {
        std::cerr << "he1->vtx == he2->vtx" << std::endl;
        return;
    }
    // Solid *s = he1->wloop->lface->fsolid;
    // Edge *de = he1->edg;
    // Vertex *dv = he1->vtx;
    // HalfEdge *he = he1;
    std::shared_ptr<Solid> s = he1->wloop.lock()->lface.lock()->fsolid.lock();
    std::shared_ptr<Vertex> dv = he1->vtx;
    std::shared_ptr<Edge> de = he1->edg;
    std::shared_ptr<HalfEdge> he = he1;
    do
    {
        he->vtx = he2->vtx;
    } while ((he = he->mate()->nxt) != he1);
    he2->vtx->vedge = he2->nxt;
    delhe(he1);
    delhe(he2);

    de->RemoveListFromSolid(s);
    // delete de;
    dv->RemoveListFromSolid(s);
    // delete dv;
}

int kev(Id s, Id f, Id v1, Id v2)
{
    std::shared_ptr<Solid> oldsolid = nullptr;
    std::shared_ptr<Face>  oldface = nullptr;
    std::shared_ptr<Edge> e = nullptr;
    std::shared_ptr<HalfEdge> h1 = nullptr, h2 = nullptr;
    if ((oldsolid = getsolid(s)) == NIL)
    {
        std::cerr << "kev : solid " << s << " not found" << std::endl;
        return ERROR; 
    }
    if ((oldface = fface(oldsolid.get(), f)) == NIL)
    {
        std::cerr << "kev : face " << f << " not found in solid " << s << std::endl;
        return ERROR; 
    }
    std::shared_ptr<Loop> l = oldface->floops;
    bool flag = false;
    while(l)
    {
        h1 = oldface->floops->ledg;
        do
        {
            if (h1->vtx->vertexno == v1 && h1->nxt->vtx->vertexno == v2)
            {
                h2 = h1->mate();
                flag = true;
                break;
            }
            else
            {
                if (h1->vtx->vertexno == v2 && h1->nxt->vtx->vertexno == v1)
                {
                   h2 = h1->mate();
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
    // Solid *s;
    std::shared_ptr<Vertex> v = nullptr;
    std::shared_ptr<Face> f = nullptr;

    std::shared_ptr<Solid> s = getsolid(sn);
    for (v = s->svertes, maxv = 0; v != NIL; v = v->nextv)
        if (v->vertexno > maxv) maxv = v->vertexno;
    for (f = s->sfaces, maxf = 0; f != NIL; f = f->nextf)
        if (f->faceno > maxf) maxf = f->faceno;
}

void merge(std::shared_ptr<Solid> s1, std::shared_ptr<Solid> s2)
{
    while (s2->sedges)
    {
        std::shared_ptr<Edge> e = s2->sedges;
        // removelist(EDGE, (Node*)s2->sedges, (Node*)s2);
        s2->sedges->RemoveListFromSolid(s2);
        // addlist(EDGE, (Node*)e, (Node*)s1);
        s1->addlist(e);
    }
    while (s2->sfaces)
    {
        std::shared_ptr<Face> f = s2->sfaces;
        // removelist(FACE, (Node*)s2->sfaces, (Node*)s2);
        s2->sfaces->RemoveListFromSolid(s2);
        s1->addlist(f);
        // addlist(FACE, (Node*)f, (Node*)s1);
    }
    while (s2->svertes)
    {
        std::shared_ptr<Vertex> v = s2->svertes;
        // removelist(VERTEX, (Node*)s2->svertes, (Node*)s2);
        s2->svertes->RemoveListFromSolid(s2);
        s2->addlist(v);
        // addlist(VERTEX, (Node*)v, (Node*)s1);
    }
    s2->RemoveListFromSolid();
    // delete s2;
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
    // if (distancetwovector(h1->vtx->vcoord, h2->vtx->vcoord) < 1e-4)
    //     return true;
    // return false;
    if ((h1->vtx->vcoord - h2->vtx->vcoord).squaredNorm() < EPS * EPS)
        return true;
    return false;
}

void loopglue(std::shared_ptr<Face> f)
{
    std::shared_ptr<HalfEdge> h1 = nullptr, h2 = nullptr, h1next = nullptr;
    h1 = f->floops->ledg;
    h2 = f->floops->nextl->ledg;
    while (!match(h1.get(), h2.get()))
    {
        h2 = h2->nxt;
    }
    if (h1->wloop.lock() == f->flout)
    {
        lmekr(h1, h2);
    }
    else
    {
        lmekr(h2, h1);
    }
    lkev(h1->prv.lock(), h2->prv.lock());
    while (h1->nxt != h2)
    {
        h1next = h1->nxt;
        lmef(h1->nxt, h1->prv.lock(), -1);
        lkev(h1->nxt, h1->nxt->mate());
        lkef(h1->mate(), h1);
        h1 = h1next;
    }
    lkef(h1->mate(), h1);
}

void glue(std::shared_ptr<Solid> s1, std::shared_ptr<Solid> s2, std::shared_ptr<Face> f1, std::shared_ptr<Face> f2)
{
    merge(s1, s2);
    lkfmrh(f1, f2);
    loopglue(f1);
}


