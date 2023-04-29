#include "build.h"
#include <string>
#include <fstream>
#include <Eigen/Dense>

void arc(Id s, Id f, Id v, double cx, double cy, double rad, double h, double phi1, double phi2, int n)
{
    double x, y, angel, inc;
    Id prev;
    int i;
    angel = phi1 * PI / 180.0;
    inc = (phi2 - phi1) * PI / (180.0 * n);
    prev = v;
    getmaxnames(s);
    for (i = 0; i < n; ++i)
    {
        angel += inc;
        x = cx + cos(angel) * rad;
        y = cy + sin(angel) * rad;
        smev(s, f, prev, ++maxv, x, y, h);
        prev = maxv;
    }
}

std::shared_ptr<Solid> circle(Id sn, double cx, double cy, double rad, double h, int n)
{
    std::shared_ptr<Solid> s;
    s = mvfs(sn, 1, 1, cx + rad, cy, h);
    arc(sn, 1, 1, cx, cy, rad, h, 0.0, ((n-1) * 360.0 / n), n-1);
    smef(sn, 1, n, 1, 2);
    return s;
}

void sweep(std::shared_ptr<Face> fac, double dx, double dy, double dz)
{
    // Loop *l;
    // HalfEdge *first, *scan;
    // Vertex *v;
    std::shared_ptr<Loop> l = nullptr;
    std::shared_ptr<HalfEdge> first = nullptr, scan = nullptr;
    std::shared_ptr<Vertex> v = nullptr;
    getmaxnames(fac->fsolid.lock()->solidno);
    l = fac->floops;
    while (l)
    {
        first = l->ledg;
        scan = first->nxt;
        v = scan->vtx;
        lmev(scan, scan, ++maxv, v->vcoord[0] + dx,
            v->vcoord[1] + dy, v->vcoord[2] + dz);
        while (scan != first)
        {
            v = scan->nxt->vtx;
            lmev(scan->nxt, scan->nxt, ++maxv, v->vcoord[0] + dx,
                v->vcoord[1] + dy, v->vcoord[2] + dz);
            lmef(scan->prv.lock(), scan->nxt->nxt, ++maxf);
            scan = scan->nxt->mate()->nxt;
        }
        lmef(scan->prv.lock(), scan->nxt->nxt, ++maxf);
        l = l->nextl;
    }
}

std::shared_ptr<Solid> block(Id sn, double dx, double dy, double dz)
{
    // std::shared_ptr<Solid> s;

    std::shared_ptr<Solid> s = mvfs(sn, 1, 1, 0.0, 0.0, 0.0);
    smev(sn, 1, 1, 2, dx, 0.0, 0.0);
    smev(sn, 1, 2, 3, dx, dy, 0.0);
    smev(sn, 1, 3, 4, 0.0, dy, 0.0);
    smef(sn, 1, 1, 4, 2);
    sweep(fface(s.get(), 2), 0.0, 0.0, dz);
    return s;
}

std::shared_ptr<Solid>  cyl(Id sn, double rad, double h, int n)
{
    std::shared_ptr<Solid> s = circle(sn, 0.0, 0.0, rad, 0.0, n);
    sweep(fface(s.get(), 2), 0.0, 0.0, h);
    return s;
}

std::shared_ptr<Solid> rsweep(std::shared_ptr<Solid> s, int nfaces)
{
    // HalfEdge *first, *cfirst, *last, *scan;
    std::shared_ptr<HalfEdge> first = nullptr, cfirst = nullptr, last = nullptr, scan = nullptr;
    // Face *tailf, *headf;
    std::shared_ptr<Face> tailf = nullptr, headf = nullptr;
    Eigen::Matrix4d m;
    Eigen::Vector3d v;
    int closed_figure = 0;

    getmaxnames(s->solidno);
    if (s->sfaces->nextf)
    {
        closed_figure = 1;
        std::shared_ptr<HalfEdge> h = s->sfaces->floops->ledg;
        Eigen::Vector3d vec = h->vtx->vcoord;
        lmev(h, h->mate()->nxt, ++maxv, vec[0], vec[1], vec[2]);
        lkef(h->prv.lock(), h->prv.lock()->mate());
        headf = s->sfaces;
    }


    first = s->sfaces->floops->ledg;
    while (first->edg != first->nxt->edg)
    {
        first = first->nxt;
    }
    last = first->nxt;
    while (last->edg != last->nxt->edg)
    {
        last = last->nxt;
    }
    cfirst = first;
    Eigen::AngleAxisd rotateVector(2 * PI / (nfaces), Eigen::Vector3d(-1, 0, 0));
    while (--nfaces)
    {
        v = rotateVector * cfirst->nxt->vtx->vcoord;
        lmev(cfirst->nxt, cfirst->nxt, ++maxv, v[0], v[1], v[2]);
        scan = cfirst->nxt;
        while (scan != last->nxt)
        {
            std::shared_ptr<HalfEdge> scanprv = scan->prv.lock();
            v = rotateVector * scanprv->vtx->vcoord;
            lmev(scanprv, scanprv, ++maxv, v[0], v[1], v[2]);
            lmef(scanprv->prv.lock(), scan->nxt, ++maxf);
            scan = scan->nxt->nxt->mate();
        }
        last = scan;
        cfirst = cfirst->nxt->nxt->mate();
    }
    tailf = lmef(cfirst->nxt, first->mate(), ++maxf);
    while (cfirst != scan)
    {
        lmef(cfirst, cfirst->nxt->nxt->nxt, ++maxf);
        cfirst = cfirst->prv.lock()->mate()->prv.lock();
    }
    
    if (closed_figure == 1)
    {
        lkfmrh(headf, tailf);
        loopglue(headf);
    }
    return nullptr;
}

std::shared_ptr<Solid> torus(Id sn, double r1, double r2, int nf1, int nf2)
{
    std::shared_ptr<Solid> s = circle(sn, 0.0, r1, r2, 0.0, nf1);
    // Solid *s;
    // s = mvfs(sn, 1, 1, 0.0 + r2, r1, 0.0);
    // arc(sn, 1, 1, 0.0, r1, r2, 0.0, 0.0, ((nf1-1) * 180.0 / nf1), nf1-1);
    rsweep(s, nf2);
    return s;
}