#include "build.h"
#include <string>
#include <fstream>
// #include  "../../../3rd/Eigen3/include/eigen3/Eigen/Geometry"
// #include "../../../3rd/eigen/Eigen/Geometry"
// #include <../../../3rd/eigen/include/eigen3/Eigen/Geometry>
// #include <Eigen/Geometry>
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

Solid *circle(Id sn, double cx, double cy, double rad, double h, int n)
{
    Solid *s;
    s = mvfs(sn, 1, 1, cx + rad, cy, h);
    arc(sn, 1, 1, cx, cy, rad, h, 0.0, ((n-1) * 360.0 / n), n-1);
    smef(sn, 1, n, 1, 2);
    return s;
}

void sweep(Face *fac, double dx, double dy, double dz)
{
    Loop *l;
    HalfEdge *first, *scan;
    Vertex *v;
    getmaxnames(fac->fsolid->solidno);
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
            lmef(scan->prv, scan->nxt->nxt, ++maxf);
            scan = mate(scan->nxt)->nxt;
        }
        lmef(scan->prv, scan->nxt->nxt, ++maxf);
        l = l->nextl;
    }
}

Solid *block(Id sn, double dx, double dy, double dz)
{
    Solid *s;

    s = mvfs(sn, 1, 1, 0.0, 0.0, 0.0);
    smev(sn, 1, 1, 2, dx, 0.0, 0.0);
    smev(sn, 1, 2, 3, dx, dy, 0.0);
    smev(sn, 1, 3, 4, 0.0, dy, 0.0);
    smef(sn, 1, 1, 4, 2);
    sweep(fface(s, 2), 0.0, 0.0, dz);
    return s;
}

Solid *cyl(Id sn, double rad, double h, int n)
{
    Solid *s = circle(sn, 0.0, 0.0, rad, 0.0, n);
    sweep(fface(s, 2), 0.0, 0.0, h);
    return s;
}

Solid *rsweep(Solid *s, int nfaces)
{
    HalfEdge *first, *cfirst, *last, *scan;
    Face *tailf, *headf;
    Eigen::Matrix4d m;
    Eigen::Vector3d v;
    int closed_figure = 0;

    getmaxnames(s->solidno);
    if (s->sfaces->nextf)
    {
        closed_figure = 1;
        HalfEdge *h = s->sfaces->floops->ledg;
        Eigen::Vector3d vec = h->vtx->vcoord;
        lmev(h, mate(h)->nxt, ++maxv, vec[0], vec[1], vec[2]);
        lkef(h->prv, mate(h->prv));
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
            v = rotateVector * scan->prv->vtx->vcoord;
            lmev(scan->prv, scan->prv, ++maxv, v[0], v[1], v[2]);
            lmef(scan->prv->prv, scan->nxt, ++maxf);
            scan = mate(scan->nxt->nxt);
        }
        last = scan;
        cfirst = mate(cfirst->nxt->nxt);
    }
    tailf = lmef(cfirst->nxt, mate(first), ++maxf);
    while (cfirst != scan)
    {
        lmef(cfirst, cfirst->nxt->nxt->nxt, ++maxf);
        cfirst = mate(cfirst->prv)->prv;
    }

    // Face *fac = s->sfaces;
    // HalfEdge *hf = mate(fac->floops->ledg);
    // headf = s->sfaces->nextf;
    // while (--nfaces)
    // {
    //     Loop *l = fac->floops;
    //     while (l)
    //     {
    //         first = l->ledg;
    //         cfirst = first;
    //         scan = first->nxt;
    //         v = rotateVector * first->nxt->vtx->vcoord;
    //         lmev(scan, scan, ++maxv, v[0], v[1], v[2]);
    //         while (scan != first)
    //         {
    //             v = rotateVector * scan->nxt->vtx->vcoord;
    //             lmev(scan->nxt, scan->nxt, ++maxv, v[0], v[1], v[2]);
    //             lmef(scan->prv, scan->nxt->nxt, ++maxf);
    //             scan = mate(scan->nxt)->nxt;
    //         }
    //         lmef(scan->prv, scan->nxt->nxt, ++maxf);
    //         l = l->nextl;
    //         cfirst = mate(cfirst->nxt->nxt);
    //     }
    // }
    // lkfmrh(headf, fac);
    // loopglue(headf);
    

    if (closed_figure == 1)
    {
        lkfmrh(headf, tailf);
        loopglue(headf);
    }
    return nullptr;
}

Solid *torus(Id sn, double r1, double r2, int nf1, int nf2)
{
    Solid *s = circle(sn, 0.0, r1, r2, 0.0, nf1);
    // Solid *s;
    // s = mvfs(sn, 1, 1, 0.0 + r2, r1, 0.0);
    // arc(sn, 1, 1, 0.0, r1, r2, 0.0, 0.0, ((nf1-1) * 180.0 / nf1), nf1-1);
    rsweep(s, nf2);
    return s;
}
// int main()
// {
//     Solid *s = circle(1, 0, 0, 10, 0, 10);
//     // Solid *s = block(1, 10, 10, 10);
//     // triangel t = delaunay(s);

//     // std::cout << "-----------" << std::endl;

//     // //write doc
//     // std::string dir("/home/yuanbk/Documents/workbench/view.obj");
//     // std::ofstream outfile(dir);

//     // //wirte vetex
//     // auto it = t.vtxarry.begin();
//     // auto end = t.vtxarry.end();
//     // for (; it != end; ++it)
//     // {
//     //     outfile << "v " << (*it)->vcoord[0] << " " <<
//     //         (*it)->vcoord[1] << " " << (*it)->vcoord[2] << std::endl;
//     // }

//     // //write face
//     // auto itf = t.face.begin();
//     // auto endf = t.face.end();
//     // for (; itf != endf; ++itf)
//     // {
//     //     outfile << "f";
//     //     std::vector<int> facevecter = *itf;
//     //     for (auto item : facevecter)
//     //     {
//     //         outfile << " " << item;
//     //     }
//     //     outfile << std::endl;
//     // }

//     // outfile.close();
//     // return 0;
    
//     listsolid(s);
//     std::cout << 1 <<std::endl;
// }