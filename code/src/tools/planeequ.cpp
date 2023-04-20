#include "planeequ.h"
#include "loop.h"
#include "halfedge.h"
#include "vertex.h"
#include "params.h"
#include "curve.h"
#include "edge.h"


HalfEdge *hithe;
Vertex *hitvetex;

int faceeq(Loop *l, Eigen::Vector4d eq)
{
    HalfEdge *he;
    double a, b, c, norm;
    double xi, yi, zi, xj, yj, zj, xc, yc, zc;
    int len;
    a = b = c = xc = yc = zc = 0.0;
    len = 0;
    he = l->ledg;
    do
    {
        xi = he->vtx->vcoord[0];
        yi = he->vtx->vcoord[1];
        zi = he->vtx->vcoord[2];
        a += (yi - yj) * (zi + zj);
        b += (zi - zj) * (xi + xj);
        c += (xi - xj) * (yi + yj);
        xc += xi;
        yc += yi;
        zc += zi;
        ++len;
    } while ((he = he->nxt) != l->ledg);
    Eigen::Vector3d vec(a, b, c);
    if (vec.squaredNorm() < 1e-4)
        return ERROR;
    vec.normalize();
    Eigen::Vector3d k(xc, yc, zc);
    c = (vec.transpose().dot(k)) / (-len);
    eq << vec, c;
}

int compare(double a, double b, double tol)
{
    double delta = std::abs(a - b);
    if (delta < tol)
        return 0;
    else if (a > b)
        return 1;
    else
        return -1;
}

int contvv(Eigen::Vector3d &v1, Eigen::Vector3d &v2)
{
    double diff = (v2 - v1).squaredNorm();
    return compare(diff, 0.0, EPS * EPS);
}

int contvv(Vertex *v1, Vertex *v2)
{
    double diff = (v1->vcoord - v2->vcoord).squaredNorm();
    return compare(diff, 0.0, EPS * EPS);
}

int intrev(Eigen::Vector3d &v1, Eigen::Vector3d &v2, Eigen::Vector3d &v3, double *t)
{
    Eigen::Vector3d v1v2 = v2 - v1;
    double v1v2LenLen = v1v2.squaredNorm();
    if (v1v2LenLen < EPS * EPS)
    {
        *t = 0.0;
        return compare((v2 - v1).squaredNorm(), 0.0, EPS);
    }
    else
    {
        Eigen::Vector3d v1v3 = v3 - v1;
        *t = v1v3.dot(v1v2) / v1v2LenLen;
        Eigen::Vector3d projectPoint = v1 + (*t) * v1v2;
        return contvv(projectPoint, v3);
    }
}

int intrev(Vertex *v1, Vertex *v2, Vertex *v3, double *t)
{
    return intrev(v1->vcoord, v2->vcoord, v3->vcoord, t);
}

int intrev(Edge *edg, Vertex *v, double *t)
{
    Eigen::Vector3d v1 = edg->cur->get_start_point();
    Eigen::Vector3d v2 = edg->cur->get_end_point();
    Eigen::Vector3d v3 = v->vcoord;
    return intrev(v1, v2, v3, t);
}

int contev(Vertex *v1, Vertex *v2, Vertex *v3)
{
    double t;
    if (intrev(v1, v2, v3, &t))
        if (t > (-EPS) && t < (1.0 + EPS))
            return 1;
    return 0;
}

int bndrlv(Loop *l, Vertex *v)
{
    HalfEdge *he;
    he = l->ledg;
    do
    {
        if (contvv(he->vtx, v))
        {
            hitvetex = he->vtx;
            hithe = nullptr;
            return 3;
        }
    } while ((he = he->nxt) != l->ledg);
    
    he = l->ledg;
    do
    {
        if (contev(he->vtx, he->nxt->vtx, v))
        {
            hithe = he;
            hitvetex = nullptr;
            return 2;
        }
    } while ((he = he->nxt) != l->ledg);
    
    return -1;
}

int int2ee(Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, int drop, double *t1, double *t2)
{
    Eigen::Vector3d v1v2 = v2->vcoord - v1->vcoord;
    Eigen::Vector3d v3v4 = v4->vcoord - v3->vcoord;
    Eigen::Matrix<double, 3, 2> mat;
    mat << v1v2, -v3v4;

    Eigen::Vector3d vec = v3->vcoord - v1->vcoord;
    Eigen::Matrix<double, 3, 3> externMat;
    externMat << v1v2, v3v4, vec;
    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> externMatSvd(externMat);
    int rankMat = matSvd.rank();
    int rankExternMat = externMatSvd.rank();
    if (rankMat == 2 && rankExternMat == 2)
    {
        Eigen::Vector2d v = matSvd.solve(vec);
        *t1 = v[0];
        *t2 = v[1];
    }
    else
        return 0;
}

int contlv(Loop *l, Vertex *v, int drop)
{
    HalfEdge *he1, *he2;
    Vertex *v1, *v2, aux;
    double t1, t2;
    int count, intr, c1, c2;

    if ((intr == bndrlv(l, v)) > 0)
        return intr;
retry:
    v1 = he2->vtx;
    v2 = he2->nxt->vtx;
    aux.vcoord = (v1->vcoord + v2->vcoord) / 2.0;
    he1 = l->ledg;
    count = 0;
    do
    {
        intr = int2ee(v, &aux, v1, v2, drop, &t1, &t2);
        if (intr == 1)
        {
            c1 = compare(t2, 0.0, EPS);
            c2 = compare(t2, 1.0, EPS);
            if (c1 == 0 || c2 == 0)
            {
                he2 = he2->nxt;
                if (he2 == l->ledg)
                    return ERROR;
                goto retry;
            }
            if (c1 == 1 && c2 == -1)
                if (t1 >= 0.0)
                    ++count;
        }
    } while ((he1 = he1->nxt) != l->ledg);
    
    count %= 2;
    return count;
}
