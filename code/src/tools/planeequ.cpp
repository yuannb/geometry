#include "planeequ.h"
#include "loop.h"
#include "halfedge.h"
#include "vertex.h"
#include "params.h"
#include "curve.h"
#include "edge.h"
#include "face.h"
#include "euler.h"
#include <vector>
#include <unordered_set>
#include "solid.h"

std::shared_ptr<HalfEdge> hithe;
std::shared_ptr<Vertex> hitvetex;


int faceeq(Loop *l, Eigen::Vector4d &eq)
{
	std::shared_ptr<HalfEdge> he = nullptr;
	double a, b, c;
	double xi, yi, zi, xj, yj, zj, xc, yc, zc;
	int len = 0;
	a = b = c = xc = yc = zc = 0.0;
	he = l->ledg;
	do
	{
		xi = he->vtx->vcoord[0];
		yi = he->vtx->vcoord[1];
		zi = he->vtx->vcoord[2];
		xj = he->nxt->vtx->vcoord[0];
		yj = he->nxt->vtx->vcoord[1];
		zj = he->nxt->vtx->vcoord[2];
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
	eq = Eigen::Vector4d(vec[0], vec[1], vec[2], c);
	return SUCCESS;
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

int contvv(Eigen::Vector3d& v1, Eigen::Vector3d& v2)
{
	double diff = (v2 - v1).squaredNorm();
	return compare(diff, 0.0, EPS * EPS);
}

int contvv(Vertex *v1, Vertex *v2)
{
	double diff = (v1->vcoord - v2->vcoord).squaredNorm();
	return compare(diff, 0.0, EPS * EPS);
}

int intrev(Eigen::Vector3d& v1, Eigen::Vector3d& v2, Eigen::Vector3d& v3, double* t)
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

int intrev(Vertex* v1, Vertex* v2, Vertex* v3, double* t)
{
	return intrev(v1->vcoord, v2->vcoord, v3->vcoord, t);
}

int intrev(Edge* edg, Vertex* v, double* t)
{
	Eigen::Vector3d v1 = edg->cur->get_start_point();
	Eigen::Vector3d v2 = edg->cur->get_end_point();
	Eigen::Vector3d v3 = v->vcoord;
	return intrev(v1, v2, v3, t);
}

int contev(Vertex* v1, Vertex* v2, Vertex* v3)
{
	double t;
	if (intrev(v1, v2, v3, &t))
		if (t > (-EPS) && t < (1.0 + EPS))
			return 1;
	return 0;
}

int bndrlv(Loop* l, Vertex* v)
{
	HalfEdge *he = l->ledg.get();
	HalfEdge *first = he;
	do
	{
		if (contvv(he->vtx.get(), v))
		{
			hitvetex = he->vtx;
			hithe = nullptr;
			return 3;
		}
	} while ((he = he->nxt.get()) != first);

	he = first;
	do
	{
		if (contev(he->vtx.get(), he->nxt->vtx.get(), v))
		{
			hithe = he->shared_from_this();
			hitvetex = nullptr;
			return 2;
		}
	} while ((he = he->nxt.get()) != first);

	return -1;
}

int int2ee(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4, double* t1, double* t2)
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
		return 1;
	}
	else
		return 0;
}

int contlv(Loop* l, Vertex* v, int drop)
{
	HalfEdge* he1 = nullptr, * he2 = nullptr;
	Vertex* v1 = nullptr, * v2 = nullptr, aux;
	double t1, t2;
	int count, intr, c1, c2;

	if ((intr = bndrlv(l, v)) > 0)
		return intr;
retry:
	v1 = he2->vtx.get();
	v2 = he2->nxt->vtx.get();
	aux.vcoord = (v1->vcoord + v2->vcoord) / 2.0;
	he1 = l->ledg.get();
	HalfEdge *first = he1;
	count = 0;
	do
	{
		intr = int2ee(v, &aux, v1, v2, &t1, &t2);
		if (intr == 1)
		{
			c1 = compare(t2, 0.0, EPS);
			c2 = compare(t2, 1.0, EPS);
			if (c1 == 0 || c2 == 0)
			{
				he2 = he2->nxt.get();
				if (he2 == first)
					return ERROR;
				goto retry;
			}
			if (c1 == 1 && c2 == -1)
				if (t1 >= 0.0)
					++count;
		}
	} while ((he1 = he1->nxt.get()) != first);

	count %= 2;
	return count;
}

double svolume(Solid* s)
{
	double res = 0.0;
	Face* f = s->sfaces.get();
	while (f)
	{
		Loop* l = f->floops.get();
		while (l)
		{
			HalfEdge* he1 = l->ledg.get();
			HalfEdge* he2 = he1->nxt.get();
			do
			{
				Eigen::Vector3d c = he1->vtx->vcoord.cross(he2->vtx->vcoord);
				res += c.dot(he2->nxt->vtx->vcoord);
			} while ((he2 = he2->nxt.get()) != he1);
			l = l->nextl.get();
		}
		f = f->nextf.get();
	}
	return res / 6.0;
}

double larea(Loop* l)
{
	Eigen::Vector4d faceq;
	Eigen::Vector4d dd(0, 0, 0, 0);
	faceeq(l, faceq);
	HalfEdge* he = l->ledg.get();
	HalfEdge *first = he;
	Vertex* v = he->vtx.get();
	he = he->nxt.get();
	do
	{
		Vertex* v2 = he->vtx.get();
		Vertex* v3 = he->nxt->vtx.get();
		Eigen::Vector4d aa, bb;
		aa << v2->vcoord - v->vcoord, 0;
		bb << v3->vcoord - v->vcoord, 0;
		Eigen::Vector4d cc = aa.cross3(bb);
		dd += cc;
	} while ((he = he->nxt.get()) != first);
	return 0.5 * (dd.dot(faceq));
}

int contfp(Face* f, double x, double y, double z)
{
	Vertex v(x, y, z);
	return contfv(f, &v);
}

void laringmv(std::shared_ptr<Face> f1, std::shared_ptr<Face> f2)
{
	std::shared_ptr<Loop> l = f1->flout;
	std::shared_ptr<Loop> ml = f1->floops;
	while (ml)
	{
		if (l != ml)
		{
			std::shared_ptr<Vertex> v = ml->ledg->vtx;
			if (1 != contlv(l.get(), v.get(), 0))
			{
				lringmv(l, f2, 0);
			}
		}
		ml = ml->nextl;
	}
}


int contfv(Face* f, Vertex* v)
{
	Loop* l = f->flout.get();
	int r = contlv(l, v, 0);
	if (r > 1 || r < 1)
		return r;
	if (r == 1)
	{
		l = f->floops.get();
		while (l)
		{
			if (l == f->flout.get())
			{
				l = l->nextl.get();
				continue;
			}
			r = contlv(l, v, 0);
			if (r == 1)
				return 0;
			if (r > 1)
				return r;
			l = l->nextl.get();
		}
	}
	return 1;
}

bool three_plane_intersect(const Eigen::Vector4d &p1, const Eigen::Vector4d &p2, const Eigen::Vector4d &p3, Eigen::Vector3d &intVtx)
{
	Eigen::Matrix3d mat;
	mat << p1[0], p1[1], p1[2],
		   p2[0], p2[1], p2[2],
		   p3[0], p3[1], p3[2];
	double det = mat.determinant();
	if (det < EPS)
	{
		std::cout << "three_plane_intersect : has no intersect point" << std::endl;
		return false;
	}
	Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
	Eigen::Vector3d vec(p1[3], p2[3], p3[3]);
	intVtx = matSvd.solve(vec);
	return true;
}
