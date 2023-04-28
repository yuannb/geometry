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

HalfEdge* hithe;
Vertex* hitvetex;
// std::unordered_set<Vertex*> soov;
// //int nvtx = 0;
// std::unordered_set<Edge*> sone;
// //int nedg = 0;
// std::vector<Face*> sonf;
//int nfac = 0;

static std::vector<HalfEdge*> ends;

static std::vector<HalfEdge*> loosehe;

// #define ABOVE 1
// #define BELOW -1
// #define ON 0

// struct hefrel
// {
// 	HalfEdge* sector;
// 	int cl;
// };

// std::vector<hefrel> nbr;

// std::vector<Edge*> sortnulledges()
// {
// 	std::vector<Edge*> s;
// 	s.reserve(sone.size());
// 	auto end = sone.end();
// 	for (auto it = sone.begin(); it != end; ++it)
// 	{
// 		s.push_back(*it);
// 	}
// 	std::sort(s.begin(), s.end(), compareVector);
// 	return s;
// }

int faceeq(Loop* l, Eigen::Vector4d& eq)
{
	HalfEdge* he;
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

int contvv(Vertex* v1, Vertex* v2)
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
	HalfEdge* he;
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

int int2ee(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4, int drop, double* t1, double* t2)
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
	HalfEdge* he1, * he2;
	Vertex* v1, * v2, aux;
	double t1, t2;
	int count, intr, c1, c2;

	if ((intr = bndrlv(l, v)) > 0)
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

double svolume(Solid* s)
{
	double res = 0.0;
	Face* f = s->sfaces;
	while (f)
	{
		Loop* l = f->floops;
		while (l)
		{
			HalfEdge* he1 = l->ledg;
			HalfEdge* he2 = he1->nxt;
			do
			{
				Eigen::Vector3d c = he1->vtx->vcoord.cross(he2->vtx->vcoord);
				res += c.dot(he2->nxt->vtx->vcoord);
			} while ((he2 = he2->nxt) != he1);
			l = l->nextl;
		}
		f = f->nextf;
	}
	return res / 6.0;
}

double larea(Loop* l)
{
	Eigen::Vector4d faceq;
	Eigen::Vector4d dd(0, 0, 0, 0);
	faceeq(l, faceq);
	HalfEdge* he = l->ledg;
	Vertex* v = he->vtx;
	he = he->nxt;
	do
	{
		Vertex* v2 = he->vtx;
		Vertex* v3 = he->nxt->vtx;
		Eigen::Vector4d aa, bb;
		aa << v2->vcoord - v->vcoord, 0;
		bb << v3->vcoord - v->vcoord, 0;
		Eigen::Vector4d cc = aa.cross3(bb);
		dd += cc;
	} while ((he = he->nxt) != l->ledg);
	return 0.5 * (dd.dot(faceq));
}

int contfp(Face* f, double x, double y, double z)
{
	Vertex v(x, y, z);
	return contfv(f, &v);
}

void laringmv(Face* f1, Face* f2)
{
	Loop* l = f1->flout;
	Loop* ml = f1->floops;
	while (ml)
	{
		if (l != ml)
		{
			Vertex* v = ml->ledg->vtx;
			if (1 != contlv(l, v, 0))
			{
				lringmv(l, f2, 0);
			}
		}
		ml = ml->nextl;
	}
}

// int checkwideness(HalfEdge* he)
// {
// 	Eigen::Vector3d nor = ((plane*)mate(he)->wloop->lface->surf)->get_normal();
// 	HalfEdge* nxthe = he->nxt;
// 	Eigen::Vector3d v0 = he->vtx->vcoord;
// 	Eigen::Vector3d v1 = nxthe->vtx->vcoord;
// 	Eigen::Vector3d v2 = mate(he)->nxt->nxt->vtx->vcoord;
// 	Eigen::Vector3d v0v1 = v1 - v0;
// 	Eigen::Vector3d v0v2 = v2 - v0;
// 	Eigen::Vector3d c = v0v1.cross(v0v2);
// 	double d = nor.dot(c);
// 	// if (std::abs(d) < EPS)
// 	// 	return ERROR;
// 	return d < -EPS ? 0 : 1;
// }

// void getneighborhood(const Vertex* v, const Eigen::Vector4d& SP)
// {
// 	nbr.clear();
// 	HalfEdge* he = v->vedge;
// 	do
// 	{
// 		double d = dist(he->nxt->vtx->vcoord, SP);
// 		int cl = compare(d, 0.0, EPS);
// 		hefrel rel{ he, cl };
// 		nbr.push_back(rel);
// 		if (checkwideness(he) != 1)
// 		{
// 			Eigen::Vector3d bisec;
// 			bisector(he, bisec);
// 			d = dist(bisec, SP);
// 			int cl2 = compare(d, 0.0, EPS);
// 			hefrel rel2{ he, cl2 };
// 			nbr.push_back(rel2);
// 		}
// 	} while ((he = mate(he)->nxt) != v->vedge);
// }

// void reclassifyonsectors(const Eigen::Vector4d& SP)
// {
// 	auto end = nbr.end();
// 	for (auto it = nbr.begin(); it != end; ++it)
// 	{
// 		Face* f = it->sector->wloop->lface;
// 		Eigen::Vector4d eq = ((plane*)f->surf)->equation;
// 		eq[3] = 0.0;
// 		Eigen::Vector4d c = eq.cross3(SP);
// 		double d = c.squaredNorm();
// 		int cl = compare(d, 0.0, EPS * EPS);
// 		if (cl == 0)
// 		{
// 			d = eq.dot(SP);
// 			int cl2 = compare(d, 0.0, EPS);
// 			if (cl2 == 1)
// 			{
// 				it->cl = BELOW;
// 				auto itnxt = (it + 1) == end ? nbr.begin() : (it + 1);
// 				//may have bug?
// 				itnxt->cl = BELOW;
// 			}
// 			if (cl2 == -1)
// 			{
// 				it->cl = ABOVE;
// 				auto itnxt = (it + 1) == end ? nbr.begin() : (it + 1);
// 				//may have bug?
// 				itnxt->cl = ABOVE;
// 			}
// 			else
// 			{
// 				std::cerr << "reclassifyonsectors : some unkown error occur" << std::endl;
// 			}
// 		}
// 	}
// }

// void reclassifyonedges()
// {
// 	int nnbr = (int)nbr.size();
// 	for (int i = 0; i != nnbr; ++i)
// 	{
// 		if (nbr[i].cl == ON)
// 		{
// 			if (nbr[(nnbr + i - 1) % nnbr].cl == BELOW)
// 			{
// 				if (nbr[(i + 1) % nnbr].cl == BELOW)
// 				{
// 					nbr[i].cl = ABOVE;
// 				}
// 				else
// 				{
// 					nbr[i].cl = BELOW;
// 				}
// 			}
// 			else
// 			{
// 				nbr[i].cl = BELOW;
// 			}
// 		}
// 	}
// }

// void insertnulledges()
// {
// 	int i = 0;
// 	int nnbr = (int)nbr.size();
// 	while (!(nbr[i].cl == BELOW && nbr[(i + 1) % nnbr].cl == ABOVE))
// 	{
// 		if (++i == nnbr) return;
// 	}
// 	int start = i;
// 	//book may have bug
	
// 	while (1)
// 	{
// 		HalfEdge* head = nbr[(i + 1) % nnbr].sector;
// 		if  (nbr[i].sector == nbr[(i + 1) % nnbr].sector)
// 			head = nbr[(i + 2) % nnbr].sector;
// 		while (!(nbr[i].cl == ABOVE && nbr[(i + 1) % nnbr].cl == BELOW))
// 		{
// 			i = (i + 1) % nnbr;
// 		}
// 		HalfEdge* tail = nbr[(i + 1) % nnbr].sector;
// 		if (nbr[i].sector == nbr[(i + 1) % nnbr].sector)
// 			tail = mate(nbr[i].sector)->nxt;
// 		lmev(head, tail, ++maxv, head->vtx->vcoord[0],
// 			head->vtx->vcoord[1],
// 			head->vtx->vcoord[2]);
// 		sone.insert(head->prv->edg);
// 		while (!(nbr[i].cl == BELOW && nbr[(i + 1) % nnbr].cl == ABOVE))
// 		{
// 			i = (i + 1) % nnbr;
// 			if (i == start) return;
// 		}
// 	}
// }

// inline void bisector(HalfEdge* he, Eigen::Vector3d& bisect)
// {
// 	bisect = -(he->vtx->vcoord + he->prv->vtx->vcoord) / 2;
// }

// double dist(const Eigen::Vector3d& v, const Eigen::Vector4d& vec)
// {
// 	// select one point on vec
// 	Eigen::Vector3d point;
// 	Eigen::Vector3d normal(vec[0], vec[1], vec[2]);
// 	double normalLen = normal.squaredNorm();
// 	if (normalLen < EPS * EPS)
// 	{
// 		std::cerr << "SP is not valid" << std::endl;
// 		return 0;
// 	}

// 	normal.normalize();
// 	if ((vec[0] * vec[0]) > 0.4)
// 	{
// 		//let y = z = 0
// 		double x = -vec[3] / vec[0];
// 		point = Eigen::Vector3d(x, 0, 0);
// 	}
// 	else if ((vec[1] * vec[1]) > 0.4)
// 	{
// 		//let x = z = 0
// 		double y = -vec[3] / vec[1];
// 		point = Eigen::Vector3d(0, y, 0);
// 	}
// 	else if ((vec[2] * vec[2]) > 0.4)
// 	{
// 		// let x = y = 0
// 		double z = -vec[3] / vec[2];
// 		point = Eigen::Vector3d(0, 0, z);
// 	}
// 	else
// 	{
// 		// some unkown error occur
// 		std::cerr << "some unkown error occur" << std::endl;
// 		return 0;
// 	}
// 	Eigen::Vector3d dir = v - point;
// 	double distance = dir.dot(normal);
// 	return distance;
// }

// void splitgenerate(Solid* S, Eigen::Vector4d& SP)
// {
// 	//nvtx = 0;
// 	for (Edge* e = S->sedges; e != nullptr; e = e->nexte)
// 	{
// 		Vertex* v1 = e->he1->vtx;
// 		Vertex* v2 = e->he2->vtx;
// 		double d1 = dist(v1->vcoord, SP);
// 		double d2 = dist(v2->vcoord, SP);
// 		int s1 = compare(d1, 0.0, EPS);
// 		int s2 = compare(d2, 0.0, EPS);
// 		if ((s1 == -1 && s2 == 1) || (s1 == 1 && s2 == -1))
// 		{
// 			double t = d1 / (d1 - d2);
// 			Eigen::Vector3d point = v1->vcoord + t * (v2->vcoord - v1->vcoord);
// 			HalfEdge* he = nullptr;
// 			lmev(e->he1, (he = e->he2->nxt), ++maxv, point[0], point[1], point[2]);
// 			soov.insert(he->prv->vtx);
// 		}
// 		else
// 		{
// 			if (s1 == 0) soov.insert(v1);
// 			if (s2 == 0) soov.insert(v2);
// 		}
// 	}
// }

// void splitclassify(Eigen::Vector4d SP)
// {
// 	auto end = soov.end();
// 	int debug = 0;
// 	for (auto it = soov.begin(); it != end; ++it)
// 	{
// 		++debug;
// 		getneighborhood(*it, SP);
// 		reclassifyonsectors(SP);
// 		reclassifyonedges();
// 		insertnulledges();
// 	}
// }

// void splitconnect()
// {
// 	std::vector<Edge*> nulledges = sortnulledges();
// 	auto end = nulledges.end();
// 	// for (auto it = nulledges.begin(); it != end; ++it)
// 	// {
// 	// 	loosehe.push_back((*it)->he1);
// 	// 	loosehe.push_back((*it)->he2);
// 	// }
// 	for (auto it = nulledges.begin(); it != end; ++it)
// 	{
// 		HalfEdge* h1 = canjoine((*it)->he1);
// 		if (h1)
// 		{
// 			join(h1, (*it)->he1);
// 			auto loophalfedgeit = std::find(ends.begin(), ends.end(), mate(h1));
// 			if (loophalfedgeit == ends.end()) cut(h1);
// 		}
// 		HalfEdge* h2 = canjoine((*it)->he2);
// 		if (h2)
// 		{
// 			join(h2, (*it)->he2);
// 			auto loophalfedgeit = std::find(ends.begin(), ends.end(), mate(h2));
// 			if (loophalfedgeit == ends.end()) cut(h2);
// 		}
// 		if (h1 && h2)
// 			cut((*it)->he1);
// 	}
// }

// void splitfinish(Solid* S, Solid** Above, Solid** Below)
// {
// 	int nfac = (int)sonf.size();
// 	for (int i = 0; i != nfac; ++i)
// 	{
// 		Loop *l = sonf[i]->floops;
// 		if (l == sonf[i]->flout)
// 			l = sonf[i]->floops->nextl;
// 		Face* f = lmfkrh(l, ++maxf);
// 		sonf.push_back(f);
// 	}
// 	*Above = new Solid();
// 	*Below = new Solid();
// 	classify(S, *Above, *Below);
// 	cleanup(*Above);
// 	cleanup(*Below);
// }

// void split(Solid* S, Eigen::Vector4d& SP, Solid** Above, Solid** Below)
// {
// 	for (Face* f = S->sfaces; f != nullptr; f = f->nextf)
// 	{
// 		Eigen::Vector4d eq;
// 		faceeq(f->flout, eq);
// 		// Eigen::Vector4d ans = ((plane*)f->surf)->equation;
// 		// ans = eq;
// 		plane *p = new plane(eq);
// 		f->surf = (surface*)p;
// 		// ((plane*)f->surf)->set_eqution(eq);
// 	}
// 	getmaxnames(S->solidno);
// 	splitgenerate(S, SP);
// 	splitclassify(SP);
// 	splitconnect();
// 	splitfinish(S, Above, Below);
// }

int contfv(Face* f, Vertex* v)
{
	Loop* l = f->flout;
	int r = contlv(l, v, 0);
	if (r > 1 || r < 1)
		return r;
	if (r == 1)
	{
		l = f->floops;
		while (l)
		{
			if (l == f->flout)
			{
				l = l->nextl;
				continue;
			}
			r = contlv(l, v, 0);
			if (r == 1)
				return 0;
			if (r > 1)
				return r;
			l = l->nextl;
		}
	}
	return 1;
}

// static bool isin(Edge *e, Edge* egelist)
// {
// 	while (egelist)
// 	{
// 		if (e == egelist) return true;
// 		egelist = egelist->nexte;
// 	}
// 	return false;
// }
// static bool isin(Vertex *v, Vertex* vertexlist)
// {
// 	while (vertexlist)
// 	{
// 		if (v == vertexlist) return true;
// 		vertexlist = vertexlist->nextv;
// 	}
// 	return false;
// }

// void cleanup(Solid* s)
// {
// 	Face* f = s->sfaces;
// 	while (f)
// 	{
// 		Loop* l = f->floops;
// 		while (l)
// 		{
// 			HalfEdge* he = l->ledg;
// 			do
// 			{
// 				Edge* e = he->edg;
// 				Edge* elist = s->sedges;
// 				if (!isin(e, elist))
// 				{
// 					addlist(EDGE, (Node*)e, (Node*)s);
// 				}
// 				Vertex* v = he->vtx;
// 				Vertex* vlist = s->svertes;
// 				if (!isin(v, vlist))
// 				{
// 					addlist(VERTEX, (Node*)v, (Node*)s);
// 				}
// 			} while ((he = he->nxt) != l->ledg);
// 			l = l->nextl;
// 		}
// 		f = f->nextf;
// 	}
// }

// void movefac(Face* f, Solid* s)
// {
// 	removelist(FACE, (Node*)f, (Node*)f->fsolid);
// 	addlist(FACE, (Node*)f, (Node*)s);
// 	Loop* l = f->floops;
// 	while (l)
// 	{
// 		HalfEdge* he = l->ledg;
// 		do
// 		{
// 			Face* f2 = mate(he)->wloop->lface;
// 			if (f2->fsolid != s)
// 				movefac(f2, s);
// 		} while ((he = he->nxt) != l->ledg);
// 		l = l->nextl;
// 	}
// }

// void classify(Solid* S, Solid* Above, Solid* Below)
// {
// 	int nfac = sonf.size() / 2;
// 	for (int i = 0; i != nfac; ++i)
// 	{
// 		movefac(sonf[i], Above);
// 		movefac(sonf[i + nfac], Below);
// 	}
// }

// int neighbor(HalfEdge* h1, HalfEdge* h2)
// {
// 	return (h1->wloop->lface == h2->wloop->lface && (
// 		(h1 == h1->edg->he1 && h2 == h2->edg->he2) || (h1 == h1->edg->he2 && h2 == h2->edg->he1)
// 		));
// }

// void cut(HalfEdge* he)
// {
// 	if (he->edg->he1->wloop == he->edg->he2->wloop)
// 	{
// 		sonf.push_back(he->wloop->lface);
// 		lkemr(he->edg->he1, he->edg->he2);
// 	}
// 	else
// 	{
// 		lkef(he->edg->he1, he->edg->he2);
// 	}
// }

// HalfEdge* canjoine(HalfEdge* he)
// {
// 	auto end = ends.end();
// 	for (auto it = ends.begin(); it != end; ++it)
// 	{
// 		if (neighbor(he, *it))
// 		{
// 			HalfEdge* ret = *it;
// 			ends.erase(it);
// 			// std::remove(loosehe.begin(), loosehe.end(), *it);
// 			// std::remove(loosehe.begin(), loosehe.end(), he);
// 			return ret;
// 		}
// 	}
// 	ends.push_back(he);
// 	return nullptr;
// }

// void join(HalfEdge* h1, HalfEdge* h2)
// {
// 	Face* oldface = h1->wloop->lface;
// 	Face* newface = nullptr;
// 	if (h1->wloop == h2->wloop)
// 	{
// 		if (h1->prv->prv != h2)
// 			newface = lmef(h1, h2->nxt, ++maxf);
// 		else lmekr(h1, h2->nxt);
// 		if (h1->nxt->nxt != h2)
// 		{
// 			lmef(h2, h1->nxt, ++maxf);
// 			if (newface) // && oldface->floops->nextl)
// 				laringmv(oldface, newface);
// 		}
// 	}
// }

// bool compareVector(const Edge* e1, const Edge* e2)
// {
// 	Eigen::Vector3d& v1 = e1->he1->vtx->vcoord;
// 	Eigen::Vector3d &v2 = e2->he1->vtx->vcoord;
// 	int xc = compare(v1[0], v2[0], EPS);
// 	if (xc == -1)
// 		return false;
// 	if (xc == 1)
// 		return true;
// 	if (xc == 0)
// 	{
// 		int yc = compare(v1[1], v2[1], EPS);
// 		if (yc == -1)
// 			return false;
// 		if (yc == 1)
// 			return true;
// 		if (yc == 0)
// 		{
// 			int zc = compare(v1[2], v2[2], EPS);
// 			if (zc == -1)
// 				return false;
// 		}
// 	}
// 	return true;
// }

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
