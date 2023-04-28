#include "split.h"
#include "face.h"
#include "solid.h"
#include "loop.h"
#include "edge.h"
#include "halfedge.h"
#include "iostream"
#include "params.h"
#include "plane.h"
#include "planeequ.h"

// static std::vector<HalfEdge*> ends;
static std::map<Vertex*, std::vector<Edge*>> nulledges;
// std::unordered_set<HalfEdge*> ends;
std::unordered_set<HalfEdge*> hasconnect;

std::vector<Face*> sonf;

double dist(const Eigen::Vector3d& v, const Eigen::Vector4d& vec)
{
	// select one point on vec
	Eigen::Vector3d point;
	Eigen::Vector3d normal(vec[0], vec[1], vec[2]);
	double normalLen = normal.squaredNorm();
	if (normalLen < EPS * EPS)
	{
		std::cerr << "SP is not valid" << std::endl;
		return 0;
	}

	normal.normalize();
	if ((vec[0] * vec[0]) > 0.4)
	{
		//let y = z = 0
		double x = -vec[3] / vec[0];
		point = Eigen::Vector3d(x, 0, 0);
	}
	else if ((vec[1] * vec[1]) > 0.4)
	{
		//let x = z = 0
		double y = -vec[3] / vec[1];
		point = Eigen::Vector3d(0, y, 0);
	}
	else if ((vec[2] * vec[2]) > 0.4)
	{
		// let x = y = 0
		double z = -vec[3] / vec[2];
		point = Eigen::Vector3d(0, 0, z);
	}
	else
	{
		// some unkown error occur
		std::cerr << "some unkown error occur" << std::endl;
		return 0;
	}
	Eigen::Vector3d dir = v - point;
	double distance = dir.dot(normal);
	return distance;
};

int neighbor(HalfEdge* h1, HalfEdge* h2)
{
	return (h1->wloop->lface == h2->wloop->lface && (
		(h1 == h1->edg->he1 && h2 == h2->edg->he2) || (h1 == h1->edg->he2 && h2 == h2->edg->he1)
		));
}

void cut(HalfEdge* he)
{
	if (he->edg->he1->wloop == he->edg->he2->wloop)
	{
		sonf.push_back(he->wloop->lface);
		lkemr(he->edg->he1, he->edg->he2);
	}
	else
	{
		lkef(he->edg->he1, he->edg->he2);
	}
}

HalfEdge* canjoine(HalfEdge* he, std::unordered_set<Edge*> nulledgs)
{
	size_t cnt = hasconnect.count(he);
	if (cnt == 1)
		return nullptr;
	auto end = nulledgs.end();
	for (auto it = nulledgs.begin(); it != end; ++it)
	{
		if (neighbor(he, (*it)->he1))
		{
			hasconnect.insert((*it)->he1);
			hasconnect.insert(he);
			return (*it)->he1;
		}
		if (neighbor(he, (*it)->he2))
		{
			hasconnect.insert((*it)->he2);
			hasconnect.insert(he);
			return (*it)->he2;
		}
	}
	return nullptr;
}

void join(HalfEdge* h1, HalfEdge* h2)
{
	Face* oldface = h1->wloop->lface;
	Face* newface = nullptr;
	if (h1->wloop == h2->wloop)
	{
		if (h1->prv->prv != h2)
			newface = lmef(h1, h2->nxt, ++maxf);
	}
	else lmekr(h1, h2->nxt);
	if (h1->nxt->nxt != h2)
	{
		lmef(h2, h1->nxt, ++maxf);
		if (newface) // && oldface->floops->nextl)
			laringmv(oldface, newface);
	}
}

bool comparevertex(intVertex *v1, intVertex *v2, Eigen::Vector3d dir)
{
	Eigen::Vector3d v1v2 = v2->pos - v1->pos;
	if (v1v2.dot(dir) >= -EPS)
		return true;
	return false;
}

bool sortintvertex(intEdge* intEdg, Eigen::Vector3d dir)
{	
	std::sort(intEdg->stl.begin(), intEdg->stl.end(), [&](intVertex *v1, intVertex *v2) {return (v2->pos - v1->pos).dot(dir) > 0; });
	std::sort(intEdg->edl.begin(), intEdg->edl.end(), [&](intVertex *v1, intVertex *v2) {return (v2->pos - v1->pos).dot(dir) > 0; });

	//debug : chekcout vailid
	size_t stcount = intEdg->stl.size();
	size_t edcount = intEdg->edl.size();
	if (stcount != edcount)
	{
		std::cerr << "sortintvertex : stl number != edl number" << std::endl;
	}
	for (size_t i = 0; i != stcount; ++i)
	{
		if (i == 0)
		{
			if (!comparevertex(intEdg->stl[0], intEdg->edl[0], dir))
			{
				std::cerr << "sortintvertex : intEdg->stl[0] >  intEdg->edl[0]" << std::endl;
			}
		}
		else if (i == stcount - 1) //除掉了stcount == 1的情况
		{
			if (!comparevertex(intEdg->edl[i - 1], intEdg->stl[i], dir))
			{
				std::cerr << "sortintvertex : intEdg->edl[i - 1] > intEdg->stl[i]" << std::endl;
			}
		}
	}
	return true;
}

std::map<intVertex*, std::vector<intVertex*>> makeconnectgraphic(std::map<vertexName, intVertex*> &intVertexMap, 
												std::map<edgeName, intEdge*> &intEdgeMap, Solid *S, Eigen::Vector3d normal)
{
	std::map<intVertex*, std::vector<intVertex*>> ans;
	//sort intvertex
	auto end = intEdgeMap.end();
	for (auto it = intEdgeMap.begin(); it != end; ++it)
	{
		Eigen::Vector3d v1 = S->getface(it->second->intEdgeName.first.first)->get_normal();
		Eigen::Vector3d dir = v1.cross(normal);
		sortintvertex(it->second, dir);
		size_t vertexcount = it->second->stl.size();
		for (size_t index = 0; index != vertexcount; ++index)
		{
			auto ansit = ans.find(it->second->stl[index]);// std::find(ans.begin(), ans.end(), it->second->stl[index]);
			if (ansit != ans.end())
			{
				ansit->second.push_back(it->second->edl[index]);
			}
			else
			{
				std::vector<intVertex*> endpoint{it->second->edl[index]};
				ans.emplace(it->second->stl[index], endpoint);
			}
		}
	}
	return ans;
}

std::unordered_set<Vertex*> splitgenerate(Solid* S, Eigen::Vector4d& SP)
{
	//nvtx = 0;
    std::unordered_set<Vertex*> soov;
	for (Edge* e = S->sedges; e != nullptr; e = e->nexte)
	{
		Vertex* v1 = e->he1->vtx;
		Vertex* v2 = e->he2->vtx;
		double d1 = dist(v1->vcoord, SP);
		double d2 = dist(v2->vcoord, SP);
		int s1 = compare(d1, 0.0, EPS);
		int s2 = compare(d2, 0.0, EPS);
		if ((s1 == -1 && s2 == 1) || (s1 == 1 && s2 == -1))
		{
			double t = d1 / (d1 - d2);
			Eigen::Vector3d point = v1->vcoord + t * (v2->vcoord - v1->vcoord);
			HalfEdge* he = nullptr;
			lmev(e->he1, (he = e->he2->nxt), ++maxv, point[0], point[1], point[2]);
			soov.insert(he->prv->vtx);
		}
		else
		{
			if (s1 == 0) soov.insert(v1);
			if (s2 == 0) soov.insert(v2);
		}
	}
    return soov;
}

int checkwideness(HalfEdge* he)
{
	Eigen::Vector3d nor = ((plane*)mate(he)->wloop->lface->surf)->get_normal();
	HalfEdge* nxthe = he->nxt;
	Eigen::Vector3d v0 = he->vtx->vcoord;
	Eigen::Vector3d v1 = nxthe->vtx->vcoord;
	Eigen::Vector3d v2 = he->prv->vtx->vcoord;
	Eigen::Vector3d v0v1 = v1 - v0;
	Eigen::Vector3d v2v0 = v0 - v2;
	Eigen::Vector3d c = v2v0.cross(v0v1);
	double d = nor.dot(c);
	// if (std::abs(d) > EPS)
	// {
	// 	Eigen::Vector3d bisec = (v0v1 + v0v2) / 2.0;
	// 	Vertex *v = new Vertex(bisec[0], bisec[1], bisec[2]);;
	// 	int flag = contlv(he->wloop, v, 0);
	// 	if (flag != 1)
	// 		return 0;
	// 	return 1;
	// }
	// if (std::abs(d) < EPS)
	// 	return ERROR;
	return d < -EPS ? 0 : 1;
	// return 1;
}
inline void bisector(HalfEdge* he, Eigen::Vector3d& bisect)
{
	bisect = -(he->vtx->vcoord + he->prv->vtx->vcoord) / 2;
};

std::vector<hefrel> getneighborhood(const Vertex* v, const Eigen::Vector4d& SP)
{
	std::vector<hefrel> nbr;
	HalfEdge* he = v->vedge;
	do
	{
		double d = dist(he->nxt->vtx->vcoord, SP);
		int cl = compare(d, 0.0, EPS);
		hefrel rel{ he, cl };
		nbr.push_back(rel);
		if (checkwideness(he) != 1)
		{
			Eigen::Vector3d bisec;
			bisector(he, bisec);
			d = dist(bisec, SP);
			int cl2 = compare(d, 0.0, EPS);
			hefrel rel2{ he, cl2 };
			nbr.push_back(rel2);
		}
	} while ((he = mate(he)->nxt) != v->vedge);
    return nbr;
}

void reclassifyonsectors(const Eigen::Vector4d& SP, std::vector<hefrel> &nbr)
{
	auto end = nbr.end();
	for (auto it = nbr.begin(); it != end; ++it)
	{
		Face* f = it->sector->wloop->lface;
		Eigen::Vector4d eq = ((plane*)f->surf)->equation;
		eq[3] = 0.0;
		Eigen::Vector4d c = eq.cross3(SP);
		double d = c.squaredNorm();
		int cl = compare(d, 0.0, EPS * EPS);
		if (cl == 0)
		{
			d = eq.dot(SP);
			int cl2 = compare(d, 0.0, EPS);
			if (cl2 == 1)
			{
				it->cl = BELOW;
				auto itnxt = (it + 1) == end ? nbr.begin() : (it + 1);
				//may have bug?
				itnxt->cl = BELOW;
			}
			if (cl2 == -1)
			{
				it->cl = ABOVE;
				auto itnxt = (it + 1) == end ? nbr.begin() : (it + 1);
				//may have bug?
				itnxt->cl = ABOVE;
			}
			else
			{
				std::cerr << "reclassifyonsectors : some unkown error occur" << std::endl;
			}
		}
	}
}

void reclassifyonedges(std::vector<hefrel> &nbr)
{
	int nnbr = (int)nbr.size();
	for (int i = 0; i != nnbr; ++i)
	{
		if (nbr[i].cl == ON)
		{
			if (nbr[(nnbr + i - 1) % nnbr].cl == BELOW)
			{
				if (nbr[(i + 1) % nnbr].cl == BELOW)
				{
					nbr[i].cl = ABOVE;
				}
				else
				{
					nbr[i].cl = BELOW;
				}
			}
			else
			{
				nbr[i].cl = BELOW;
			}
		}
	}
}

int findfirstwidesector(std::vector<hefrel> &nbr)
{
    int count = (int)nbr.size();
    for (int index = 0; index != count; ++index)
    {
        if (nbr[index].sector == nbr[index + 1].sector)
            return index;
    }
    return 0;
}

bool createintvertex(HalfEdge *he, std::map<vertexName, intVertex*> &intVertexMap, std::map<edgeName, intEdge*> &intEdgeMap, Eigen::Vector4d SP)
{
	faceSolidMap f0No{0, 0};
	int soildN0 = he->edg->he1->wloop->lface->fsolid->solidno;
	Face *f1 = he->wloop->lface;
	Face *f2 = mate(he)->wloop->lface;
	// Face *f1 = he->edg->he1->wloop->lface;
	// Face *f2 = he->edg->he2->wloop->lface;
	faceSolidMap f1No{f1->faceno, soildN0};
	faceSolidMap f2No{f2->faceno, soildN0};
	vertexName vn;
	edgeName e1;
	edgeName e2;
	Eigen::Vector3d pos = he->vtx->vcoord;
	if (f1 == f2)
	{
		//处理两个面平行的情况
		vn = vertexName(f0No, f1No, f2No);
		e1 = edgeName(f2No, f0No);
		e2 = edgeName(f1No, f0No); 
	}
	Eigen::Vector3d f0Nor(SP[0], SP[1], SP[2]);
	Eigen::Vector3d f1Nor = f1->get_normal();
	Eigen::Vector3d f2Nor = f2->get_normal();
	Eigen::Vector3d c = f1Nor.cross(f2Nor);
	
	Eigen::Vector3d testv = he->nxt->nxt->vtx->vcoord - he->vtx->vcoord;
	double x  = testv.dot(c);
	int check = x < -EPS ? -1 : 1; //testv.dot(c) > 0 ? 1 : -1;

	double mixproduct = f0Nor.dot(c) * check;
	if (std::abs(mixproduct) < EPS)
	{
		//不应该出现这种情况
		std::cerr << "createintvertex: null edge on split plane" << std::endl;
		return false;
	}
	else if(mixproduct > EPS)
	{
		vn = vertexName(f0No, f1No, f2No);
		e1 = edgeName(f2No, f0No);
		e2 = edgeName(f1No, f0No);
	}
	else
	{
		vn = vertexName(f0No, f2No, f1No);
		e1 = edgeName(f1No, f0No);
		e2 = edgeName(f2No, f0No);	
	}
	setintVertexintEdgerelation(vn, e1, true, pos, he->vtx->vertexno, intVertexMap, intEdgeMap);
	setintVertexintEdgerelation(vn, e2, false, pos, he->vtx->vertexno, intVertexMap, intEdgeMap);
	return true;
}

bool insertnulledges(std::vector<hefrel> &nbr, std::map<vertexName, intVertex*> &intVertexMap, std::map<edgeName, intEdge*> &intEdgeMap, Eigen::Vector4d SP)
{
    int i = 0;
	int nnbr = (int)nbr.size();
	while (!(nbr[i].cl == BELOW && nbr[(i + 1) % nnbr].cl == ABOVE))
	{
		if (++i == nnbr) return true;
	}
	int start = i;
	//book may have bug
	std::vector<Edge*> nulledgs;
	while (1)
	{
		HalfEdge* head = nbr[(i + 1) % nnbr].sector;
		if  (nbr[i].sector == nbr[(i + 1) % nnbr].sector)
			head = nbr[(i + 2) % nnbr].sector;
		while (!(nbr[i].cl == ABOVE && nbr[(i + 1) % nnbr].cl == BELOW))
		{
			i = (i + 1) % nnbr;
		}
		HalfEdge* tail = nbr[(i + 1) % nnbr].sector;
		if (nbr[i].sector == nbr[(i + 1) % nnbr].sector)
			tail = mate(nbr[i].sector)->nxt;
		lmev(head, tail, ++maxv, head->vtx->vcoord[0],
			head->vtx->vcoord[1],
			head->vtx->vcoord[2]);
		// sone.insert(head->prv->edg);
		nulledgs.push_back(head->prv->edg);
		createintvertex(head->prv, intVertexMap, intEdgeMap, SP);
		while (!(nbr[i].cl == BELOW && nbr[(i + 1) % nnbr].cl == ABOVE))
		{
			i = (i + 1) % nnbr;
			if (i == start) 
			{
				if (!nulledgs.empty())
				{
					nulledges[head->prv->vtx] = nulledgs;
				}
				return true;
			}
		}
	}
	return true;
}

bool splitclassify(Eigen::Vector4d SP, std::unordered_set<Vertex*> &soov, std::map<vertexName, intVertex*> &intVertexMap, std::map<edgeName, intEdge*> &intEdgeMap)
{
    auto end = soov.end();
	for (auto it = soov.begin(); it != end; ++it)
	{
		std::vector<hefrel> nbr = getneighborhood(*it, SP);
		reclassifyonsectors(SP, nbr);
		reclassifyonedges(nbr);
        insertnulledges(nbr, intVertexMap, intEdgeMap, SP);
        // createintvertex()
        // creatintVertex(nbr)
		// insertnulledges();
	}
}

bool splitconnect(std::map<intVertex*, std::vector<intVertex*>> wire, Solid *S)
{
	auto wireend = wire.end();
	for (auto wireit = wire.begin(); wireit != wireend; ++wireit)
	{
		Vertex *v = S->getvertex(wireit->first->vertexno);
		std::vector<Edge*> edgs = nulledges[v];
		auto edgsend = edgs.end();
		std::unordered_set<Edge*> canjoineedge;
		auto intvetxend = wireit->second.end();
		for (auto intvetxeit = wireit->second.begin(); intvetxeit != intvetxend; ++intvetxeit)
		{
			Vertex *vs = S->getvertex((*intvetxeit)->vertexno);
			canjoineedge.insert(nulledges[vs].begin(), nulledges[vs].end());
		}
		for (auto edgsit = edgs.begin(); edgsit != edgsend; ++edgsit)
		{
			HalfEdge* h1 = canjoine((*edgsit)->he1, canjoineedge);
			if (h1)
			{
				join(h1, (*edgsit)->he1);
				size_t cnt = hasconnect.count(mate(h1));
				if (cnt == 1) cut(h1);
				cnt = hasconnect.count(mate((*edgsit)->he1));
				if (cnt == 1) cut((*edgsit)->he1);
				// auto loophalfedgeit = std::find(ends.begin(), ends.end(), mate(h1));
				// if (loophalfedgeit == ends.end()) cut(h1);
			}
			HalfEdge* h2 = canjoine((*edgsit)->he2, canjoineedge);
			if (h2)
			{
				join(h2, (*edgsit)->he2);
				size_t cnt = hasconnect.count(mate(h2));
				if (cnt == 1) cut(h2);
				cnt = hasconnect.count(mate((*edgsit)->he2));
				if (cnt == 1) cut((*edgsit)->he2);
				// auto loophalfedgeit = std::find(ends.begin(), ends.end(), mate(h2));
				// if (loophalfedgeit == ends.end()) cut(h2);
			}
			if (h1 && h2)
			{
				std::cerr << "splitconnect : some unkown error occured" << std::endl;
			}
			// bool flag1 = ends.count(h1);
			// bool flag2 = ends.count(h2);
			// if (flag1)
			// 	cut((*edgsit)->he1);
			// if (flag2)
			// 	cut((*edgsit)->he2);
		}
	}
	return true;
	
}

void movefac(Face* f, Solid* s)
{
	// Face *f1 = s->sfaces;
	// int count = 0;
	// while (f1 && f1->nextf)
	// {
	// 	++count;
	// 	f1 = f1->nextf;
	// }
	// if (f1 && f1->faceno != 10204)
	// {
	// 	std::cerr << "--------" << count << "-------" << std::endl;
	// }
	removelist(FACE, (Node*)f, (Node*)f->fsolid);
	// f1 = s->sfaces;
	// while (f1 && f1->nextf)
	// {
	// 	f1 = f1->nextf;
	// }
	// if (f1 && f1->faceno != 10204)
	// {
	// 	std::cerr << "---------------" << std::endl;
	// }
	// count = 0;
	// f1 = s->sfaces;
	// while (f1 && f1->nextf)
	// {
	// 	++count;
	// 	f1 = f1->nextf;
	// }
	// if (f1 && f1->faceno != 10204)
	// {
	// 	std::cerr << "--------" << count << "-------" << std::endl;
	// }
	addlist(FACE, (Node*)f, (Node*)s);
	// count = 0;
	// f1 = s->sfaces;
	// while (f1 && f1->nextf)
	// {
	// 	++count;
	// 	f1 = f1->nextf;
	// }
	// if (f1 && f1->faceno != 10204)
	// {
	// 	std::cerr << "--------" << count << "-------" << std::endl;
	// }
	// f1 = s->sfaces;
	// while (f1 && f1->nextf)
	// {
	// 	f1 = f1->nextf;
	// }
	// if (f1 && f1->faceno != 10204)
	// {
	// 	std::cerr << "---------------" << std::endl;
	// }
	Loop* l = f->floops;
	while (l)
	{
		HalfEdge* he = l->ledg;
		do
		{
			Face* f2 = mate(he)->wloop->lface;
			if (f2->fsolid != s)
				movefac(f2, s);
		} while ((he = he->nxt) != l->ledg);
		l = l->nextl;
	}
}

void classify(Solid* S, Solid* Above, Solid* Below)
{
	int nfac = sonf.size() / 2;
	for (int i = 0; i != nfac; ++i)
	{
		//暂时不支持内环，稍微修改一下即可
		movefac(sonf[i], Above);
		movefac(sonf[i + nfac], Below);
	//debug
	// Face *af1 = Above->sfaces;
	// Face *f = af1;
	// while (f->nextf)
	// {
	// 	f = f->nextf;
	// }
	// Face *af2 = f;
	// Face *bf1 = Below->sfaces;
	// f = bf1;
	// while (f->nextf)
	// {
	// 	f = f->nextf;
	// }
	// Face *bf2 = f;
	}

}

static bool isin(Edge *e, Edge* egelist)
{
	while (egelist)
	{
		if (e == egelist) return true;
		egelist = egelist->nexte;
	}
	return false;
}
static bool isin(Vertex *v, Vertex* vertexlist)
{
	while (vertexlist)
	{
		if (v == vertexlist) return true;
		vertexlist = vertexlist->nextv;
	}
	return false;
}

void cleanup(Solid* s)
{
	Face* f = s->sfaces;
	while (f)
	{
		Loop* l = f->floops;
		while (l)
		{
			HalfEdge* he = l->ledg;
			do
			{
				Edge* e = he->edg;
				Edge* elist = s->sedges;
				if (!isin(e, elist))
				{
					addlist(EDGE, (Node*)e, (Node*)s);
				}
				Vertex* v = he->vtx;
				Vertex* vlist = s->svertes;
				if (!isin(v, vlist))
				{
					addlist(VERTEX, (Node*)v, (Node*)s);
				}
			} while ((he = he->nxt) != l->ledg);
			l = l->nextl;
		}
		f = f->nextf;
	}
}

void splitfinish(Solid* S, Solid** Above, Solid** Below)
{
	int nfac = (int)sonf.size();
	for (int i = 0; i != nfac; ++i)
	{
		Loop *l = sonf[i]->floops;
		if (l == sonf[i]->flout)
			l = sonf[i]->floops->nextl;
		Face* f = lmfkrh(l, ++maxf);
		sonf.push_back(f);
	}
	*Above = new Solid();
	*Below = new Solid();
	classify(S, *Above, *Below);
	cleanup(*Above);
	cleanup(*Below);
}

void split(Solid* S, Eigen::Vector4d& SP, Solid** Above, Solid** Below)
{
    std::map<vertexName, intVertex*> intVertexMap;
    std::map<edgeName, intEdge*> intEdgeMap;
	for (Face* f = S->sfaces; f != nullptr; f = f->nextf)
	{
		Eigen::Vector4d eq;
		faceeq(f->flout, eq);
		plane *p = new plane(eq);
		f->surf = (surface*)p;
	}
	getmaxnames(S->solidno);
	std::unordered_set<Vertex*> soov = splitgenerate(S, SP);
    splitclassify(SP, soov, intVertexMap, intEdgeMap);
	Eigen::Vector3d normal(SP[0], SP[1], SP[2]);
	std::map<intVertex*, std::vector<intVertex*>> wire = makeconnectgraphic(intVertexMap, intEdgeMap, S, normal);
	splitconnect(wire, S);
	splitfinish(S, Above, Below);
}
// segments::~segments()
// {
// }

void setintVertexintEdgerelation(vertexName vtxName, edgeName edgName, bool inout, Eigen::Vector3d &point, Id vertexno,
									std::map<vertexName, intVertex*> &intVertexMap, std::map<edgeName, intEdge*> &intEdgeMap)
{
    auto vtxIt = intVertexMap.find(vtxName);
    auto edgIt = intEdgeMap.find(edgName);
    intVertex *intVtx = nullptr;
    intEdge *intEdg = nullptr;
    if (vtxIt == intVertexMap.end())
    {
        intVtx = new intVertex{vtxName, point, vertexno};
        intVertexMap.emplace(vtxName, intVtx);
    }
    else
    {
        intVtx = vtxIt->second;
    }

    if (edgIt == intEdgeMap.end())
    {
        intEdg = new intEdge{edgName};
        intEdgeMap.emplace(edgName, intEdg);
    }
    else
    {
        intEdg = edgIt->second;
    }

    if (inout)
    {
        //edge是vertex的入边
        intEdg->edl.push_back(intVtx);
        intVtx->inEdge.push_back(intEdg);
    }
    else
    {
        //edge是vertex的出边
        intEdg->stl.push_back(intVtx);
        intVtx->outEdge.push_back(intEdg);
    }
}

// onFaceVertex* segments::get_on_face_vertex(const Eigen::Vector4d &SP, Face *f2, Face *f3)
// {
//     faceSolidMap fs1{0, 0};
//     faceSolidMap fs2{f2->faceno, f2->fsolid->solidno};
//     faceSolidMap fs3{f3->faceno, f3->fsolid->solidno};

//     Eigen::Vector3d face2Normal = f2->get_normal();
//     Eigen::Vector3d face3Normal = f3->get_normal();
//     Eigen::Vector3d face2face3Cross = face2Normal.cross(face3Normal);
//     face2face3Cross.normalize();
//     if (face2face3Cross.squaredNorm() < EPS * EPS)
//     {
//         // 暂时不支持边的相邻两个面平行(或者根本不用支持)
//         // std::cerr << "get_intVertex : f1 is parallel with f2" << std::endl;
//         return;
//     }
//     Eigen::Vector3d spNormal{SP[0], SP[1], SP[2]};
//     double mixproduct = spNormal.dot(face1face2Cross);
//     faceSolidMap f1M{f1->faceno, s->solidno};
//     faceSolidMap f2M{f2->faceno, s->solidno};
//     edgeName inEdgeName;
//     edgeName outEdgeName;
//     vertexName vtxName;
//     if (mixproduct > EPS)
//     {
//         vtxName = vertexName(f0M, f1M, f2M);
//         inEdgeName = edgeName(f2M, f0M);
//         outEdgeName = edgeName(f1M, f0M);
//     }
//     else if (mixproduct < -EPS)
//     {
//         vtxName = vertexName(f0M, f2M, f1M);
//         inEdgeName = edgeName(f1M, f0M);
//         outEdgeName = edgeName(f2M, f0M);
//     }
//     else
//     {
//         // std::cerr << "get_intVertex : edge(f1, f2) is on SP" << std::endl;
//         continue;
//     }

// }

// void segments::get_intVertex(const Solid *S, const Eigen::Vector4d &SP)
// {
//     std::vector<onFaceVertex*> onFaceVertexVector;
//     faceSolidMap f1{0, 0};
//     for (Edge* e = S->sedges; e != nullptr; e = e->nexte)
// 	{
// 		Vertex* v1 = e->he1->vtx;
// 		Vertex* v2 = e->he2->vtx;
// 		double d1 = dist(v1->vcoord, SP);
// 		double d2 = dist(v2->vcoord, SP);
// 		int s1 = compare(d1, 0.0, EPS);
// 		int s2 = compare(d2, 0.0, EPS);
//         faceSolidMap f2 {e->he1->wloop->lface->faceno, S->solidno};
//         faceSolidMap f3 {e->he2->wloop->lface->faceno, S->solidno};
// 		if ((s1 == -1 && s2 == 1) || (s1 == 1 && s2 == -1))
// 		{
// 			double t = d1 / (d1 - d2);
// 			Eigen::Vector3d point = v1->vcoord + t * (v2->vcoord - v1->vcoord);
// 			HalfEdge* he = nullptr;
// 			lmev(e->he1, (he = e->he2->nxt), ++maxv, point[0], point[1], point[2]);
//             Face *f2 = e->he1->wloop->lface;
//             Face *f3 = e->he2->wloop->lface;
// 			soov.insert(he->prv->vtx);
// 		}
// 		else
// 		{
// 			if (s1 == 0) soov.insert(v1);
// 			if (s2 == 0) soov.insert(v2);
// 		}
// 	}
//     // std::vector<intVertex*> intVectexVector;
//     faceSolidMap f0M{0, 0};
//     Edge *e = s->sedges;
//     while (e)
//     {
//         Face *f1 = e->he1->wloop->lface;
//         Face *f2 = e->he2->wloop->lface;
//         Eigen::Vector4d f1Equ = ((plane*)f1->surf)->equation;
//         Eigen::Vector4d f2Equ = ((plane*)f2->surf)->equation;
//         Eigen::Vector3d intPos(0.0, 0.0, 0.0);
//         bool flag = three_plane_intersect(SP, f1Equ, f2Equ, intPos);
//         if (flag == false)
//             continue;

//         if (f1 == f2)
//         {
//             // 暂时不支持一个边的相邻两面平行
//             // std::cerr << "get_intVertex : f1 == f2" << std::endl;
//             continue;;
//         }
//         Eigen::Vector3d face1Normal = f1->get_normal();
//         Eigen::Vector3d face2Normal = f2->get_normal();
//         Eigen::Vector3d face1face2Cross = face1Normal.cross(face2Normal);
//         face1face2Cross.normalize();
//         if (face1face2Cross.squaredNorm() < EPS * EPS)
//         {
//             // 暂时不支持边的相邻两个面平行(或者根本不用支持)
//             // std::cerr << "get_intVertex : f1 is parallel with f2" << std::endl;
//             continue;
//         }
//         Eigen::Vector3d spNormal{SP[0], SP[1], SP[2]};
//         double mixproduct = spNormal.dot(face1face2Cross);
//         faceSolidMap f1M{f1->faceno, s->solidno};
//         faceSolidMap f2M{f2->faceno, s->solidno};
//         edgeName inEdgeName;
//         edgeName outEdgeName;
//         vertexName vtxName;
//         if (mixproduct > EPS)
//         {
//             vtxName = vertexName(f0M, f1M, f2M);
//             inEdgeName = edgeName(f2M, f0M);
//             outEdgeName = edgeName(f1M, f0M);
//         }
//         else if (mixproduct < -EPS)
//         {
//             vtxName = vertexName(f0M, f2M, f1M);
//             inEdgeName = edgeName(f1M, f0M);
//             outEdgeName = edgeName(f2M, f0M);
//         }
//         else
//         {
//             // std::cerr << "get_intVertex : edge(f1, f2) is on SP" << std::endl;
//             continue;
//         }
//         set_intVertex_intEdge_relation(vtxName, inEdgeName, true, intPos);
//         set_intVertex_intEdge_relation(vtxName, outEdgeName, false, intPos);
//         e = e->nexte;
//     }
// }

// void segments::make_clusters()
// {
//     auto end = intVertexMap.end();
//     for (auto it = intVertexMap.begin(); it != end; ++it)
//     {
//         intVertex *intV = it->second;
//         Eigen::Vector3d itp = intV->pos;
//         auto endc = vtxClust.end();
//         bool flag = false;
//         for (auto itc = vtxClust.begin(); itc != endc; ++itc)
//         {
//             //算不知道什么均值，可能有bug
//             Eigen::Vector3d itcp = (*itc)->pos;
//             if ((itcp - itp).squaredNorm() < EPS * EPS)
//             {
//                 (*itc)->intVtxSet.push_back(it->second);
//                 (*itc)->pos = (itcp + itp) / 2.0;
//                 flag = true;
//                 VtxVtxclusRel.emplace(it->second, *itc); 
//                 break;
//             }
//         }
//         if (flag == false)
//         {
//             intVetexClusters *intvc = new intVetexClusters{itp};
//             intvc->intVtxSet.push_back(intV);
//         }

//     }
// }
