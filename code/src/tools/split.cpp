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

#define ABOVE 1
#define BELOW -1
#define ON 0

std::map<std::shared_ptr<Vertex>, std::vector<std::shared_ptr<Edge>>> nulledges;
std::unordered_set<std::shared_ptr<HalfEdge>> hasconnect;

std::vector<std::shared_ptr<Face>> sonf;

struct hefrel
{
    std::shared_ptr<HalfEdge> sector;
	int cl;
};

typedef std::pair<Id, Id> faceSolidMap;

typedef std::pair<faceSolidMap, faceSolidMap> edgeName;
typedef std::tuple<faceSolidMap, faceSolidMap, faceSolidMap> vertexName;

struct intVertex;

struct intEdge
{
    edgeName intEdgeName;
    std::vector<std::shared_ptr<intVertex>> stl;
    std::vector<std::shared_ptr<intVertex>> edl;
    intEdge() = default;
    intEdge(edgeName en): intEdgeName(en) {}
};

struct intVertex
{
    vertexName intVtxName;
    Eigen::Vector3d pos;
    Id vertexno;
    std::vector<std::shared_ptr<intEdge>> inEdge;
    std::vector<std::shared_ptr<intEdge>> outEdge;
    intVertex(vertexName xintVtxName, Eigen::Vector3d xpos, Id xvertexno): intVtxName(xintVtxName)
                    ,pos(xpos), vertexno(xvertexno) {};
    // int edgeNo;
    // int faceNo;
};

int neighbor(std::shared_ptr<HalfEdge> h1, std::shared_ptr<HalfEdge> h2)
{
	return (h1->wloop.lock()->lface.lock() == h2->wloop.lock()->lface.lock() && (
		(h1 == h1->edg->he1.lock() && h2 == h2->edg->he2.lock()) || (h1 == h1->edg->he2.lock() && h2 == h2->edg->he1.lock())
		));
}

void cut(std::shared_ptr<HalfEdge> he)
{
	std::shared_ptr<HalfEdge> edghe1 = he->edg->he1.lock();
	std::shared_ptr<HalfEdge> edghe2 = he->edg->he2.lock();
	if (edghe1->wloop.lock() == edghe2->wloop.lock())
	{
		sonf.push_back(he->wloop.lock()->lface.lock());
		lkemr(edghe1, edghe2);
	}
	else
	{
		lkef(edghe1, edghe2);
	}
}

std::shared_ptr<HalfEdge> canjoine(std::shared_ptr<HalfEdge> he, std::unordered_set<std::shared_ptr<Edge>> nulledgs)
{
	size_t cnt = hasconnect.count(he);
	if (cnt == 1)
		return nullptr;
	auto end = nulledgs.end();
	for (auto it = nulledgs.begin(); it != end; ++it)
	{
		if (neighbor(he, (*it)->he1.lock()))
		{
			hasconnect.insert((*it)->he1.lock());
			hasconnect.insert(he);
			return (*it)->he1.lock();
		}
		if (neighbor(he, (*it)->he2.lock()))
		{
			hasconnect.insert((*it)->he2.lock());
			hasconnect.insert(he);
			return (*it)->he2.lock();
		}
	}
	return nullptr;
}

void join(std::shared_ptr<HalfEdge> h1, std::shared_ptr<HalfEdge> h2)
{
	std::shared_ptr<Face> oldface = h1->wloop.lock()->lface.lock();
	std::shared_ptr<Face> newface = nullptr;
	if (h1->wloop.lock() == h2->wloop.lock())
	{
		if (h1->prv.lock()->prv.lock() != h2)
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
	double product = v1v2.dot(dir);
	if (product > 0)
		return true;
	return false;
}

bool sortintvertex(std::shared_ptr<intEdge> intEdg, Eigen::Vector3d dir)
{	
	std::sort(intEdg->stl.begin(), intEdg->stl.end(), [&](std::shared_ptr<intVertex> v1, std::shared_ptr<intVertex> v2) {return (v2->pos - v1->pos).dot(dir) > 0; });
	std::sort(intEdg->edl.begin(), intEdg->edl.end(), [&](std::shared_ptr<intVertex> v1, std::shared_ptr<intVertex> v2) {return (v2->pos - v1->pos).dot(dir) > 0; });

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
			if (!comparevertex(intEdg->stl[0].get(), intEdg->edl[0].get(), dir))
			{
				std::cerr << "sortintvertex : intEdg->stl[0] >  intEdg->edl[0]" << std::endl;
			}
		}
		else if (i == stcount - 1) //除掉了stcount == 1的情况
		{
			if (!comparevertex(intEdg->edl[i - 1].get(), intEdg->stl[i].get(), dir))
			{
				std::cerr << "sortintvertex : intEdg->edl[i - 1] > intEdg->stl[i]" << std::endl;
			}
		}
	}
	return true;
}

std::map<std::shared_ptr<intVertex>, std::vector<std::shared_ptr<intVertex>>> makeconnectgraphic(std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, 
												std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap, Solid *S, Eigen::Vector3d normal)
{
	std::map<std::shared_ptr<intVertex>, std::vector<std::shared_ptr<intVertex>>> ans;
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
				std::vector<std::shared_ptr<intVertex>> endpoint{it->second->edl[index]};
				ans.emplace(it->second->stl[index], endpoint);
			}
		}
	}
	return ans;
}

std::unordered_set<std::shared_ptr<Vertex>> splitgenerate(std::shared_ptr<Solid> S, Eigen::Vector4d& SP)
{
	//nvtx = 0;
    std::unordered_set<std::shared_ptr<Vertex>> soov;
	for (std::shared_ptr<Edge> e = S->sedges; e != nullptr; e = e->nexte)
	{
		std::shared_ptr<Vertex> v1 = e->he1.lock()->vtx;
		std::shared_ptr<Vertex> v2 = e->he2.lock()->vtx;
		double d1 = dist(v1->vcoord, SP);
		double d2 = dist(v2->vcoord, SP);
		int s1 = compare(d1, 0.0, EPS);
		int s2 = compare(d2, 0.0, EPS);
		if ((s1 == -1 && s2 == 1) || (s1 == 1 && s2 == -1))
		{
			double t = d1 / (d1 - d2);
			Eigen::Vector3d point = v1->vcoord + t * (v2->vcoord - v1->vcoord);
			std::shared_ptr<HalfEdge> he = e->he2.lock()->nxt;
			lmev(e->he1.lock(), he, ++maxv, point[0], point[1], point[2]);
			soov.insert(he->prv.lock()->vtx);
		}
		else
		{
			if (s1 == 0) soov.insert(v1);
			if (s2 == 0) soov.insert(v2);
		}
	}
    return soov;
}

int checkwideness(std::shared_ptr<HalfEdge> he)
{
	Eigen::Vector3d nor = he->mate()->wloop.lock()->lface.lock()->surf->get_normal();
	he = he->mate()->nxt;
	std::shared_ptr<HalfEdge> nxthe = he->nxt;
	Eigen::Vector3d v0 = he->vtx->vcoord;
	Eigen::Vector3d v1 = nxthe->vtx->vcoord;
	Eigen::Vector3d v2 = he->prv.lock()->vtx->vcoord;
	Eigen::Vector3d v0v1 = v1 - v0;
	Eigen::Vector3d v2v0 = v0 - v2;
	Eigen::Vector3d c = v2v0.cross(v0v1);
	double d = nor.dot(c);
	return d < -EPS ? 0 : 1;
}
inline void bisector(HalfEdge* he, Eigen::Vector3d& bisect)
{
	bisect = -(he->vtx->vcoord + he->prv.lock()->vtx->vcoord) / 2;
};

std::vector<hefrel> getneighborhood(const Vertex* v, const Eigen::Vector4d& SP)
{
	std::vector<hefrel> nbr;
	std::shared_ptr<HalfEdge> he = v->vedge.lock();
	std::shared_ptr<HalfEdge> first = he;
	// HalfEdge* he = v->vedge;
	do
	{
		double d = dist(he->nxt->vtx->vcoord, SP);
		int cl = compare(d, 0.0, EPS);
		hefrel rel{ he, cl };
		nbr.push_back(rel);
		if (checkwideness(he) != 1)
		{
			Eigen::Vector3d bisec;
			bisector(he.get(), bisec);
			d = dist(bisec, SP);
			int cl2 = compare(d, 0.0, EPS);
			hefrel rel2{ he, cl2 };
			nbr.push_back(rel2);
		}
	} while ((he = he->mate()->nxt) != first);
    return nbr;
}

void reclassifyonsectors(const Eigen::Vector4d& SP, std::vector<hefrel> &nbr)
{
	auto end = nbr.end();
	for (auto it = nbr.begin(); it != end; ++it)
	{
		Face* f = it->sector->wloop.lock()->lface.lock().get();
		Eigen::Vector4d eq = ((plane*)f->surf.get())->equation;
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

void setintVertexintEdgerelation(vertexName vtxName, edgeName edgName, bool inout, Eigen::Vector3d &point, Id vertexno,
		std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap)
{
    auto vtxIt = intVertexMap.find(vtxName);
    auto edgIt = intEdgeMap.find(edgName);
    std::shared_ptr<intVertex> intVtx = nullptr;
    std::shared_ptr<intEdge> intEdg = nullptr;
    if (vtxIt == intVertexMap.end())
    {
        intVtx = std::make_shared<intVertex> (vtxName, point, vertexno);
        intVertexMap.emplace(vtxName, intVtx);
    }
    else
    {
        intVtx = vtxIt->second;
    }

    if (edgIt == intEdgeMap.end())
    {
        intEdg = std::make_shared<intEdge> (edgName);
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


bool createintvertex(HalfEdge *he, std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, 
						std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap, Eigen::Vector4d SP)
{
	faceSolidMap f0No{0, 0};
	Face *f1 = he->wloop.lock()->lface.lock().get();
	Face *f2 = he->mate()->wloop.lock()->lface.lock().get();
	int soildN0 = f1->fsolid.lock()->solidno;
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
	Eigen::Vector3d testv = he->nxt->nxt->vtx->vcoord - he->vtx->vcoord;
	double mixproduct = testv.dot(f0Nor);
	if (std::abs(mixproduct) < EPS)
	{
		//不应该出现这种情况(除非是非流形或者精度问题)
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

bool insertnulledges(std::vector<hefrel> &nbr, std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, 
						std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap, Eigen::Vector4d SP)
{
    int i = 0;
	int nnbr = (int)nbr.size();
	while (!(nbr[i].cl == BELOW && nbr[(i + 1) % nnbr].cl == ABOVE))
	{
		if (++i == nnbr) return true;
	}
	int start = i;
	//book may have bug
	std::vector<std::shared_ptr<Edge>> nulledgs;
	while (1)
	{
		std::shared_ptr<HalfEdge> head = nbr[(i + 1) % nnbr].sector;
		if  (nbr[i].sector == nbr[(i + 1) % nnbr].sector)
			head = nbr[(i + 2) % nnbr].sector;
		while (!(nbr[i].cl == ABOVE && nbr[(i + 1) % nnbr].cl == BELOW))
		{
			i = (i + 1) % nnbr;
		}
		std::shared_ptr<HalfEdge> tail = nbr[(i + 1) % nnbr].sector;
		if (nbr[i].sector == nbr[(i + 1) % nnbr].sector)
			tail = nbr[i].sector->mate()->nxt;
		lmev(head, tail, ++maxv, head->vtx->vcoord[0],
			head->vtx->vcoord[1],
			head->vtx->vcoord[2]);
		nulledgs.push_back(head->prv.lock()->edg);
		createintvertex(head->prv.lock().get(), intVertexMap, intEdgeMap, SP);
		while (!(nbr[i].cl == BELOW && nbr[(i + 1) % nnbr].cl == ABOVE))
		{
			i = (i + 1) % nnbr;
			if (i == start) 
			{
				if (!nulledgs.empty())
				{
					nulledges[head->prv.lock()->vtx] = nulledgs;
				}
				return true;
			}
		}
	}
	return true;
}

bool splitclassify(Eigen::Vector4d SP, std::unordered_set<std::shared_ptr<Vertex>> &soov, std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, 
						std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap)
{
    auto end = soov.end();
	for (auto it = soov.begin(); it != end; ++it)
	{
		std::vector<hefrel> nbr = getneighborhood(it->get(), SP);
		reclassifyonsectors(SP, nbr);
		reclassifyonedges(nbr);
        insertnulledges(nbr, intVertexMap, intEdgeMap, SP);
	}
	return true;
}

bool splitconnect(std::map<std::shared_ptr<intVertex>, std::vector<std::shared_ptr<intVertex>>> wire, std::shared_ptr<Solid> S)
{
	auto wireend = wire.end();
	for (auto wireit = wire.begin(); wireit != wireend; ++wireit)
	{
		std::shared_ptr<Vertex> v = S->getvertex(wireit->first->vertexno);
		std::vector<std::shared_ptr<Edge>> edgs = nulledges[v];
		auto edgsend = edgs.end();
		std::unordered_set<std::shared_ptr<Edge>> canjoineedge;
		auto intvetxend = wireit->second.end();
		for (auto intvetxeit = wireit->second.begin(); intvetxeit != intvetxend; ++intvetxeit)
		{
			std::shared_ptr<Vertex> vs = S->getvertex((*intvetxeit)->vertexno);
			canjoineedge.insert(nulledges[vs].begin(), nulledges[vs].end());
		}
		for (auto edgsit = edgs.begin(); edgsit != edgsend; ++edgsit)
		{
			std::shared_ptr<HalfEdge> h1 = canjoine((*edgsit)->he1.lock(), canjoineedge);
			if (h1)
			{
				join(h1, (*edgsit)->he1.lock());
				size_t cnt = hasconnect.count(h1->mate());
				if (cnt == 1) cut(h1);
				cnt = hasconnect.count((*edgsit)->he1.lock()->mate());
				if (cnt == 1) cut((*edgsit)->he1.lock());
			}
			std::shared_ptr<HalfEdge> h2 = canjoine((*edgsit)->he2.lock(), canjoineedge);
			if (h2)
			{
				join(h2, (*edgsit)->he2.lock());
				size_t cnt = hasconnect.count(h2->mate());
				if (cnt == 1) cut(h2);
				cnt = hasconnect.count((*edgsit)->he2.lock()->mate());
				if (cnt == 1) cut((*edgsit)->he2.lock());
			}
			if (h1 && h2)
			{
				std::cerr << "splitconnect : some unkown error occured" << std::endl;
			}
		}
	}
	return true;
	
}

void  movefac(std::shared_ptr<Face> f, std::shared_ptr<Solid> s)
{
	std::vector<std::shared_ptr<Face>> fvec{f};
	size_t index = 0;
	while (index != fvec.size())
	{
		fvec[index]->RemoveListFromSolid(fvec[index]->fsolid.lock());
		s->addlist(fvec[index]);
		std::shared_ptr<Loop> l = fvec[index]->floops;
		while (l)
		{
			std::shared_ptr<HalfEdge> he = l->ledg;
			std::shared_ptr<HalfEdge> first = he;
		    do
			{
	 			std::shared_ptr<Face> f2 = he->mate()->wloop.lock()->lface.lock();
					if (f2->fsolid.lock() != s)
					{
						// f->RemoveListFromSolid(f->fsolid.lock());
						// s->addlist(f);
						fvec.push_back(f2);
					}
	 		} while ((he = he->nxt) != first);
			l = l->nextl;
		}
		++index;
	}
	
}

void classify(std::shared_ptr<Solid> S, std::shared_ptr<Solid> &Above, std::shared_ptr<Solid> &Below)
{
	int nfac = sonf.size() / 2;
	for (int i = 0; i != nfac; ++i)
	{
		//暂时不支持内环，稍微修改一下即可
		movefac(sonf[i], Above);
	}
	for (int i = 0; i != nfac; ++i)
	{
		movefac(sonf[i + nfac], Below);
	}

}

void cleanup(std::shared_ptr<Solid> s)
{
	std::shared_ptr<Face> f = s->sfaces;
	std::unordered_set<std::shared_ptr<Edge>> edgeset;
	std::unordered_set<std::shared_ptr<Vertex>> vtxset;
	while (f)
	{
		std::shared_ptr<Loop> l = f->floops;
		while (l)
		{
			std::shared_ptr<HalfEdge> he = l->ledg;
			do
			{
				std::shared_ptr<Edge> e = he->edg;
				int cnt = edgeset.count(e);
				if (cnt == 0)
				{
					s->addlist(e);
					edgeset.insert(e);
				}
				std::shared_ptr<Vertex> v = he->vtx;
				cnt = vtxset.count(v);
				if (cnt == 0)
				{
					s->addlist(v);
					vtxset.insert(v);
				}
			} while ((he = he->nxt) != l->ledg);
			l = l->nextl;
		}
		f = f->nextf;
	}
}

void splitfinish(std::shared_ptr<Solid> S, std::shared_ptr<Solid> &Above, std::shared_ptr<Solid> &Below)
{
	int nfac = (int)sonf.size();
	for (int i = 0; i != nfac; ++i)
	{
		std::shared_ptr<Loop> l = sonf[i]->floops;
		if (l == sonf[i]->flout)
			l = sonf[i]->floops->nextl;
		std::shared_ptr<Face> f = lmfkrh(l, ++maxf);
		sonf.push_back(f);
	}
	Above = std::make_shared<Solid>();
	Above->addlist();
	Below = std::make_shared<Solid>();
	Below->addlist();
	classify(S, Above, Below);
	cleanup(Above);
	cleanup(Below);
}

void split(std::shared_ptr<Solid> S, Eigen::Vector4d &SP, std::shared_ptr<Solid> &Above, std::shared_ptr<Solid> &Below)
{
    std::map<vertexName, std::shared_ptr<intVertex>> intVertexMap;
    std::map<edgeName, std::shared_ptr<intEdge>> intEdgeMap;
	for (std::shared_ptr<Face> f = S->sfaces; f != nullptr; f = f->nextf)
	{
		Eigen::Vector4d eq;
		faceeq(f->flout.get(), eq);
		std::shared_ptr<plane> p = std::make_shared<plane>(eq);
		f->surf = std::static_pointer_cast<surface>(p);
	}
	getmaxnames(S->solidno);
	std::unordered_set<std::shared_ptr<Vertex>> soov = splitgenerate(S, SP);
    splitclassify(SP, soov, intVertexMap, intEdgeMap);
	Eigen::Vector3d normal(SP[0], SP[1], SP[2]);
	std::map<std::shared_ptr<intVertex>, std::vector<std::shared_ptr<intVertex>>> wire = makeconnectgraphic(intVertexMap, intEdgeMap, S.get(), normal);
	splitconnect(wire, S);
	splitfinish(S, Above, Below);
}