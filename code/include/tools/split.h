/*
对于流形. 在任意一点同胚于平面，因此局部上总是可以化成下图的形状

                A           B          A
                  x         x         x
                   x        x        x
                    x       x       x
                     x      x      x
                      x     x     x
                       x    x    x
                        x   x   x
                         x  x  x
                          x x x
                           xxx
         xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx B
                           xxx
                          x x x
                         x  x  x
                        x   x   x
                       x    x    x
                      x     x     x
                     x      x      x
                    x       x       x
                   x        x        x
                  x         x         x
                 x          x          x
                x           x           x
               x            x            x 
              A             A            B

*/
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "params.h"
#include "unordered_map"
#include "map"
#include "unordered_set"

#define ABOVE 1
#define BELOW -1
#define ON 0

struct hefrel
{
	HalfEdge* sector;
	int cl;
};

// struct faceSolidMap
// {
//     int faceNumber;
//     int solidNumber;
// };

typedef std::pair<Id, Id> faceSolidMap;

typedef std::pair<faceSolidMap, faceSolidMap> edgeName;
typedef std::tuple<faceSolidMap, faceSolidMap, faceSolidMap> vertexName;

struct intVertex;

struct intEdge
{
    edgeName intEdgeName;
    std::vector<intVertex*> stl;
    std::vector<intVertex*> edl;
    intEdge() = default;
    intEdge(edgeName en): intEdgeName(en) {}
};

struct intVertex
{
    vertexName intVtxName;
    Eigen::Vector3d pos;
    Id vertexno;
    std::vector<intEdge*> inEdge;
    std::vector<intEdge*> outEdge;
    intVertex(vertexName xintVtxName, Eigen::Vector3d xpos, Id xvertexno): intVtxName(xintVtxName)
                    ,pos(xpos), vertexno(xvertexno) {};
    // int edgeNo;
    // int faceNo;
};

// struct intEdgeClusters;

// struct intVetexClusters
// {
//     Eigen::Vector3d pos;
//     std::vector<intVertex*> intVtxSet;
//     std::vector<intEdgeClusters*> intEdgeSet;
//     std::vector<intEdgeClusters*> outEdgeSet;
// };

// struct intEdgeClusters
// {
//     edgeName intEdgeName;
//     std::vector<intVetexClusters*> stlc;
//     std::vector<intVetexClusters*> edlc;
// };



// struct onFaceVertex
// {
//     intVertex* intInfo;
//     vertex *v;
// };

double dist(const Eigen::Vector3d& v, const Eigen::Vector4d& vec);

// static bool isin(Edge *e, Edge* egelist);

// static bool isin(Vertex *v, Vertex* vertexlist);
int checkwideness(HalfEdge* he);
inline void bisector(HalfEdge* he, Eigen::Vector3d& bisect);
void cleanup(Solid* s);

void classify(Solid* S, Solid* Above, Solid* Below);

void movefac(Face* f, Solid* s);

int neighbor(HalfEdge* h1, HalfEdge* h2);

void cut(HalfEdge* he);

void splitfinish(Solid* S, Solid** Above, Solid** Below);

HalfEdge* canjoine(HalfEdge* he, std::unordered_set<Edge*> nulledgs);

void join(HalfEdge* h1, HalfEdge* h2);

bool splitconnect(std::map<intVertex*, std::vector<intVertex*>> wire, Solid *S);

bool sortintvertex(intEdge* intEdg, Eigen::Vector3d dir);

std::map<intVertex*, std::vector<intVertex*>> makeconnectgraphic(std::map<vertexName, intVertex*> &intVertexMap, 
                                                                std::map<edgeName, intEdge*> &intEdgeMap, Solid *S, Eigen::Vector3d normal);
void setintVertexintEdgerelation(vertexName vtxName, edgeName edgName, bool inout, Eigen::Vector3d &point, Id vertexno,
									std::map<vertexName, intVertex*> &intVertexMap, std::map<edgeName, intEdge*> &intEdgeMap);
bool createintvertex(HalfEdge *he, std::map<vertexName, intVertex*> &intVertexMap, std::map<edgeName, intEdge*> &intEdgeMap, Eigen::Vector4d SP);
int findfirstwidesector(std::vector<hefrel> &nbr);
std::unordered_set<Vertex*> splitgenerate(Solid* S, Eigen::Vector4d& SP);

std::vector<hefrel> getneighborhood(const Vertex* v, const Eigen::Vector4d& SP);

void reclassifyonsectors(const Eigen::Vector4d& SP, std::vector<hefrel> &nbr);

void reclassifyonedges(std::vector<hefrel> &nbr);

bool splitclassify(Eigen::Vector4d SP, std::unordered_set<Vertex*> &soov, std::map<vertexName, intVertex*> &intVertexMap, std::map<edgeName, intEdge*> &intEdgeMap);

bool insertnulledges(std::vector<hefrel> &nbr, std::map<vertexName, intVertex*> &intVertexMap, std::map<edgeName, intEdge*> &intEdgeMap, Eigen::Vector4d SP);

// class segments
// {
// private:
//     // std::vector<intVertex*> intVs;
//     // std::vector<intEdge*> intEdgs;
//     std::unordered_map<edgeName, intEdge*> intEdgeMap;
//     std::unordered_map<vertexName, intVertex*> intVertexMap;
//     std::unordered_map<intVertex*, intVetexClusters*> VtxVtxclusRel;
//     std::vector<intVetexClusters*> vtxClust;
//     // std::vector<intEdgeClusters*> edgeClust;
// public:
//     segments(/* args */) = default;
//     ~segments();
//     //inout: edge是vertex的入边(true)或者出边(false)
//     void set_intVertex_intEdge_relation(vertexName vtxName, edgeName edgName, bool inout, Eigen::Vector3d &point);
//     onFaceVertex* get_on_face_vertex(const Eigen::Vector4d &SP, Face *f2, Face *f3);
//     void get_intVertex(const Solid *s, const Eigen::Vector4d &SP);

//     void make_clusters();
// };

void split(Solid* S, Eigen::Vector4d &SP, Solid** Above, Solid** Below);


