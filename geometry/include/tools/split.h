/*
对于流形. 在任意一点同胚于平面，因此局部上在某一点总是可以画成下图的形状

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

// struct hefrel
// {
//     std::shared_ptr<HalfEdge> sector;
// 	int cl;
// };

// typedef std::pair<Id, Id> faceSolidMap;

// typedef std::pair<faceSolidMap, faceSolidMap> edgeName;
// typedef std::tuple<faceSolidMap, faceSolidMap, faceSolidMap> vertexName;

// struct intVertex;

// struct intEdge
// {
//     edgeName intEdgeName;
//     std::vector<std::shared_ptr<intVertex>> stl;
//     std::vector<std::shared_ptr<intVertex>> edl;
//     intEdge() = default;
//     intEdge(edgeName en): intEdgeName(en) {}
// };

// struct intVertex
// {
//     vertexName intVtxName;
//     Eigen::Vector3d pos;
//     Id vertexno;
//     std::vector<std::shared_ptr<intEdge>> inEdge;
//     std::vector<std::shared_ptr<intEdge>> outEdge;
//     intVertex(vertexName xintVtxName, Eigen::Vector3d xpos, Id xvertexno): intVtxName(xintVtxName)
//                     ,pos(xpos), vertexno(xvertexno) {};
//     // int edgeNo;
//     // int faceNo;
// };


// double dist(const Eigen::Vector3d& v, const Eigen::Vector4d& vec);

// int checkwideness(std::shared_ptr<HalfEdge> he);
// inline void bisector(HalfEdge* he, Eigen::Vector3d& bisect);
// void cleanup(std::shared_ptr<Solid> s);

// void classify(std::shared_ptr<Solid> S, std::shared_ptr<Solid> &Above, std::shared_ptr<Solid> &Below);

// void  movefac(std::shared_ptr<Face> f, std::shared_ptr<Solid> s);

// int neighbor(std::shared_ptr<HalfEdge> h1, std::shared_ptr<HalfEdge> h2);

// void cut(std::shared_ptr<HalfEdge> he);

// void splitfinish(std::shared_ptr<Solid> S, std::shared_ptr<Solid> &Above, std::shared_ptr<Solid> &Below);

// std::shared_ptr<HalfEdge> canjoine(std::shared_ptr<HalfEdge> he, std::unordered_set<std::shared_ptr<Edge>> nulledgs);

// void join(std::shared_ptr<HalfEdge> h1, std::shared_ptr<HalfEdge> h2);

// bool splitconnect(std::map<std::shared_ptr<intVertex>, std::vector<std::shared_ptr<intVertex>>> wire, std::shared_ptr<Solid> S);

// bool sortintvertex(std::shared_ptr<intEdge> intEdg, Eigen::Vector3d dir);

// std::map<std::shared_ptr<intVertex>, std::vector<std::shared_ptr<intVertex>>> makeconnectgraphic(std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, 
// 												std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap, Solid *S, Eigen::Vector3d normal);
// void setintVertexintEdgerelation(vertexName vtxName, edgeName edgName, bool inout, Eigen::Vector3d &point, Id vertexno,
// 		std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap);

// bool createintvertex(HalfEdge *he, std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, 
// 						std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap, Eigen::Vector4d SP);

// std::unordered_set<std::shared_ptr<Vertex>> splitgenerate(std::shared_ptr<Solid> S, Eigen::Vector4d& SP);

// std::vector<hefrel> getneighborhood(const Vertex* v, const Eigen::Vector4d& SP);

// void reclassifyonsectors(const Eigen::Vector4d& SP, std::vector<hefrel> &nbr);

// void reclassifyonedges(std::vector<hefrel> &nbr);

// bool splitclassify(Eigen::Vector4d SP, std::unordered_set<std::shared_ptr<Vertex>> &soov, std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, 
// 						std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap);

// bool insertnulledges(std::vector<hefrel> &nbr, std::map<vertexName, std::shared_ptr<intVertex>> &intVertexMap, 
// 						std::map<edgeName, std::shared_ptr<intEdge>> &intEdgeMap, Eigen::Vector4d SP);

void split(std::shared_ptr<Solid> S, Eigen::Vector4d &SP, std::shared_ptr<Solid> &Above, std::shared_ptr<Solid> &Below);


