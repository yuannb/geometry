#pragma once
#include "../Geo/plane.h"
#include "../Topo/globalSymbol.h"

extern HalfEdge *hithe;
extern Vertex *hitvetex;

int faceeq(Loop *l, Eigen::Vector4d eq);

//a == b : return 0; a < b : return -1; a > b : return 1
int compare(double a, double b, double tol);

int contvv(Vertex *v1, Vertex *v2);

//evaluate project point
int intrev(Eigen::Vector3d &v1, Eigen::Vector3d &v2, Eigen::Vector3d &v3, double *t);

int intrev(Vertex *v1, Vertex *v2, Vertex *v3, double *t);

int intrev(Edge *edg, Vertex *v, double *t);

int contev(Vertex *v1, Vertex *v2, Vertex *v3);

int bndrlv(Loop *l, Vertex *v);

int int2ee(Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, int drop, double *t1, double *t2);

int contlv(Loop *l, Vertex *v, int drop);
