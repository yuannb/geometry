#pragma once
#include "../Geo/plane.h"
#include "../Topo/globalSymbol.h"
#include "params.h"
#include "solid.h"
#include <set>
#include <vector>
#include "euler.h"
#include "memory"


bool three_plane_intersect(const Eigen::Vector4d &p1, const Eigen::Vector4d &p2, const Eigen::Vector4d &p3, Eigen::Vector3d &intVtx);


int faceeq(Loop *l, Eigen::Vector4d &eq);

//a == b : return 0; a < b : return -1; a > b : return 1
int compare(double a, double b, double tol);

int contvv(Vertex *v1, Vertex *v2);

//evaluate project point
int intrev(Eigen::Vector3d &v1, Eigen::Vector3d &v2, Eigen::Vector3d &v3, double *t);

int intrev(Vertex *v1, Vertex *v2, Vertex *v3, double *t);

int intrev(Edge *edg, Vertex *v, double *t);

int contev(Vertex *v1, Vertex *v2, Vertex *v3);

int bndrlv(Loop* l, Vertex* v);

int int2ee(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4, double* t1, double* t2);

int contlv(Loop *l, Vertex *v, int drop);

double svolume(Solid* s);

double larea(Loop* l);

int contfv(Face* f, Vertex* v);

int contfp(Face* f, double x, double y, double z);

void laringmv(std::shared_ptr<Face> f1, std::shared_ptr<Face> f2);