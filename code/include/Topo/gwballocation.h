#pragma once
#include "geodef.h"
#include <cstdlib>

void addlist(int what, Node *which, Node *where);
// Node *gnew(int what, Node *where);
void removelist(int what, Node *which, Node *where);
// void Free(Solid *dsolid);
// void Free(Face *face);
// void Free(Loop *dl);
// void Free(Edge *edg);
// void Free(Vertex *vtx);
// void Free(HalfEdge *he);
// HalfEdge *findOtherHalfEdge(HalfEdge *he);

// // TODO : free storage
// void gdel(int what, Node *which, Node *where);