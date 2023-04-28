#pragma once
#include "geodef.h"
#include "gwballocation.h"
#include <iostream>
#include <cmath>
#include "query.h"
// #include <cstdlib>
HalfEdge *addhe(Edge *e, Vertex *v, HalfEdge *where, int sign);

HalfEdge *delhe(HalfEdge *he);

Solid *getsolid(Id sn);

Face    *fface(Solid *s, Id fn);

HalfEdge *fhe(Face *f, Id vn1, Id vn2);

HalfEdge *fhe(Face *f, Id vn1);

Vertex *getvertex(Solid *s, Id vno);

Solid *mvfs(Id s, Id f, Id v, double x, double y, double z);

void lmev(HalfEdge *he1, HalfEdge *he2, Id v, double x, double y, double z);
void lmev2(HalfEdge *he1, HalfEdge *he2, Id v, double x, double y, double z);

Face *lmef(HalfEdge *he1, HalfEdge *he2, Id f);

void lkef(HalfEdge *he1, HalfEdge *he2);

int kef(Id s, Id f, Id v1, Id v2);

void lkemr(HalfEdge *h1, HalfEdge *h2);

int kemr(Id s, Id f, Id v1, Id v2);

void lmekr(HalfEdge *h1, HalfEdge *h2);

int mekr(Id s, Id f, Id v1, Id v2, Id v3, Id v4);

int smekr(Id s, Id f, Id v1, Id v3);

void lkfmrh(Face *fac1, Face *fac2);

int kfmrh(Id s, Id f1, Id f2);

Face* lmfkrh(Loop *l, Id f);

int mfkrh(Id s, Id f1, Id v1, Id v2, Id f2);

void lringmv(Loop *l, Face *toface, int inout);

int ringmv(Solid *s, Id f1, Id f2, Id v1, Id v2, int inout);

int mev(Id s, Id f1, Id f2, Id v1, Id v2, Id v3, Id v4, double x, double y, double z);

int smev(Id s, Id f1, Id v1, Id v4, double x, double y, double z);

int smef(Id s, Id f1, Id v1, Id v3, Id f2);

void lkvfs(Solid *s);

void kvfs(Id s);

void lkev(HalfEdge *he1, HalfEdge *he2);

int kev(Id s, Id f, Id v1, Id v2);

void getmaxnames(Id sn);

void merge(Solid *s1, Solid *s2);

double distancetwovector(vector v1, vector v2);

bool match(HalfEdge *h1, HalfEdge *h2);

void loopglue(Face *f);

void glue(Solid *s1, Solid *s2, Face *f1, Face *f2);

