#pragma once
#include "geodef.h"
#include "gwballocation.h"
#include <iostream>
#include <cmath>
#include "query.h"
// #include <cstdlib>
std::shared_ptr<HalfEdge> addhe(std::shared_ptr<Edge> e, std::shared_ptr<Vertex> v, std::shared_ptr<HalfEdge> where, int sign);

std::shared_ptr<HalfEdge> delhe(std::shared_ptr<HalfEdge> he);

std::shared_ptr<Solid> getsolid(Id sn);

std::shared_ptr<Face>  fface(Solid *s, Id fn);

std::shared_ptr<HalfEdge> fhe(Face *f, Id vn1, Id vn2);

std::shared_ptr<HalfEdge> fhe(Face *f, Id vn1);

std::shared_ptr<Vertex> getvertex(Solid *s, Id vno);

std::shared_ptr<Solid> mvfs(Id s, Id f, Id v, double x, double y, double z);

void lmev(std::shared_ptr<HalfEdge> he1, std::shared_ptr<HalfEdge> he2, Id v, double x, double y, double z);

std::shared_ptr<Face> lmef(std::shared_ptr<HalfEdge> he1, std::shared_ptr<HalfEdge> he2, Id f);

void lkef(std::shared_ptr<HalfEdge> he1, std::shared_ptr<HalfEdge> he2);

int kef(Id s, Id f, Id v1, Id v2);

void lkemr(std::shared_ptr<HalfEdge> h1, std::shared_ptr<HalfEdge> h2);

int kemr(Id s, Id f, Id v1, Id v2);

void lmekr(std::shared_ptr<HalfEdge> h1, std::shared_ptr<HalfEdge> h2);

int mekr(Id s, Id f, Id v1, Id v2, Id v3, Id v4);

int smekr(Id s, Id f, Id v1, Id v3);

void lkfmrh(std::shared_ptr<Face> fac1, std::shared_ptr<Face> fac2);

int kfmrh(Id s, Id f1, Id f2);

std::shared_ptr<Face> lmfkrh(std::shared_ptr<Loop> l, Id f);

int mfkrh(Id s, Id f1, Id v1, Id v2, Id f2);

void lringmv(std::shared_ptr<Loop> l, std::shared_ptr<Face> toface, int inout);

int ringmv(std::shared_ptr<Solid> s, Id f1, Id f2, Id v1, Id v2, int inout);

int mev(Id s, Id f1, Id f2, Id v1, Id v2, Id v3, Id v4, double x, double y, double z);

int smev(Id s, Id f1, Id v1, Id v4, double x, double y, double z);

int smef(Id s, Id f1, Id v1, Id v3, Id f2);

void lkvfs(std::shared_ptr<Solid> s);

void kvfs(Id s);

void lkev(std::shared_ptr<HalfEdge> h1, std::shared_ptr<HalfEdge> h2);

int kev(Id s, Id f, Id v1, Id v2);

void getmaxnames(Id sn);

void merge(std::shared_ptr<Solid> s1, std::shared_ptr<Solid> s2);

// double distancetwovector(vector v1, vector v2);

bool match(HalfEdge *h1, HalfEdge *h2);

void loopglue(std::shared_ptr<Face> f);

void glue(std::shared_ptr<Solid> s1, std::shared_ptr<Solid> s2, std::shared_ptr<Face> f1, std::shared_ptr<Face> f2);

