#pragma once
#include "euler.h"
#include "delaunay.h"
#include "memory"
#include "solid.h"
#include "params.h"


void arc(Id s, Id f, Id v, double cx, double cy, double rad, double h, double phi1, double phi2, int n);

std::shared_ptr<Solid> circle(Id sn, double cx, double cy, double rad, double h, int n);

void sweep(std::shared_ptr<Face> fac, double dx, double dy, double dz);

std::shared_ptr<Solid> block(Id sn, double dx, double dy, double dz);

std::shared_ptr<Solid>  cyl(Id sn, double rad, double h, int n);

std::shared_ptr<Solid> rsweep(std::shared_ptr<Solid> s, int nfaces);

std::shared_ptr<Solid> torus(Id sn, double r1, double r2, int nf1, int nf2);
