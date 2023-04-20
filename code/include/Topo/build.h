#pragma once
#include "euler.h"
#include "delaunay.h"

void arc(Id s, Id f, Id v, double cx, double cy, double rad, double h, double phi1, double phi2, int n);

Solid *circle(Id sn, double cx, double cy, double rad, double h, int n);

void sweep(Face *fac, double dx, double dy, double dz);

Solid *block(Id sn, double dx, double dy, double dz);

Solid *cyl(Id sn, double rad, double h, int n);

Solid *rsweep(Solid *s, int nfaces);

Solid *torus(Id sn, double r1, double r2, int nf1, int nf2);
