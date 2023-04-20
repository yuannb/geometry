#pragma once
#include "globalSymbol.h"

// return values and misc constants
# define    ERROR       -1
# define    SUCCESS     -2
# define    NIL         0
# define    PI          3.1415926535589793

// node type parameters for memory allocation routines
# define    SOLID       0
# define    FACE        1
# define    LOOP        2
# define    HALFEDGE    3
# define    EDGE        4
# define    VERTEX      5

// coordinate plane names
# define    XPLANE      0
# define    YPLANE      1
# define    ZPLANE      2

// orientations
# define    PLUS    0
# define    MINUS   1

//tolerance
#define  EPS        1e-4      // tolerance for geometric tests
#define  BIGEPS     1e-8      // a more permissive tolerance 

// macros
# define    mate(he)    (((he) == (he)->edg->he1) ? \
                        (he)->edg->he2 : (he)->edg->he1)
// # define    max(x, y)   ((x) > (y) ? (x) : (y))
// # define     abs(x)      ((x) > 0.0 ? (x) : (-y))

// global variables
extern  Solid           *firsts;    // head of the list of all solids
extern  Id      maxs;       // largest solid no
extern  Id      maxf;       // largest face no
extern  Id      maxv;       // largest vertes no
