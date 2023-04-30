#pragma once

#include "memory"

typedef int                Id;
typedef class  solid       Solid;
typedef class  face        Face;
typedef class  loop        Loop;
typedef class  halfedge    HalfEdge;
typedef class  vertex      Vertex;
typedef class  edge        Edge;
typedef union  nodes       Node;

// return values and misc constants
# define    ERROR       -1
# define    SUCCESS     -2
# define    PI          3.1415926535589793

// node type parameters for memory allocation routines
# define    SOLID       0
# define    FACE        1
# define    LOOP        2
# define    HALFEDGE    3
# define    EDGE        4
# define    VERTEX      5

// orientations
# define    PLUS    0
# define    MINUS   1

//tolerance
#define  EPS        1e-4      // tolerance for geometric tests
#define  BIGEPS     1e-8      // a more permissive tolerance 

// global variables
extern  std::weak_ptr<Solid>   firsts;    // head of the list of all solids
extern  Id      maxs;       // largest solid no
extern  Id      maxf;       // largest face no
extern  Id      maxv;       // largest vertes no
