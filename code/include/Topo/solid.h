#pragma once
#include "globalSymbol.h"

class solid
{
public:
    solid();
    //have to call it before deconstruct
    bool RemoveListFromSolid();
    Face *getface(Id faceno);
    Vertex *getvertex(Id vertexno);
    ~solid(){ }
public:
    Id      solidno;  //solid identifier
    Face    *sfaces;  //pointer to list of faces
    Edge    *sedges;  //pointer to list of edges
    Vertex  *svertes; //pointer to list of vertex
    Solid   *nexts;   //pointer to next solid
    Solid   *prevs;   //pointer to previous solid
};