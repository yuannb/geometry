#pragma once
#include "globalSymbol.h"
#include <memory>
    #include <fstream>
#include <algorithm>
#include <iosfwd>

class solid : public std::enable_shared_from_this<solid>
{
public:
    solid();
    //have to call it before deconstruct
    bool RemoveListFromSolid();
    bool addlist();
    std::shared_ptr<Face> getface(Id faceno);
    std::shared_ptr<Vertex> getvertex(Id vertexno);
    ~solid(){ }
    void addlist(std::shared_ptr<Edge> edg);
    void addlist(std::shared_ptr<Face> fac);
    void addlist(std::shared_ptr<Vertex> vtx);

public:
    Id      solidno;  //solid identifier
    std::shared_ptr<Face> sfaces; //pointer to list of faces
    // Face    *sfaces;  //pointer to list of faces
    std::shared_ptr<Edge> sedges; //pointer to list of edges
    // Edge    *sedges;  //pointer to list of edges
    std::shared_ptr<Vertex> svertes; //pointer to list of vertex
    // Vertex  *svertes; //pointer to list of vertex
    std::shared_ptr<Solid> nexts; //pointer to next solid
    // Solid   *nexts;   //pointer to next solid
    std::weak_ptr<Solid> prevs;   //pointer to previous solid
    // Solid   *prevs;   //pointer to previous solid
};