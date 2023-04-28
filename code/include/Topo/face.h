#pragma once
#include "globalSymbol.h"
#include <Eigen/Core>
#include <Eigen/Dense>

class surface;
class face
{
public:
    face(Solid *s);
    
    //have to call it before deconstruct
    bool RemoveListFromSolid(Solid *s);
    Eigen::Vector3d get_normal();
    
    ~face() { };
public:
    Id      faceno;  // face identifier;
    Solid   *fsolid; // back pointer to solid
    Loop    *flout;  // pointer to outer loop
    Loop    *floops; // pointer to list of loops
    surface *surf;     // face equation
    Face    *nextf;  // pointer to next face
    Face    *prevf;  // pointer to previous face
};
