#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <memory>
#include <fstream>
#include <algorithm>
#include <iosfwd>

class surface;
class face : public std::enable_shared_from_this<face>
{
public:
    face();
    
    //have to call it before deconstruct
    bool RemoveListFromSolid(std::shared_ptr<Solid> s);
    Eigen::Vector3d get_normal();
    void addlist(std::shared_ptr<Loop> l);
    
    ~face();
public:
    Id      faceno;  // face identifier;
    std::weak_ptr<Solid> fsolid; // back pointer to solid
    std::shared_ptr<Loop> flout;  // pointer to outer loop
    std::shared_ptr<Loop> floops; // pointer to list of loops
    std::shared_ptr<surface> surf; //face equation
    std::shared_ptr<Face> nextf;  // pointer to next face
    std::weak_ptr<Face> prevf;  // pointer to previous face
};
