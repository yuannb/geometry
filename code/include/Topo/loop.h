#pragma once

#include <memory>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iosfwd>
#include "params.h"

class loop : public std::enable_shared_from_this<loop>
{
public:
    loop();
    //have to call it before deconstruct
    bool RemoveListFromFace(std::shared_ptr<Face> f);

    ~loop();
public:
    std::shared_ptr<HalfEdge> ledg; // ptr to ring of halfedges
    std::weak_ptr<Face> lface; // back pointer to face
    std::shared_ptr<Loop> nextl; // pointer to next loop
    std::weak_ptr<Loop> prevl; // pointer to previous loop
};
