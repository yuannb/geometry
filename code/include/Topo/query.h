#pragma once
#include    "geodef.h"
# include <iostream>
#include "memory"
#include "params.h"

void listsolid(std::shared_ptr<Solid> s);

void listneighbors(std::shared_ptr<Vertex> v);