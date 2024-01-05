#pragma once
#include <json.hpp>
#include "NodeGraph.hpp"
using json = nlohmann::json;

bool CreateNodeGraphFromJson(const json &js, NodeGraph &ng);
