#pragma once
#include "json.hpp"
#include "NodeGraph.hpp"

bool CreateNodeGraph(const nlohmann::ordered_json& js, NodeGraph& ng);


