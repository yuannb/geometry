#pragma once
#include "parameter.hpp"
#include "json.hpp"
class Node;
class AbstractNodeData;

struct serialize_utils
{
    static AbstractNodeData* to_nodedata(nlohmann::ordered_json node_data);

    static bool to_path(const std::string& s, char delimiter, Path& path);

    static bool to_intface(nlohmann::ordered_json intf_json, NodeInterface& intf);

};
