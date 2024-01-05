#pragma once
#include "node.hpp"
#include <map>
#include <string>

class NodeGraph
{
public:
    std::string m_id;
    std::map<std::string, Node*>  m_nodes;
    

    std::map<std::string, tnurbs::curve3d*> m_curves;
    std::map<std::string, tnurbs::mesh<1>*> m_curves_disc;

    // std::vector<std::unique_ptr<tnurbs::curve<double, 3>>>  m_curves;

    std::vector<std::string> m_nodes_sort;

    NodeGraph(const std::string &id): m_id(id) { }
    ~NodeGraph(); 
    // {
    //     auto it = m_nodes.begin();
    //     auto end = m_nodes.end();
    //     for (; it != end; ++it)
    //     {
    //         delete it->second;
    //     }
    // }

    bool calculate();
    // {
    //     auto end = m_nodes.end();
    //     for (auto it = m_nodes.begin(); it != end; ++it)
    //     {
    //         it->second->collect_data();
    //         it->second->calculate();
    //     }
    //     return true;
    // }

};



