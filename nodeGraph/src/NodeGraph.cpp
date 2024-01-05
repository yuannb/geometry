#include "NodeGraph.hpp"


NodeGraph::~NodeGraph() 
{
    for (auto it = m_nodes.begin(); it != m_nodes.end(); ++it)
    {
        delete it->second;
    }

    for (auto it = m_curves.begin(); it != m_curves.end(); ++it)
    {
        delete it->second;
    }

    for (auto it = m_curves_disc.begin(); it != m_curves_disc.end(); ++it)
    {
        delete it->second;
    }

    // std::map<std::string, tnurbs::curve3d*> m_curves;
    // std::map<std::string, tnurbs::mesh<1>*> m_curves_disc;
}

bool NodeGraph::calculate()
{
    for (auto node_id : m_nodes_sort)
    {
        Node *node = m_nodes[node_id];
        node->collect_data();
        node->calculate();
    }
    return true;
}


