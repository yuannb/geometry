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

    for (auto it = m_node_data.begin(); it != m_node_data.end(); ++it)
    {
        delete *it;
    }
}

bool NodeGraph::calculate()
{
    for (auto node_id : m_nodes_sort)
    {
        if (m_nodes.contains(node_id))
        {
			Node* node = m_nodes[node_id];
			node->collect_data();
			node->calculate();
        }
    }
    return true;
}


