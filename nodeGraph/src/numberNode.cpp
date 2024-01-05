#include "numberNode.hpp"
#include "NodeGraph.hpp"

NumberNode::NumberNode(const std::string &id, double num)
{
    m_inputs.clear();
    m_outputs.resize(1);

    m_id = id;
    NodeInterface *interface = new NodeInterface();
    m_outputs[0] = interface;
    NumberData *data = new NumberData(num);
    std::vector<std::vector<unsigned>> path_data{ {0} };
    Path path(path_data);
    interface->m_data[path_data] = std::vector<AbstractNodeData*> { data };
    interface->set_kind_data(KindData::AsItem);
}

bool NumberNode::collect_data()
{
    return true;
}

bool NumberNode::calculate()
{
    return true;
}


