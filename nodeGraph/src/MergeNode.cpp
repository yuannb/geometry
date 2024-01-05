#include "MergeNode.hpp"
#include "NodeGraph.hpp"


MergeNode::MergeNode(const std::string &id, unsigned count): Node(id)
{
    m_inputs.resize(count);
    m_outputs.resize(1);
    for (unsigned index = 0; index < count; ++index)
    {
        NodeInterface *node_interface = new NodeInterface();
        node_interface->set_kind_data(KindData::AsList);
        m_inputs[index] = node_interface;        
    }

    m_outputs.resize(1);
    NodeInterface *out_interface = new NodeInterface();
    out_interface->set_kind_data(KindData::AsList);
    m_outputs[0] = out_interface;
}

bool MergeNode::calculate()
{
    std::vector<DataRelation> rels;
    calculate_path(rels);
    size_t count = rels.size();
    size_t inputs_count = m_inputs.size();
    for (size_t index = 0; index < count; ++index)
    {
        DataRelation &rel = rels[index];
        std::vector<AbstractNodeData*> datas;
        for (size_t input_index = 0; input_index < inputs_count; ++input_index)
        {
            std::vector<AbstractNodeData*> inputs = m_inputs[input_index]->m_data[rel.m_paths[input_index]];
            datas.insert(datas.end(), inputs.begin(), inputs.end());
        }

        m_outputs[0]->set_data(rel, datas);
    }
    return true;
}

