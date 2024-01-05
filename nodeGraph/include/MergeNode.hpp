#pragma once
#include "node.hpp"


class MergeNode : public Node
{
public:
    MergeNode(const std::string &id, unsigned count);
    // {
    //     m_inputs.resize(4);
    //     m_inputs[0]->set_kind_data(KindData::AsList);
    //     m_inputs[1]->set_kind_data(KindData::AsList);
    //     m_inputs[2]->set_kind_data(KindData::AsList);
    //     m_inputs[3]->set_kind_data(KindData::AsList);
    //     m_outputs.resize(1);
    //     m_inputs[0]->set_kind_data(KindData::AsList);
    // }
    virtual ~MergeNode() { }

    

    virtual bool calculate() override;
    // {

    //     std::vector<DataRelation> rels;
    //     calculate_path(rels);
    //     size_t count = rels.size();
    //     size_t inputs_count = m_inputs.size();
    //     for (size_t index = 0; index < count; ++index)
    //     {
    //         DataRelation &rel = rels[index];
    //         std::vector<AbstractNodeData*> datas;
    //         for (size_t input_index = 0; input_index < inputs_count; ++index)
    //         {
    //             std::vector<AbstractNodeData*> inputs = m_inputs[input_index]->m_data[rel.m_paths[input_index]];
    //             datas.insert(datas.end(), inputs.begin(), inputs.end());
    //         }

    //         m_outputs[0]->set_data(rel, datas);
    //         return true;
    //     }
    // }
    
};

