#pragma once
#include "node.hpp"

class ConstructPoint : public Node
{
public:
    ConstructPoint(const std::string &id);
    // {
    //     m_inputs.resize(3);
    //     m_outputs.resize(1);
    //     for (int index = 0; index < 3; ++index)
    //     {
    //         m_inputs[index]->set_kind_data(KindData::AsItem);
    //     }
    //     m_outputs[0]->set_kind_data(KindData::AsItem);
    // }
    virtual ~ConstructPoint() { }

    virtual bool calculate() override;
    // {
    //     std::vector<DataRelation> rels;
    //     calculate_path(rels);
    //     size_t count = rels.size();
    //     for (size_t index = 0; index < count; ++index)
    //     {
    //         DataRelation &rel = rels[index];
    //         size_t item_count = rels[index].m_list_indexs.size();
    //         std::vector<AbstractNodeData*> points;
    //         points.resize(item_count);
    //         for (size_t item_index = 0; item_index < item_count; ++item_index)
    //         {
    //             AbstractNodeData *x_data = m_inputs[0]->m_data[rel.m_paths[0]][rel.m_indexs[item_index]];
    //             AbstractNodeData *y_data = m_inputs[1]->m_data[rel.m_paths[1]][rel.m_indexs[item_index]];
    //             AbstractNodeData *z_data = m_inputs[2]->m_data[rel.m_paths[2]][rel.m_indexs[item_index]];
    //             if (x_data->m_type == "NumberData" && y_data->m_type == "NumberData" && z_data->m_type == "NumberData")
    //             {
    //                 double x = static_cast<NumberData*>(x_data)->m_data;
    //                 double y = static_cast<NumberData*>(y_data)->m_data;
    //                 double z = static_cast<NumberData*>(z_data)->m_data;
    //                 PointData *point = new PointData(x, y, z);
    //                 points[item_index] = point;
    //             }
    //             else
    //             {
    //                 std::cout << "ConstructPoint : error " << item_count << std::endl;
    //             }
    //             m_outputs[0]->set_data(rel, points);
    //         }
    //     }
    //     return true;
    // }
    
};

