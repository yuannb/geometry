#pragma once
#include "node.hpp"
#include "fit_nurbs.h"

class InterpolateCurve :  public Node
{
public:
    InterpolateCurve(const std::string &id);
    // {
    //     m_inputs.resize(1);
    //     m_inputs[0]->set_kind_data(KindData::AsList);
    //     m_outputs.resize(1);
    //     m_inputs[0]->set_kind_data(KindData::AsList);
    // }
    virtual ~InterpolateCurve() { }

    virtual bool calculate() override;
    // {
    //     std::vector<DataRelation> rels;
    //     calculate_path(rels);
    //     size_t count = rels.size();
    //     for (size_t index = 0; index < count; ++index)
    //     {
    //         DataRelation &rel = rels[index];
    //         // size_t item_count = rel.m_list_indexs.size();
    //         std::vector<tnurbs::curve3d*> curve;
    //         size_t points_count = m_inputs[0]->m_data[rel.m_paths[0]].size();
    //         int degree = 1;
    //         if (points_count < 2)
    //         {
    //             std::cout << "too less points" << std::endl;
    //             return false;
    //         }
    //         if (points_count == 3)
    //             degree = 2;
    //         if (points_count > 3)
    //             degree = 3;
    //         Eigen::Matrix3Xd points(3, points_count);
    //         for (size_t i = 0; i < points_count; ++i)
    //         {
    //             AbstractNodeData *point = m_inputs[0]->m_data[rel.m_paths[0]][i];
    //             if (point->m_type == "PointData")
    //             {
    //                 PointData *pt = static_cast<PointData*>(point);
    //                 points.col(i) = pt->m_data;
    //             }
    //             else
    //             {
    //                 std::cout << "InterpolateCurve : error " << i << std::endl;
    //             }
    //         }
    //         tnurbs::nunbs3d *interpolate_curve = new tnurbs::nunbs3d();
    //         tnurbs::global_curve_interpolate<double, 3, tnurbs::ENPARAMETERIEDTYPE::CHORD>(points, degree, *interpolate_curve);

    //         m_parent->m_curves.push_back(std::unique_ptr<tnurbs::nunbs3d>(interpolate_curve));
    //     }
    //     return true;
    // }

};
