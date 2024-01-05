#include "InterpolateCurve.hpp"
#include "NodeGraph.hpp"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <functional>
#include <random>
std::string create_uuid()
{
    std::stringstream stream;
    auto random_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 seed_engine(random_seed);
    std::uniform_int_distribution<std::size_t> random_gen ;
    std::size_t value = random_gen(seed_engine);
    stream << std::hex << value;
    
    return stream.str();
}


InterpolateCurve::InterpolateCurve(const std::string &id): Node(id)
{
    m_inputs.resize(1);
    NodeInterface *node_interface = new NodeInterface();
    node_interface->set_kind_data(KindData::AsList);
    m_inputs[0] = node_interface;

    m_outputs.resize(1);
    NodeInterface *out_interface = new NodeInterface();
    out_interface->set_kind_data(KindData::AsList);
    m_outputs[0] = out_interface;
}

bool InterpolateCurve::calculate()
{
    std::vector<DataRelation> rels;
    calculate_path(rels);
    size_t count = rels.size();
    for (size_t index = 0; index < count; ++index)
    {
        DataRelation &rel = rels[index];
        // size_t item_count = rel.m_list_indexs.size();
        std::vector<tnurbs::curve3d*> curve;
        size_t points_count = m_inputs[0]->m_data[rel.m_paths[0]].size();
        int degree = 1;
        if (points_count < 2)
        {
            std::cout << "too less points" << std::endl;
            return false;
        }
        if (points_count == 3)
            degree = 2;
        if (points_count > 3)
            degree = 3;
        Eigen::Matrix3Xd points(3, points_count);
        for (size_t i = 0; i < points_count; ++i)
        {
            AbstractNodeData *point = m_inputs[0]->m_data[rel.m_paths[0]][i];
            if (point->m_type == "PointData")
            {
                PointData *pt = static_cast<PointData*>(point);
                points.col(i) = pt->m_data;
            }
            else
            {
                std::cout << "InterpolateCurve : error " << i << std::endl;
            }
        }
        tnurbs::nunbs3d *interpolate_curve = new tnurbs::nunbs3d();
        tnurbs::global_curve_interpolate<double, 3, tnurbs::ENPARAMETERIEDTYPE::CHORD>(points, degree, *interpolate_curve);
        std::string gid = create_uuid();
        m_parent->m_curves[gid] = interpolate_curve;
        tnurbs::curve_mesh_helper<tnurbs::nunbs3d> mh;
        tnurbs::disc_curve(interpolate_curve, mh, 10.0, 10.0);
        tnurbs::mesh<1> *curve_mesh = new tnurbs::mesh<1>();
        tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
        m_parent->m_curves_disc[gid] = curve_mesh;
        GeometryData *data = new GeometryData(gid);
        std::vector<AbstractNodeData*> node_vec { data };

        m_outputs[0]->set_data(rel, node_vec);
    }
    return true;
}

