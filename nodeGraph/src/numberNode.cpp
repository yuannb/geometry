#include "numberNode.hpp"
#include "NodeGraph.hpp"
#include "parameter.hpp"
#include "contruct_primitive_nurbs_surface.h"
NumberNode::NumberNode(const std::string &id)
{
    m_inputs.clear();
    m_outputs.resize(1);

    m_id = id;
    NodeInterface& interface = m_outputs[0];
    // NumberData *data = new NumberData(num);
    // std::vector<size_t> path_data{0};
    // Path path(path_data);
    // interface.m_data[path_data] = std::vector<AbstractNodeData*> { data };
    interface.set_kind_data(KindData::AsItem);
}

bool NumberNode::collect_data()
{
    return true;
}

bool NumberNode::calculate()
{
    if (true == is_calculate)
    {
        return true;
    }
    //m_outputs.clear();
    for (auto elem : m_inputs[0].m_data)
    {
        m_outputs[0].m_data.insert(std::make_pair(elem.first, elem.second));
    }
    // m_outputs.insert(m_outputs.end(), m_inputs.begin(), m_inputs.end());
    return true;
}

Extrude::Extrude(const std::string& id)
{
    m_inputs.resize(3);
    m_outputs.resize(1);

    m_id = id;
    for (int index = 0; index < 3; ++index)
    {
		auto& node_interface = m_outputs[index];
		node_interface.set_kind_data(KindData::AsItem);

    }
    NodeInterface& interface = m_outputs[0];
	interface.set_kind_data(KindData::AsItem);
}


std::vector<std::vector<AbstractNodeData*>> Extrude::calculate_inner(const NodeIterator& node_iter)
{
    auto curve_it = node_iter.m_node_iter[0].m_iter;
    auto curve_data = curve_it.m_path_iter->second[curve_it.m_index];
    auto it = node_iter.m_node_iter[1].m_iter;
    auto dir_data = it.m_path_iter->second[it.m_index];

    auto dist_it = node_iter.m_node_iter[2].m_iter;
    auto dist_data = dist_it.m_path_iter->second[dist_it.m_index];
    
    std::string curve_id = static_cast<GeometryData*>(curve_data)->m_data;
    tnurbs::curve3d* curve = m_parent->m_curves[curve_id];

    Eigen::Vector3d dir = static_cast<PointData*>(dir_data)->m_data;
    
    double dist = static_cast<NumberData*>(dist_data)->m_data;
    tnurbs::surface<double, 3>* result_surf = nullptr;
	tnurbs::mesh<2>* sruf_mesh = new tnurbs::mesh<2>();
    if (curve->get_type() == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
    {
        tnurbs::nurbs3d* bcurve = static_cast<tnurbs::nurbs3d*>(curve);
        tnurbs::nurbs_surf3d* extrude_surf = new tnurbs::nurbs_surf3d();
        tnurbs::create_cylinder_surface<double, 3, true>(*bcurve, dist, dir, *extrude_surf);
        result_surf = extrude_surf;
		tnurbs::surface_mesh_helper<tnurbs::nurbs_surf3d> mh;
        tnurbs::disc_surface(extrude_surf, mh);
		tnurbs::mesh_help_to_mesh(mh, *sruf_mesh);

    }
    else if (curve->get_type() == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
    {
        tnurbs::nunbs3d* bcurve = static_cast<tnurbs::nunbs3d*>(curve);
        tnurbs::nunbs_surf3d* extrude_surf = new tnurbs::nunbs_surf3d();
        tnurbs::create_cylinder_surface<double, 3, false>(*bcurve, dist, dir, *extrude_surf);
        result_surf = extrude_surf;
		
        tnurbs::surface_mesh_helper<tnurbs::nunbs_surf3d> mh;
        tnurbs::disc_surface(extrude_surf, mh, tnurbs::TDEFAULT_ERROR<double>::value, 0.5, 0.1, 1.0);
		tnurbs::mesh_help_to_mesh(mh, *sruf_mesh);
    }

	std::string gid = create_uuid();
	m_parent->m_surfaces[gid] = result_surf;
	m_parent->m_surface_disc[gid] = sruf_mesh;
	m_parent->m_nid_to_gid[m_id].push_back(gid);
	GeometryData* data = new GeometryData(gid);
	// std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid);
    std::vector<std::vector<AbstractNodeData*>> result{ {nullptr} };
    result[0][0] = data;
    m_parent->m_node_data.push_back(data);
    return result;
}

Ranges::Ranges(const std::string& id)
{
    m_inputs.resize(3);
    m_outputs.resize(1);

    m_id = id;
    for (int index = 0; index < 3; ++index)
    {
		auto& node_interface = m_outputs[index];
		node_interface.set_kind_data(KindData::AsItem);

    }
    NodeInterface& interface = m_outputs[0];
	interface.set_kind_data(KindData::AsList);
}

std::vector<std::vector<AbstractNodeData*>> Ranges::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result;
    std::array<double, 3> params;
    for (int index = 0; index < 3; ++index)
    {
		auto it = node_iter.m_node_iter[index].m_iter;
		auto abstract_data = it.m_path_iter->second[it.m_index];
        if (abstract_data->m_type != NodeDataType::numberData)
        {
            assert(false);
        }
        auto number_data = static_cast<NumberData*>(abstract_data);
        params[index] = number_data->m_data;
    }
    int count = static_cast<int>(params[2]);
    result.resize(1);
    result[0].resize(count + 1);
    double step = (params[1] - params[0]) / (double)count;
    double current_param = params[0];
    for (int index = 0; index < count; ++index)
    {
        NumberData* data = new NumberData(current_param);
		// std::unique_ptr<NumberData> data = std::make_unique<NumberData>(current_param);
        current_param += step;
        result[0][index] = data;
        m_parent->m_node_data.push_back(data);
    }
	NumberData* data = new NumberData(params[1]);
	// std::unique_ptr<NumberData> data = std::make_unique<NumberData>(params[1]);
	result[0][count] = data;
    m_parent->m_node_data.push_back(data);
    return result;
}
