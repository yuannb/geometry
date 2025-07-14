#include "InterpolateCurve.hpp"
#include "NodeGraph.hpp"
#include "bezier_curve_int.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <functional>
#include <random>


InterpolateCurve::InterpolateCurve(const std::string &id): Node(id)
{
    m_id = id;
    m_inputs.resize(1);
    m_outputs.resize(1);
    auto& node_interface = m_inputs[0];
    node_interface.set_kind_data(KindData::AsList);

    auto& out_node_interface = m_outputs[0];
    out_node_interface.set_kind_data(KindData::AsItem);
}

std::vector<std::vector<AbstractNodeData*>> InterpolateCurve::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result{ {nullptr} };
    
    auto points = node_iter.m_node_iter[0].m_iter.m_path_iter->second;
    
    size_t points_count = points.size();

    Eigen::Matrix3Xd interpolate_points(3, points_count);
    for (size_t index = 0; index != points_count; ++index)
    {
        AbstractNodeData* data = points[index];
        if (data->m_type != "PointData")
        {
            assert(false);
        }
        PointData* point_data = static_cast<PointData*>(data);
        interpolate_points.col(index) = point_data->m_data;
    }
	int degree = 1;
	if (points_count < 2)
	{
		std::cout << "too less points" << std::endl;
		return result;
	}
	if (points_count == 3)
		degree = 2;
	else if (points_count > 3)
		degree = 3;
	tnurbs::nunbs3d* interpolate_curve = new tnurbs::nunbs3d();
	tnurbs::global_curve_interpolate<double, 3, tnurbs::ENPARAMETERIEDTYPE::CHORD>(interpolate_points, degree, *interpolate_curve);
	std::string gid = create_uuid();
	m_parent->m_curves[gid] = interpolate_curve;
	tnurbs::curve_mesh_helper<tnurbs::nunbs3d> mh;
	tnurbs::disc_curve(interpolate_curve, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
	tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
	tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
	m_parent->m_curves_disc[gid] = curve_mesh;
    m_parent->m_nid_to_gid[m_id].push_back(gid);

	GeometryData* data = new GeometryData(gid);
    data->m_geoemtryType = 1;
    result[0][0] = data;
    m_parent->m_node_data.push_back(data);
    return result;
}

CurveIntCurve::CurveIntCurve(const std::string& id)
{
    m_id = id;
    m_inputs.resize(2);
    m_outputs.resize(1);
    auto& node_interface = m_inputs[0];
    node_interface.set_kind_data(KindData::AsItem);

    auto& node_interface2 = m_inputs[1];
    node_interface2.set_kind_data(KindData::AsItem);
    
    auto& out_node_interface = m_outputs[0];
    out_node_interface.set_kind_data(KindData::AsItem);
}

static std::vector<tnurbs::curve_curve_int_point<double, 3>> curve_int_curve(tnurbs::curve3d* curve1, tnurbs::curve3d* curve2)
{
    std::vector<tnurbs::curve_curve_int_point<double, 3>>  int_result;
    if (curve1->get_type() == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
    {
        if (curve2->get_type() == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
        {
            tnurbs::curve_curve_int<tnurbs::nunbs3d, tnurbs::nunbs3d> curves_int;
            curves_int.init(static_cast<tnurbs::nunbs3d*>(curve1), static_cast<tnurbs::nunbs3d*>(curve2));
            curves_int.run_intersect();
            int_result = curves_int.m_int_points;
        }
        else if (curve2->get_type() == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
        {
            tnurbs::curve_curve_int<tnurbs::nunbs3d, tnurbs::nurbs3d> curves_int;
            curves_int.init(static_cast<tnurbs::nunbs3d*>(curve1), static_cast<tnurbs::nurbs3d*>(curve2));
            curves_int.run_intersect();
            int_result = curves_int.m_int_points;
        }
        else
        {
            assert(false);
        }
    }
    else if (curve1->get_type() == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
    {
        if (curve2->get_type() == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
        {
            tnurbs::curve_curve_int<tnurbs::nurbs3d, tnurbs::nunbs3d> curves_int;
            curves_int.init(static_cast<tnurbs::nurbs3d*>(curve1), static_cast<tnurbs::nunbs3d*>(curve2));
            curves_int.run_intersect();
            int_result = curves_int.m_int_points;
        }
        else if (curve2->get_type() == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
        {
            tnurbs::curve_curve_int<tnurbs::nurbs3d, tnurbs::nurbs3d> curves_int;
            curves_int.init(static_cast<tnurbs::nurbs3d*>(curve1), static_cast<tnurbs::nurbs3d*>(curve2));
            curves_int.run_intersect();
            int_result = curves_int.m_int_points;
        }
        else
        {
            assert(false);
        }
    }
    else
    {
        assert(false);
    }
    return int_result;
}

std::vector<std::vector<AbstractNodeData*>> CurveIntCurve::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result;
    std::array<tnurbs::curve3d*, 2> curves;
    for (size_t index = 0; index < 2; ++index)
    {
        AbstractNodeData* data = get_data_item(index, node_iter);
        if (data->m_type != NodeDataType::geometryData)
        {
            assert(false);
        }
        GeometryData* geoData = static_cast<GeometryData*>(data);
        if (geoData->m_geoemtryType != 1)
        {
            assert(false);
        }
        curves[index] = m_parent->m_curves[geoData->m_data];
    }

    std::vector<tnurbs::curve_curve_int_point<double, 3>> int_reuslt = curve_int_curve(curves[0], curves[1]);
    result.resize(1);
    for (auto& intpoint : int_reuslt)
    {

		// std::unique_ptr<PointData> data = std::make_unique<PointData>(intpoint.m_int_point);
        PointData* data = new PointData(intpoint.m_int_point);
		std::string gid = create_uuid();
		m_parent->m_points_disc[gid] = data->m_data;
		m_parent->m_nid_to_gid[m_id].push_back(gid);
        result[0].push_back(data);
		m_parent->m_node_data.push_back(data);
    }
    return result;
}
