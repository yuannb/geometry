#pragma once
#include "node.hpp"
#include "parameter.hpp"
#include <map>
#include <string>
#include <memory>
#include <vector>
#include <set>
#include "fit_nurbs.h"
// #include "logger_manager.hpp"

class NodeGraph
{
public:
    std::string m_id;
    std::map<std::string, Node*>  m_nodes;
    
    std::map<std::string, Eigen::Vector3d> m_points_disc;
    
    std::map<std::string, tnurbs::curve3d*> m_curves;
    std::map<std::string, tnurbs::mesh<1>*> m_curves_disc;

    std::map<std::string, tnurbs::surface<double, 3>*> m_surfaces;
    std::map<std::string, tnurbs::mesh<2>*> m_surface_disc;

    std::map<std::string, std::vector<std::string>> m_nid_to_gid;

    std::vector<std::string> m_nodes_sort;

    std::vector<AbstractNodeData*> m_node_data;


    bool from_database = false;
    std::set<std::string> changed_nodes;
    std::set<std::string> deleted_nodes;
    std::set<std::string> depended_nodes;
    NodeGraph(const std::string &id): m_id(id) { }
    ~NodeGraph(); 

    bool calculate();
    nlohmann::ordered_json to_json() const
    {
        nlohmann::ordered_json node_json;
        node_json["id"] = m_id;
        node_json["node_data"] = {};
        for (auto it = m_nodes.begin(); it != m_nodes.end(); ++it)
        {
            node_json["node_data"][it->first] = it->second->to_json();
        }
        
        node_json["curve_disc"] = {};
        // auto& disc = node_json["curve_disc"];
        for (auto it = m_curves_disc.begin(); it != m_curves_disc.end(); ++it)
        {
            nlohmann::ordered_json disc = {};
            disc["type"] = 1;
            disc["points"] = nlohmann::ordered_json::array();
            std::vector<double> points;
            auto& ps = it->second->m_ders;
            for (auto p : ps)
            {
                Eigen::Vector3d& tp = p(0, 0);
                disc["points"].push_back(tp[0]);
                disc["points"].push_back(tp[1]);
                disc["points"].push_back(tp[2]);
            }
            disc["indexs"] = {};
            for (auto index : it->second->m_indexs)
            {
                std::copy(index.begin(), index.end(), std::back_inserter(disc["indexs"]));
            }
            node_json["curve_disc"][it->first] = disc;
        }
        
        node_json["surf_disc"] = {};
        // auto& surf_disc = node_json["surf_disc"];
        for (auto it = m_surface_disc.begin(); it != m_surface_disc.end(); ++it)
        {
            nlohmann::ordered_json surf_disc = {};
            surf_disc["type"] = 2;
            surf_disc["points"] = nlohmann::ordered_json::array();
            std::vector<double> points;
            auto& ps = it->second->m_ders;
            for (auto p : ps)
            {
                Eigen::Vector3d& tp = p(0, 0);
                surf_disc["points"].push_back(tp[0]);
                surf_disc["points"].push_back(tp[1]);
                surf_disc["points"].push_back(tp[2]);
            }
            surf_disc["indexs"] = {};
            for (auto index : it->second->m_indexs)
            {
                surf_disc["indexs"].emplace_back(index);
            }
            node_json["surf_disc"][it->first] = surf_disc;
        }
        
        node_json["points"] = {};
        // auto& points_array = node_json["points"];
        for (auto it = m_points_disc.begin(); it != m_points_disc.end(); ++it)
        {
            nlohmann::ordered_json points_array = {};
            points_array["type"] = 0;
            auto& point = it->second;
            points_array["points"] = {point[0], point[1], point[2]};
            node_json["points"][it->first] = points_array;
        }
        return node_json;
    }

};


template<bool is_rational>
inline bool bcurve_to_json(const tnurbs::nurbs_curve<double, 3, is_rational, -1, -1>& nurbs_curve, nlohmann::ordered_json& curve_json)
{
    curve_json[ParamterType::curveType] = magic_enum::enum_name(nurbs_curve.get_type());
    const Eigen::VectorX<double>& knots_vec = nurbs_curve.get_knots_ref();
    std::vector<double> knots(knots_vec.data(), knots_vec.data() + knots_vec.size());
    curve_json[ParamterType::knots] = nlohmann::ordered_json(knots);
    curve_json[ParamterType::degree] = nurbs_curve.get_degree();

    auto control_points = nurbs_curve.get_control_points_ref();

	std::vector<double> elements(control_points.data(), control_points.data() + control_points.size());
	curve_json[ParamterType::controlPoints] = nlohmann::ordered_json(elements);
    return true;
}

nlohmann::ordered_json inline curve_to_json(const tnurbs::curve3d* curve)
{
    nlohmann::ordered_json curve_json;
    curve_json[ParamterType::geoemtryType] = 1;
    if (curve->get_type() == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
    {
        const tnurbs::nunbs3d* bcurve = static_cast<const tnurbs::nunbs3d*>(curve);
        bcurve_to_json<false>(*bcurve, curve_json);
    }
    else if (curve->get_type() == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
    {
        const tnurbs::nurbs3d* bcurve = static_cast<const tnurbs::nurbs3d*>(curve);
        bcurve_to_json<true>(*bcurve, curve_json);
    }
    else
    {
        assert(false);
    }
    return curve_json;
}

template<bool is_rational>
inline bool bsurf_to_json(const tnurbs::nurbs_surface<double, 3, -1, -1, -1, -1, is_rational>& nurbs_surf, nlohmann::ordered_json& surf_json)
{
    surf_json[ParamterType::surfType] = magic_enum::enum_name(nurbs_surf.get_type());
   
    auto u_knots_vec = nurbs_surf.get_u_knots_ref();
    auto v_knots_vec = nurbs_surf.get_v_knots_ref();
    std::vector<double> u_knots(u_knots_vec.data(), u_knots_vec.data() + u_knots_vec.size());
    surf_json[ParamterType::uKnots] = nlohmann::ordered_json(u_knots);

    std::vector<double> v_knots(v_knots_vec.data(), v_knots_vec.data() + v_knots_vec.size());
    surf_json[ParamterType::vKnots] = nlohmann::ordered_json(v_knots);
    
    surf_json[ParamterType::uDegree] = nurbs_surf.get_u_degree();
    surf_json[ParamterType::vDegree] = nurbs_surf.get_v_degree();
    
    auto control_points = nurbs_surf.get_control_points_ref();

    surf_json[ParamterType::controlPoints] = nlohmann::ordered_json::array();
    Eigen::Index cols = control_points[0].cols();
    size_t rows = control_points.size();

    for (size_t row_index = 0; row_index < rows; ++row_index)
    {
        std::vector<double> elements(control_points[row_index].data(), control_points[row_index].data() + control_points[row_index].size());
        surf_json[ParamterType::controlPoints].push_back(nlohmann::ordered_json(elements));
    }

    return true;
}

nlohmann::ordered_json inline surf_to_json(const tnurbs::surf3d* surf)
{
    nlohmann::ordered_json surf_json;
    surf_json[ParamterType::geoemtryType] = 2;
    if (surf->get_type() == tnurbs::ENGEOMETRYTYPE::NUNBS_SURF)
    {
        const tnurbs::nunbs_surf3d* bsurf = static_cast<const tnurbs::nunbs_surf3d*>(surf);
        bsurf_to_json<false>(*bsurf, surf_json);
    }
    else if (surf->get_type() == tnurbs::ENGEOMETRYTYPE::NURBS_SURF)
    {
        const tnurbs::nurbs_surf3d* bsurf = static_cast<const tnurbs::nurbs_surf3d*>(surf);
        bsurf_to_json<true>(*bsurf, surf_json);
    }
    else
    {
        assert(false);
    }
    return surf_json;
}

inline bool json_to_knots(const nlohmann::ordered_json& js, Eigen::VectorXd& knots)
{
    size_t knots_count = js.size();
    knots.resize(knots_count);
    for (size_t i = 0; i < knots_count; ++i)
    {
        knots[i] = js[i].get<double>();
    }
    return true;
}

template<bool is_rational>
inline tnurbs::nurbs_curve<double, 3, is_rational, -1, -1>* json_to_bcurve(const nlohmann::ordered_json& curve_json)
{
    constexpr int point_size = is_rational ? 4 : 3;
    using curve_type = tnurbs::nurbs_curve<double, 3, is_rational, -1, -1>;
    using CPVType = Eigen::Matrix<double, point_size, Eigen::Dynamic>;
    Eigen::VectorXd knots;
    json_to_knots(curve_json[ParamterType::knots], knots);
    nlohmann::ordered_json cp_js = curve_json[ParamterType::controlPoints];
	std::vector<double> points_vec = cp_js.get<std::vector<double>>();
	size_t cols = points_vec.size() / point_size;
	if (cols * point_size != points_vec.size())
	{
		// LoggerManager::getLogger()->error("matrix data is invalid: json_to_bsurf()");
		return nullptr;
	}
	CPVType points = Eigen::Map<CPVType>(points_vec.data(), point_size, cols);

    curve_type* curve = new curve_type(knots, points);
    return curve;
}

inline tnurbs::curve3d* json_to_curve(const nlohmann::ordered_json& curve_js)
{
    auto curve_type = magic_enum::enum_cast<tnurbs::ENGEOMETRYTYPE> (curve_js[ParamterType::curveType].get<std::string>());
    if (curve_type == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
    {
        return json_to_bcurve<false>(curve_js);
    }
    else if (curve_type == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
    {
        return json_to_bcurve<true>(curve_js);
    }
    else
    {
        assert(false);
    }
    return nullptr;
}


template<bool is_rational>
inline tnurbs::nurbs_surface<double, 3, -1, -1, -1, -1, is_rational>* json_to_bsurf(const nlohmann::ordered_json& surf_json)
{
    // surf_json[ParamterType::surfType] = magic_enum::enum_name(nurbs_surf.get_type());
    constexpr int point_size = is_rational ? 4 : 3;
    using surf_type = tnurbs::nurbs_surface<double, 3, -1, -1, -1, -1, is_rational>;
    using CPType = surf_type::CPType;
    using CPVType = Eigen::Matrix<double, point_size, Eigen::Dynamic>;
    Eigen::VectorXd u_knots, v_knots;
    json_to_knots(surf_json[ParamterType::uKnots], u_knots);
    json_to_knots(surf_json[ParamterType::vKnots], v_knots);
    CPType control_points;
    nlohmann::ordered_json cp_js = surf_json[ParamterType::controlPoints];
    size_t rows = cp_js.size();
    for (size_t i = 0; i < rows; ++i)
    {
        std::vector<double> points_vec = cp_js[i].get<std::vector<double>>();
        size_t cols = points_vec.size() / point_size;
        if (cols * point_size != points_vec.size())
        {
            // LoggerManager::getLogger()->error("matrix data is invalid: json_to_bsurf()");
            return nullptr;
        }
        CPVType points = Eigen::Map<CPVType>(points_vec.data(), point_size, cols);
        control_points.push_back(std::move(points));
    }
    
    surf_type* surf = new surf_type(u_knots, v_knots, control_points);
    return surf;
}


inline tnurbs::surf3d* json_to_surf(const nlohmann::ordered_json& surf_js)
{
    std::string surf_type_str = surf_js[ParamterType::surfType].get<std::string>();
    auto surf_type = magic_enum::enum_cast<tnurbs::ENGEOMETRYTYPE> (surf_type_str);
    if (surf_type== tnurbs::ENGEOMETRYTYPE::NUNBS_SURF)
    {
        return json_to_bsurf<false>(surf_js);
    }
    else if (surf_type == tnurbs::ENGEOMETRYTYPE::NURBS_SURF)
    {
        return json_to_bsurf<true>(surf_js);
    }
    else
    {
        assert(false);
    }
    return nullptr;
}

inline nlohmann::ordered_json point_to_json(const Eigen::Vector3d& point)
{
    nlohmann::ordered_json point_json;
    point_json[ParamterType::geoemtryType] = 0;
    point_json["point"] = { point[0], point[1], point[2] };
    return point_json;
}

inline bool json_to_point(nlohmann::ordered_json&js, Eigen::Vector3d& point)
{
    nlohmann::ordered_json point_js = js[ParamterType::point];
    for (int index = 0; index < 3; ++index)
    {
        point[index] = point_js[index];
    }
    return true;
}
// tnurbs::mesh<1>*

template<bool is_rational>
inline tnurbs::mesh<1>* bcurve_disc(const tnurbs::nurbs_curve<double, 3, is_rational, -1, -1>& bcurve)
{
    using curve_type = tnurbs::nurbs_curve<double, 3, is_rational, -1, -1>;
	tnurbs::curve_mesh_helper<curve_type> mh;
	tnurbs::disc_curve(&bcurve, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
	tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
	tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
    return curve_mesh;
}

inline tnurbs::mesh<1>* curve_disc(const tnurbs::curve3d* curve)
{
    auto curve_type = curve->get_type();
    if (curve_type == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
    {
        const tnurbs::nunbs3d* bcurve = static_cast<const tnurbs::nunbs3d*>(curve);
		return bcurve_disc<false>(*bcurve);
    }
    else if (curve_type == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
    {
        const tnurbs::nurbs3d* bcurve = static_cast<const tnurbs::nurbs3d*>(curve);
        return bcurve_disc<true>(*bcurve);
    }
    else
    {
        assert(false);
    }
    return nullptr;
}

template<bool is_rational>
inline tnurbs::mesh<2>* bsurf_disc(const tnurbs::nurbs_surface<double, 3, -1, -1, -1, -1, is_rational>& bsurf)
{
    using surf_type = tnurbs::nurbs_surface<double, 3, -1, -1, -1, -1, is_rational>;
	tnurbs::surface_mesh_helper<surf_type> mh;
	tnurbs::disc_surface(&bsurf, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
	tnurbs::mesh<2>* surf_mesh = new tnurbs::mesh<2>();
	tnurbs::mesh_help_to_mesh(mh, *surf_mesh);
    return surf_mesh;
}

inline tnurbs::mesh<2>* surf_disc(const tnurbs::surf3d* surf)
{
    auto surf_type = surf->get_type();
    if (surf_type == tnurbs::ENGEOMETRYTYPE::NUNBS_SURF)
    {
        const tnurbs::nunbs_surf3d* bsurf = static_cast<const tnurbs::nunbs_surf3d*>(surf);
		return bsurf_disc<false>(*bsurf);
    }
    else if (surf_type == tnurbs::ENGEOMETRYTYPE::NURBS_SURF)
    {
        const tnurbs::nurbs_surf3d* bsurf = static_cast<const tnurbs::nurbs_surf3d*>(surf);
		return bsurf_disc<true>(*bsurf);
    }
    else
    {
        assert(false);
    }
    return nullptr;
}
