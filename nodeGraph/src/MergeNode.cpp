#include "MergeNode.hpp"
#include "NodeGraph.hpp"
#include "bezier_curve_int.h"
#include "nurbs_build.h"
#include "debug_used.h"

MergeNode::MergeNode(const std::string &id, size_t count): Node(id)
{
    m_id = id;
    m_inputs.resize(count);
    m_outputs.resize(1);
    for (size_t index = 0; index < count; ++index)
    {
        auto& node_interface = m_inputs[index];
        node_interface.set_kind_data(KindData::AsList);
    }

	auto& out_node_interface = m_outputs[0];
	out_node_interface.set_kind_data(KindData::AsList);
}

bool MergeNode::calculate()
{
    if (is_calculate == true)
    {
        return true;
    }
    size_t count = m_inputs.size();
	auto& out_node_interface = m_outputs[0];
    for (size_t index = 0; index < count; ++index)
    {
        auto& node_interface = m_inputs[index];
        for (auto& element : node_interface.m_data)
        {
            auto& vecs = out_node_interface.m_data[element.first];
            auto& eles = element.second;
            vecs.insert(vecs.end(), eles.begin(), eles.end());
        }
    }
    for (auto& input : m_inputs)
    {
        input.m_data.clear();
    }
    return true;
}

Ellipse::Ellipse(const std::string& id)
{
    m_id = id;
    m_inputs.resize(5);
    m_outputs.resize(1);
    for (size_t index = 0; index < 5; ++index)
    {
        auto& node_interface = m_inputs[index];
        node_interface.set_kind_data(KindData::AsItem);
    }

	auto& out_node_interface = m_outputs[0];
	out_node_interface.set_kind_data(KindData::AsItem);
}

std::vector<std::vector<AbstractNodeData*>> Ellipse::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result{ {nullptr} };
    Eigen::Vector3d normal;

	auto iter = node_iter.m_node_iter[0].m_iter;
	AbstractNodeData* coord = iter.m_path_iter->second[iter.m_index];
	if (coord->m_type != NodeDataType::pointData)
	{
		assert(false);
	}
	auto origin_data = static_cast<PointData*>(coord);
    Eigen::Vector3d origin = origin_data->m_data;

	iter = node_iter.m_node_iter[1].m_iter;
	coord = iter.m_path_iter->second[iter.m_index];
	if (coord->m_type != NodeDataType::pointData)
	{
		assert(false);
	}

    Eigen::AngleAxisd rotateVector(0.9, Eigen::Vector3d(-1, 0, 0));


	auto x_dir_data = static_cast<PointData*>(coord);
    Eigen::Vector3d x_dir = x_dir_data->m_data;
    x_dir.normalize();
    if (x_dir.squaredNorm() < 0.1)
    {
        assert(false);
    }
	
    iter = node_iter.m_node_iter[2].m_iter;
	coord = iter.m_path_iter->second[iter.m_index];
	if (coord->m_type != NodeDataType::pointData)
	{
		assert(false);
	}
	auto y_dir_data = static_cast<PointData*>(coord);
    Eigen::Vector3d y_dir = y_dir_data->m_data;
    y_dir.normalize();
    if (y_dir.squaredNorm() < 0.1)
    {
        assert(false);
    }

    std::array<double, 2> radius{};
    for (size_t index = 3; index < 5; ++index)
    {
        auto iter = node_iter.m_node_iter[index].m_iter;
        AbstractNodeData* coord = iter.m_path_iter->second[iter.m_index];
        if (coord->m_type != NodeDataType::numberData)
        {
            assert(false);
        }
        auto number_data = static_cast<NumberData*>(coord);
        radius[index - 3] = number_data->m_data;
    }
    if (radius[0] < radius[1])
    {
        assert(false);
    }

    tnurbs::nurbs3d* ellipse = new tnurbs::nurbs3d();
    tnurbs::create_nurbs_ellipse(origin, x_dir, y_dir, radius[0], radius[1], 0.0, M_PI * 2, *ellipse);

	std::string gid = create_uuid();
	m_parent->m_curves[gid] = ellipse;
	m_parent->m_nid_to_gid[m_id].push_back(gid);
	tnurbs::curve_mesh_helper<tnurbs::nurbs3d> mh;
	tnurbs::disc_curve(ellipse, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
	tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
	tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
	m_parent->m_curves_disc[gid] = curve_mesh;
	// std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid);
	GeometryData* data = new GeometryData(gid);
    data->m_geoemtryType = 1;
    result[0][0] = data;
    m_parent->m_node_data.push_back(data);
    return result;
}

Rotate::Rotate(const std::string& id)
{
    m_id = id;
    m_inputs.resize(4);
    m_outputs.resize(1);
    for (size_t index = 0; index < 4; ++index)
    {
        auto& node_interface = m_inputs[index];
        node_interface.set_kind_data(KindData::AsItem);
    }

	auto& out_node_interface = m_outputs[0];
	out_node_interface.set_kind_data(KindData::AsItem);
}

// 函数功能：使点绕任意轴旋转
// 参数说明：
// - point：待旋转的点
// - axisPoint：旋转轴上的一个点
// - axisDirection：旋转轴的方向向量（需为单位向量）
// - angle：旋转角度（单位：弧度）
// 返回值：旋转后的点坐标
Eigen::Vector3d rotatePointAroundAxis(const Eigen::Vector3d& point,
    const Eigen::Vector3d& axisPoint,
    const Eigen::Vector3d& axisDirection,
    double angle) {
    // 1. 将点平移至旋转轴经过原点
    Eigen::Vector3d translatedPoint = point - axisPoint;

    // 2. 构建旋转矩阵
    angle *= M_PI;
    Eigen::AngleAxisd rotation(angle, axisDirection);

    // 3. 执行旋转操作
    Eigen::Vector3d rotatedPoint = rotation * translatedPoint;

    // 4. 将点平移回原来的位置
    return rotatedPoint + axisPoint;
}
std::vector<std::vector<AbstractNodeData*>> Rotate::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result{ {nullptr} };

    auto iter = node_iter.m_node_iter[0].m_iter;
    AbstractNodeData* geoData = iter.m_path_iter->second[iter.m_index];
    if (geoData->m_type != NodeDataType::geometryData)
    {
        assert(false);
    }
    auto geo_data = static_cast<GeometryData*>(geoData);



    iter = node_iter.m_node_iter[1].m_iter;
    AbstractNodeData* angle_data = iter.m_path_iter->second[iter.m_index];
    if (angle_data->m_type != NodeDataType::numberData)
    {
        assert(false);
    }

    NumberData* angle = static_cast<NumberData*>(angle_data);

    iter = node_iter.m_node_iter[2].m_iter;
    auto origin_data = iter.m_path_iter->second[iter.m_index];
    if (origin_data->m_type != NodeDataType::pointData)
    {
        assert(false);
    }
    auto origin = static_cast<PointData*>(origin_data);
    
    iter = node_iter.m_node_iter[3].m_iter;
    auto abstract_data = iter.m_path_iter->second[iter.m_index];
    if (abstract_data->m_type != NodeDataType::pointData)
    {
        assert(false);
    }
    auto normal_data = static_cast<PointData*>(abstract_data);
    Eigen::Vector3d normal = normal_data->m_data;
    normal.normalize();
    if (normal.squaredNorm() < 0.1)
    {
        assert(false);
    }

    std::string gid = geo_data->m_data;
    int type = geo_data->m_geoemtryType;

    if (type == 0)
    {
        Eigen::Vector3d point = m_parent->m_points_disc[gid];
        Eigen::Vector3d rotate_point = rotatePointAroundAxis(point, origin->m_data, normal, angle->m_data);
		std::string gid2 = create_uuid();
        // std::unique_ptr<GeometryData> result_data = std::make_unique<GeometryData>(gid);
        GeometryData* result_data = new GeometryData(gid2);
        result_data->m_data = gid2;
		m_parent->m_nid_to_gid[m_id].push_back(gid2);
        result_data->m_geoemtryType = 0;

		m_parent->m_points_disc[gid2] = rotate_point;
		result[0][0] = result_data;
		m_parent->m_node_data.push_back(result_data);
    }
    else if (type == 1)
    {
        tnurbs::curve3d* curve = m_parent->m_curves[gid];
        if (curve->get_type() == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
        {
            tnurbs::nunbs3d* nurbs_curve = static_cast<tnurbs::nunbs3d*>(curve);
            tnurbs::nunbs3d* rotate_curve = new tnurbs::nunbs3d(*nurbs_curve);
            auto& control_points = rotate_curve->get_control_points_ref();
            Eigen::Index points_count = control_points.cols();
            for (Eigen::Index index = 0; index < points_count; ++index)
            {
                Eigen::Vector3d point = control_points.col(index);
				Eigen::Vector3d rotate_point = rotatePointAroundAxis(point, origin->m_data, normal, angle->m_data);
                control_points.col(index) = rotate_point;
            }
			std::string gid2 = create_uuid();
			m_parent->m_curves[gid2] = rotate_curve;
			tnurbs::curve_mesh_helper<tnurbs::nunbs3d> mh;
			tnurbs::disc_curve(rotate_curve, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
			tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
			tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
			m_parent->m_curves_disc[gid2] = curve_mesh;
			// std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid2);
			m_parent->m_nid_to_gid[m_id].push_back(gid2);
			GeometryData* data = new GeometryData(gid2);
			data->m_geoemtryType = 1;
			result[0][0] = data;
            m_parent->m_node_data.push_back(data);
        }
        else if (curve->get_type() == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
        {
            tnurbs::nurbs3d* nurbs_curve = static_cast<tnurbs::nurbs3d*>(curve);
            tnurbs::nurbs3d* rotate_curve = new tnurbs::nurbs3d(*nurbs_curve);
            auto& control_points = rotate_curve->get_control_points_ref();
            Eigen::Index points_count = control_points.cols();
            for (Eigen::Index index = 0; index < points_count; ++index)
            {
                Eigen::Vector4d point4d = control_points.col(index);
                Eigen::Vector3d point3d = point4d.head<3>() / point4d[3];
				Eigen::Vector3d rotate_point3d = rotatePointAroundAxis(point3d, origin->m_data, normal, angle->m_data);
                point4d.head<3>() = (rotate_point3d * point4d[3]);
                control_points.col(index) = point4d;
            }
			std::string gid2 = create_uuid();
			m_parent->m_curves[gid2] = rotate_curve;
			tnurbs::curve_mesh_helper<tnurbs::nurbs3d> mh;
			tnurbs::disc_curve(rotate_curve, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
			tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
			tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
			m_parent->m_curves_disc[gid2] = curve_mesh;
			m_parent->m_nid_to_gid[m_id].push_back(gid2);
			GeometryData* data = new GeometryData(gid2);
			// std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid2);
			data->m_geoemtryType = 1;
			result[0][0] = data;
            m_parent->m_node_data.push_back(data);
        }
        else
        {
            assert(false);
        }
    }
    else if (type == 2)
    {
        assert(false);
    }
    else
    {
        assert(false);
    }
    return result;
}

Scale::Scale(const std::string& id)
{
    m_id = id;
    size_t count = 3;
    m_inputs.resize(count);
    m_outputs.resize(1);
    for (size_t index = 0; index < count; ++index)
    {
        auto& node_interface = m_inputs[index];
        node_interface.set_kind_data(KindData::AsItem);
    }

	auto& out_node_interface = m_outputs[0];
	out_node_interface.set_kind_data(KindData::AsItem);
}

std::vector<std::vector<AbstractNodeData*>> Scale::calculate_inner(const NodeIterator& node_iter)
{
    auto iter = node_iter.m_node_iter[0].m_iter;
    AbstractNodeData* geoData = iter.m_path_iter->second[iter.m_index];
    if (geoData->m_type != NodeDataType::geometryData)
    {
        assert(false);
    }
    auto geo_data = static_cast<GeometryData*>(geoData);

    iter = node_iter.m_node_iter[2].m_iter;
    auto origin_data = iter.m_path_iter->second[iter.m_index];
    if (origin_data->m_type != NodeDataType::pointData)
    {
        assert(false);
    }
    auto origin = static_cast<PointData*>(origin_data);
    
    iter = node_iter.m_node_iter[1].m_iter;
    auto abstract_data = iter.m_path_iter->second[iter.m_index];
    if (abstract_data->m_type != NodeDataType::numberData)
    {
        assert(false);
    }
    auto factor = static_cast<NumberData*>(abstract_data);
    Eigen::Vector3d factor_vec{ factor->m_data, factor->m_data, factor->m_data };

    Eigen::Matrix3d mat = Eigen::Matrix3d::Identity();
    std::vector<std::vector<AbstractNodeData*>> result{ {nullptr} };
    std::string gid = geo_data->m_data;
    int type = geo_data->m_geoemtryType;
    if (type == 0)
    {
        assert(false);
    }
    else if (type == 1)
    {
        tnurbs::curve3d* curve = m_parent->m_curves[gid];
        if (curve->get_type() == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
        {
            tnurbs::nunbs3d* nurbs_curve = static_cast<tnurbs::nunbs3d*>(curve);
            tnurbs::nunbs3d* rotate_curve = new tnurbs::nunbs3d(*nurbs_curve);
            rotate_curve->scale_curve(origin->m_data, mat, factor_vec);
            std::string gid2 = create_uuid();
            m_parent->m_curves[gid2] = rotate_curve;
            tnurbs::curve_mesh_helper<tnurbs::nunbs3d> mh;
            tnurbs::disc_curve(rotate_curve, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
            tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
            tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
            m_parent->m_curves_disc[gid2] = curve_mesh;
			m_parent->m_nid_to_gid[m_id].push_back(gid2);
			// std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid2);
            GeometryData* data = new GeometryData(gid2);
            data->m_geoemtryType = 1;
            result[0][0] = data;
            m_parent->m_node_data.push_back(data);
        }
        else if (curve->get_type() == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
        {
            tnurbs::nurbs3d* nurbs_curve = static_cast<tnurbs::nurbs3d*>(curve);
            tnurbs::nurbs3d* rotate_curve = new tnurbs::nurbs3d(*nurbs_curve);
            rotate_curve->scale_curve(origin->m_data, mat, factor_vec);
            std::string gid2 = create_uuid();
            m_parent->m_curves[gid2] = rotate_curve;
            tnurbs::curve_mesh_helper<tnurbs::nurbs3d> mh;
            tnurbs::disc_curve(rotate_curve, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
            tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
            tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
            m_parent->m_curves_disc[gid2] = curve_mesh;
			m_parent->m_nid_to_gid[m_id].push_back(gid2);
            GeometryData* data = new GeometryData(gid2);
			// std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid2);
            data->m_geoemtryType = 1;
            result[0][0] = data;
            m_parent->m_node_data.push_back(data);
        }
    }
	else
	{
		assert(false);
	}
    return result;
}

CurveSplit::CurveSplit(const std::string& id)
{
    m_id = id;
    size_t count = 4;
    m_inputs.resize(count);
    m_outputs.resize(1);
    for (size_t index = 0; index < count; ++index)
    {
        auto& node_interface = m_inputs[index];
        node_interface.set_kind_data(KindData::AsItem);
    }

	auto& out_node_interface = m_outputs[0];
	out_node_interface.set_kind_data(KindData::AsList);
}

template<bool is_ratioanl>
static std::vector<tnurbs::nurbs_curve<double, 3, is_ratioanl, -1, -1>*> curve_split_by_plane(tnurbs::nurbs_curve<double, 3, is_ratioanl, -1, -1>* curve, tnurbs::nunbs_surf3d* plane, double tolerance)
{
    using curve_type = tnurbs::nurbs_curve<double, 3, is_ratioanl, -1, -1>;
    std::vector<curve_type*> curves_segments;

	tnurbs::curve_surface_int<curve_type, tnurbs::nunbs_surf3d> curve_surf_intsection;
	curve_surf_intsection.init(curve, plane, tolerance);
	curve_surf_intsection.run_intersect();
	size_t count = curve_surf_intsection.m_int_points.size();
	std::vector<double> int_points_params(count);
	if (count < 2)
	{
		assert(false);
	}
	for (size_t index = 0; index < count; ++index)
	{
		int_points_params[index] = curve_surf_intsection.m_int_points[index].m_int_param[0];
	}
	std::sort(int_points_params.begin(), int_points_params.end());
	bool is_closed = curve->is_closed();
	tnurbs::Box<double, 1> subbox;
	for (size_t index = 0; index < count - 1; ++index)
	{
		subbox.Min[0] = int_points_params[index];
		subbox.Max[0] = int_points_params[index + 1];
		double end_param = int_points_params[index + 1];

		curve_type* sub_curve = new curve_type();
		curve->sub_divide(subbox, *sub_curve);
		curves_segments.push_back(sub_curve);
	}
	tnurbs::Box<double, 1> domain = curve->get_domain_box();
	subbox.Min[0] = int_points_params[count - 1];
	subbox.Max[0] = domain.Max[0];
	curve_type* sub_curve_right = nullptr;
	if (subbox.Max[0] - subbox.Min[0] > 1e-8)
	{
		sub_curve_right = new curve_type();
		curve->sub_divide(subbox, *sub_curve_right);
	}
	subbox.Min[0] = domain.Min[0];
	subbox.Max[0] = int_points_params[0];
	curve_type* sub_curve_left = nullptr;
	if (subbox.Max[0] - subbox.Min[0] > 1e-8)
	{
		sub_curve_left = new curve_type();
		curve->sub_divide(subbox, *sub_curve_left);
	}
	if (is_closed == false)
	{
        if (sub_curve_left != nullptr)
        {
            curves_segments.insert(curves_segments.begin(), sub_curve_left);
        }
        if (sub_curve_right != nullptr)
        {
            curves_segments.push_back(sub_curve_right);
        }
        return curves_segments;
	}
    // else
    if (sub_curve_left != nullptr && sub_curve_right != nullptr)
    {
		int degree = curve->get_degree();
        constexpr int point_size = is_ratioanl ? 4 : 3;
		Eigen::Matrix<double, point_size, Eigen::Dynamic> new_control_points;
		Eigen::VectorXd new_knots;
		Eigen::VectorXd right_knots = sub_curve_left->get_knots_vector();
		for (auto& ele : right_knots)
		{
			ele += domain.Max[0];
		}
		tnurbs::merge_two_curve(degree, degree, sub_curve_right->get_knots_vector(), right_knots, sub_curve_left->get_control_points(), sub_curve_right->get_control_points(), new_knots, new_control_points);
		sub_curve_right->set_knots_vector(new_knots);
		sub_curve_right->set_control_points(new_control_points);
		curves_segments.push_back(sub_curve_right);
        delete sub_curve_left;
        sub_curve_left = nullptr;
    }
    else if (sub_curve_left != nullptr)
    {
		curves_segments.push_back(sub_curve_left);
    }
    else if (sub_curve_right != nullptr)
    {
		curves_segments.push_back(sub_curve_right);
    }
	return curves_segments;
}


std::vector<std::vector<AbstractNodeData*>> CurveSplit::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result;
    result.resize(1);

    auto iter = node_iter.m_node_iter[0].m_iter;
    AbstractNodeData* geoData = iter.m_path_iter->second[iter.m_index];
    if (geoData->m_type != NodeDataType::geometryData)
    {
        assert(false);
    }
    auto geo_data = static_cast<GeometryData*>(geoData);
    if (geo_data->m_geoemtryType != 1)
    {
        assert(false);
    }

    tnurbs::curve3d* curve = m_parent->m_curves[geo_data->m_data];

    // origin
    iter = node_iter.m_node_iter[1].m_iter;
    AbstractNodeData* abstractData = iter.m_path_iter->second[iter.m_index];
    if (abstractData->m_type != NodeDataType::pointData)
    {
        assert(false);
    }
    auto pointData = static_cast<PointData*>(abstractData);


    // normal 
    iter = node_iter.m_node_iter[2].m_iter;
    abstractData = iter.m_path_iter->second[iter.m_index];
    if (abstractData->m_type != NodeDataType::pointData)
    {
        assert(false);
    }
    auto normalData = static_cast<PointData*>(abstractData);

    // tolerance
    iter = node_iter.m_node_iter[3].m_iter;
    abstractData = iter.m_path_iter->second[iter.m_index];
    if (abstractData->m_type != NodeDataType::numberData)
    {
        assert(false);
    }
    auto toleranceData = static_cast<NumberData*>(abstractData);
    tnurbs::ENGEOMETRYTYPE curve_type = curve->get_type();

    Eigen::Vector3d normal = normalData->m_data;
    normal.normalize();
    if (normal.squaredNorm() < 0.1)
    {
        assert(false);
    }
    Eigen::Vector3d udir = Eigen::Vector3d::Zero();
    if (std::abs(normal[0]) > 0.5)
    {
        udir[0] = normal[1];
        udir[1] = -normal[0];
    }
    else
    {
        udir[1] = normal[2];
        udir[2] = -normal[1];
    }

    udir.normalize();

    Eigen::Vector3d vdir = normal.cross(udir);
    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(2);
    control_points[0].resize(3, 2);
    control_points[1].resize(3, 2);
    double length = 500.0;
    control_points[0].col(0) = pointData->m_data - length * udir - length * vdir;
    control_points[0].col(1) = pointData->m_data + length * udir - length * vdir;

    control_points[1].col(0) = pointData->m_data - length * udir + length * vdir;
    control_points[1].col(1) = pointData->m_data + length * udir + length * vdir;
    Eigen::VectorXd knots_vector(4);
    knots_vector[0] = knots_vector[1] = 0.0;
    knots_vector[2] = knots_vector[3] = 1.0;

    tnurbs::nunbs_surf3d* plane = new  tnurbs::nunbs_surf3d(knots_vector, knots_vector, control_points);
    std::unique_ptr<tnurbs::nunbs_surf3d> plane_ptr(plane);

	std::vector<double> int_points_params;

    if (curve_type == tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
    {
        tnurbs::nurbs3d* nurbs = static_cast<tnurbs::nurbs3d*>(curve);
        auto curves_segments = curve_split_by_plane<true>(nurbs, plane, toleranceData->m_data);
        for (const auto curve_seg : curves_segments)
        {
			std::string gid = create_uuid();
			m_parent->m_curves[gid] = curve_seg;
			tnurbs::curve_mesh_helper<tnurbs::nurbs3d> mh;
			tnurbs::disc_curve(curve_seg, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
			tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
			tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
			m_parent->m_curves_disc[gid] = curve_mesh;
			m_parent->m_nid_to_gid[m_id].push_back(gid);
			GeometryData* data = new GeometryData(gid);
			data->m_geoemtryType = 1;
			result[0].push_back(data);
			m_parent->m_node_data.push_back(data);
        }
    }
    else if (curve_type == tnurbs::ENGEOMETRYTYPE::NUNBS_CURVE)
    {
        tnurbs::nunbs3d* nurbs = static_cast<tnurbs::nunbs3d*>(curve);
        auto curves_segments = curve_split_by_plane<false>(nurbs, plane, toleranceData->m_data);
        for (const auto curve_seg : curves_segments)
        {
			std::string gid = create_uuid();
			m_parent->m_curves[gid] = curve_seg;
			tnurbs::curve_mesh_helper<tnurbs::nunbs3d> mh;
			tnurbs::disc_curve(curve_seg, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
			tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
			tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
			m_parent->m_curves_disc[gid] = curve_mesh;
			m_parent->m_nid_to_gid[m_id].push_back(gid);
			GeometryData* data = new GeometryData(gid);
			data->m_geoemtryType = 1;
			result[0].push_back(data);
			m_parent->m_node_data.push_back(data);
        }
    }
    else
    {
        assert(false);
    }

    for (auto& elem : result)
    {
        elem.erase(elem.begin() + 1, elem.end());
    }
    return result;
}

Loft::Loft(const std::string& id)
{
    size_t count = 1;
    m_id = id;
    m_inputs.resize(count);
    m_outputs.resize(1);
    for (size_t index = 0; index < count; ++index)
    {
        auto& node_interface = m_inputs[index];
        node_interface.set_kind_data(KindData::AsList);
    }

	auto& out_node_interface = m_outputs[0];
	out_node_interface.set_kind_data(KindData::AsItem);
}

std::vector<std::vector<AbstractNodeData*>> Loft::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<AbstractNodeData*> section_curves_data = get_data_list(0, node_iter);
    size_t count = section_curves_data.size();
    std::vector<tnurbs::nurbs3d*> section_curves(count);
    for (size_t index = 0; index < count; ++index)
    {
        AbstractNodeData* data = section_curves_data[index];
        if (data->m_type != NodeDataType::geometryData)
        {
            assert(false);
        }
        GeometryData* geoData = static_cast<GeometryData*>(data);
        if (geoData->m_geoemtryType != 1)
        {
            assert(false);
        }
        tnurbs::curve3d* curve = m_parent->m_curves[geoData->m_data];
        if (curve->get_type() != tnurbs::ENGEOMETRYTYPE::NURBS_CURVE)
        {
            assert(false);
        }
        tnurbs::nurbs3d* nurbs = static_cast<tnurbs::nurbs3d*>(curve);
        section_curves[index] = nurbs;
    }
    tnurbs::nurbs_surf3d* skin_surf = new tnurbs::nurbs_surf3d();
    tnurbs::skin_surface<double, 3, true>(2, section_curves, *skin_surf);

	std::string gid = create_uuid();
	m_parent->m_surfaces[gid] = skin_surf;
	tnurbs::surface_mesh_helper<tnurbs::nurbs_surf3d> mh;
	tnurbs::disc_surface(skin_surf, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
	tnurbs::mesh<2>* surf_mesh = new tnurbs::mesh<2>();
	tnurbs::mesh_help_to_mesh(mh, *surf_mesh);
	m_parent->m_surface_disc[gid] = surf_mesh;
	m_parent->m_nid_to_gid[m_id].push_back(gid);
	GeometryData* data = new GeometryData(gid);
	// std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid);
    data->m_geoemtryType = 2;
    std::vector<std::vector<AbstractNodeData*>> result{ {data}};
    m_parent->m_node_data.push_back(data);
    return result;
}


Flatten::Flatten(const std::string& id)
{
    m_id = id;
    m_inputs.resize(1);
    m_outputs.resize(1);
    for (size_t index = 0; index < 1; ++index)
    {
        auto& node_interface = m_inputs[index];
        // node_interface.set_kind_data(KindData::AsList);
    }

	auto& out_node_interface = m_outputs[0];
	// out_node_interface.set_kind_data(KindData::AsItem);
}

bool Flatten::calculate()
{
    if (true == is_calculate)
    {
        return true;
    }
	const auto& input_interface = m_inputs[0];
	auto& output_interface = m_outputs[0];
    Path path({ 0 });
    for (auto& element : input_interface.m_data)
    {
        auto& out_parma = output_interface.m_data[path];
        out_parma.insert(out_parma.end(), element.second.begin(), element.second.end());
    }
    return false;
}

ListItem::ListItem(const std::string& id, size_t count)
{
    m_id = id;
    m_outputs.resize(count);
    m_inputs.resize(3);
    for (size_t index = 0; index < count; ++index)
    {
        auto& node_interface = m_outputs[index];
        node_interface.set_kind_data(KindData::AsList);
    }

    for (size_t index = 0; index < 2; ++index)
    {
        auto& node_interface = m_inputs[index];
        node_interface.set_kind_data(KindData::AsList);
    }
	auto& node_interface = m_inputs[2];
	node_interface.set_kind_data(KindData::AsItem);
}

std::vector<std::vector<AbstractNodeData*>> ListItem::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<AbstractNodeData*> data = get_data_list(0, node_iter);
    std::vector<AbstractNodeData*> indexsData = get_data_list(1, node_iter);
    size_t index_count = indexsData.size();
    std::vector<int> indexs(index_count);
    for (size_t i = 0; i < index_count; ++i)
    {
        if (indexsData[i]->m_type != NodeDataType::numberData)
        {
            assert(false);
        }
        NumberData* index = static_cast<NumberData*>(indexsData[i]);
        indexs.push_back(static_cast<int> (index->m_data));
    }
    // wrap == true暂时写死
    bool wrap = true;

    if (data.size() < index_count)
    {
        data.insert(data.end(), index_count - data.size(), data.back());
    }
    size_t data_count = data.size();
    size_t output_count = m_outputs.size();
    std::vector<std::vector<AbstractNodeData*>> result(output_count);
    for (size_t j = 0; j < output_count; ++j)
    {
        std::vector<AbstractNodeData*> T(index_count);
        for (size_t k = 0; k < index_count; ++k)
        {
            size_t l = (k + j) % data_count;
            T[k] = data[l];
        }
        result[j] = T;
    }

    return result;
}

Line::Line(const std::string& id)
{
    m_id = id;
    m_outputs.resize(1);
    m_inputs.resize(2);
    for (size_t index = 0; index < 2; ++index)
    {
        auto& node_interface = m_inputs[index];
        node_interface.set_kind_data(KindData::AsItem);
    }
	auto& node_interface = m_outputs[0];
	node_interface.set_kind_data(KindData::AsItem);
}

std::vector<std::vector<AbstractNodeData*>> Line::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result{ {nullptr} };
    Eigen::Matrix3Xd control_points(3, 2);

    for (size_t index = 0; index < 2; ++index)
    {
		AbstractNodeData* abstractData = get_data_item(index, node_iter);
        if (abstractData->m_type != NodeDataType::pointData)
        {
            assert(false);
        }
        PointData* point = static_cast<PointData*>(abstractData);
        control_points.col(index) = point->m_data;
    }
    tnurbs::nunbs3d* curve = new tnurbs::nunbs3d(control_points);

	std::string gid = create_uuid();
	m_parent->m_curves[gid] = curve;
	tnurbs::curve_mesh_helper<tnurbs::nunbs3d> mh;
	tnurbs::disc_curve(curve, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
	tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
	tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
	m_parent->m_curves_disc[gid] = curve_mesh;
	GeometryData* data = new GeometryData(gid);
	// std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid);
    data->m_geoemtryType = 1;
    m_parent->m_nid_to_gid[m_id].push_back(gid);
    result[0][0] = data;
    m_parent->m_node_data.push_back(data);
    return result;
}

static void save_obj3(const tnurbs::nurbs_surface<double, 3, -1, -1, -1, -1, true>& surf, const char* path)
{
    tnurbs::surface_mesh_helper<tnurbs::nurbs_surface<double, 3, -1, -1, -1, -1, true>> mh;
    disc_surface(&surf, mh, 1e-4, 10.0, 0.1, 1.0);
    tnurbs:: mesh<2> surface_mesh;
    mesh_help_to_mesh2(mh, surface_mesh);

    save_obj(&surface_mesh, path);
    return;
}
std::vector<std::vector<AbstractNodeData*>> Sweep::calculate_inner(const NodeIterator& node_iter)
{
    AbstractNodeData* abstractData = get_data_item(0, node_iter);
    if (abstractData->m_type != NodeDataType::geometryData)
    {
        assert(false);
    }
    GeometryData* geoData = static_cast<GeometryData*>(abstractData);
    if (geoData->m_geoemtryType != 1)
    {
        assert(false);
    }
    tnurbs::curve3d* curve = m_parent->m_curves[geoData->m_data];
    tnurbs::nurbs3d* spine = static_cast<tnurbs::nurbs3d*>(curve);

    abstractData = get_data_item(1, node_iter);
    if (abstractData->m_type != NodeDataType::geometryData)
    {
        assert(false);
    }
    geoData = static_cast<GeometryData*>(abstractData);
    if (geoData->m_geoemtryType != 1)
    {
        assert(false);
    }
    curve = m_parent->m_curves[geoData->m_data];
    tnurbs::nunbs3d* profile = static_cast<tnurbs::nunbs3d*>(curve);
    
    tnurbs::nurbs_surf3d* sweep = new tnurbs::nurbs_surf3d();
    tnurbs::sweep_surface(*profile, *spine, 20, *sweep);
    volatile tnurbs::nunbs3d* cc1 = profile;
    volatile tnurbs::nurbs3d* cc2 = spine;

	std::string gid = create_uuid();
	m_parent->m_surfaces[gid] = sweep;
	tnurbs::surface_mesh_helper<tnurbs::nurbs_surf3d> mh;
	tnurbs::disc_surface(sweep, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
	tnurbs::mesh<2>* surf_mesh = new tnurbs::mesh<2>();
	tnurbs::mesh_help_to_mesh(mh, *surf_mesh);
	m_parent->m_surface_disc[gid] = surf_mesh;
	m_parent->m_nid_to_gid[m_id].push_back(gid);
	GeometryData* data = new GeometryData(gid);
	// std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid);
    data->m_geoemtryType = 2;
    std::vector<std::vector<AbstractNodeData*>> result{ {data} };
    m_parent->m_node_data.push_back(data);
	
    return result;
}

Sweep::Sweep(const std::string& id)
{
    m_id = id;
    m_inputs.resize(2);
    m_outputs.resize(1);
	auto& node_interface = m_inputs[0];
	node_interface.set_kind_data(KindData::AsItem);
	auto& node_interface2 = m_inputs[1];
	node_interface2.set_kind_data(KindData::AsItem);
	auto& node_interface3 = m_outputs[0];
	node_interface3.set_kind_data(KindData::AsItem);
}
