#include "pointNode.hpp"
#include "NodeGraph.hpp"


ConstructPoint::ConstructPoint(const std::string &id)
{
    m_id = id;
    m_inputs.resize(3);
    m_outputs.resize(1);
    for (int index = 0; index < 3; ++index)
    {
        auto& interface = m_inputs[index];
        interface.set_kind_data(KindData::AsItem);
    }
    NodeInterface *interface = new NodeInterface();
    interface->set_kind_data(KindData::AsItem);
}

std::vector<std::vector<AbstractNodeData*>> ConstructPoint::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result{ {nullptr} };
    std::vector<double> coords(3);
    for (size_t index = 0; index < 3; ++index)
    {
        auto iter = node_iter.m_node_iter[index].m_iter;
        AbstractNodeData* coord = iter.m_path_iter->second[iter.m_index];
        if (coord->m_type != NodeDataType::numberData)
        {
            assert(false);
        }
        auto number_data = static_cast<NumberData*>(coord);
        coords[index] = number_data->m_data;
    }
    PointData* point = new PointData(coords[0], coords[1], coords[2]);
    // std::unique_ptr<PointData> point = std::make_unique<PointData>(coords[0], coords[1], coords[2]);

	std::string gid = create_uuid();
    Eigen::Vector3d p(coords[0], coords[1], coords[2]);
    m_parent->m_points_disc[gid] = p;
	m_parent->m_nid_to_gid[m_id].push_back(gid);
    result[0][0] = point;
    m_parent->m_node_data.push_back(point);
    return result;
}

DistancePt::DistancePt(const std::string& id)
{
    m_id = id;
    size_t inputCount = 2;
    size_t outputCount = 1;
    m_inputs.resize(inputCount);
    m_outputs.resize(outputCount);
    for (int index = 0; index < inputCount; ++index)
    {
        auto& interface = m_inputs[index];
        interface.set_kind_data(KindData::AsItem);
    }
    for (int index = 0; index < outputCount; ++index)
    {
        auto& interface = m_outputs[index];
        interface.set_kind_data(KindData::AsItem);
    }
}

std::vector<std::vector<AbstractNodeData*>> DistancePt::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result{ {nullptr} };
    std::array<Eigen::Vector3d, 2> points;
    for (size_t index = 0; index < 2; ++index)
    {
		AbstractNodeData* data = get_data_item(index, node_iter);
		if (data->m_type != NodeDataType::pointData)
		{
			assert(false);
		}
        PointData* pointData = static_cast<PointData*>(data);
        points[index] = pointData->m_data;
    }
    double distance = (points[1] - points[0]).norm();

    NumberData* distData = new NumberData(distance);
    // std::unique_ptr<NumberData> distData = std::make_unique<NumberData>(distance);
    result[0][0] = distData;
    m_parent->m_node_data.push_back(distData);
    return result;
}

EndsPoints::EndsPoints(const std::string& id)
{
    m_id = id;
    size_t inputCount = 1;
    size_t outputCount = 4;
    m_inputs.resize(inputCount);
    m_outputs.resize(outputCount);
    for (int index = 0; index < inputCount; ++index)
    {
        auto& interface = m_inputs[index];
        interface.set_kind_data(KindData::AsItem);
    }
    for (int index = 0; index < outputCount; ++index)
    {
        auto& interface = m_outputs[index];
        interface.set_kind_data(KindData::AsItem);
    }
}

std::vector<std::vector<AbstractNodeData*>> EndsPoints::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>>result(4);

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

    std::array<Eigen::Vector3d, 2> endpoints;
    curve->get_ends_point(endpoints);
    tnurbs::Box<double, 1> domain = curve->get_domain_box();

    for (size_t index = 0; index < 2; ++index)
    {
        PointData* pointData = new PointData(endpoints[index]);
        // std::unique_ptr<PointData> pointData = std::make_unique<PointData>(endpoints[index]);
		std::string gid = create_uuid();
		m_parent->m_points_disc[gid] = endpoints[index];
		m_parent->m_nid_to_gid[m_id].push_back(gid);
        result[index].push_back(pointData);
        m_parent->m_node_data.push_back(pointData);
    }
    std::array<double, 2> endParam{ domain.Min[0], domain.Max[0] };
    for (size_t index = 0; index < 2; ++index)
    {
        NumberData* numberData = new NumberData(endParam[index]);
        // std::unique_ptr<NumberData> numberData = std::make_unique<NumberData>(endParam[index]);
        result[index + 2].push_back(numberData);
        m_parent->m_node_data.push_back(numberData);
    }
    return result;
}

PrepFrame::PrepFrame(const std::string& id)
{
    m_id = id;
    size_t inputCount = 2;
    size_t outputCount = 3;
    m_inputs.resize(inputCount);
    m_outputs.resize(outputCount);
    for (int index = 0; index < inputCount; ++index)
    {
        auto& interface = m_inputs[index];
        interface.set_kind_data(KindData::AsItem);
    }
    for (int index = 0; index < outputCount; ++index)
    {
        auto& interface = m_outputs[index];
        interface.set_kind_data(KindData::AsItem);
    }
}

std::vector<std::vector<AbstractNodeData*>> PrepFrame::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result(3);
    AbstractNodeData* data = get_data_item(0, node_iter);
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
    
    data = get_data_item(1, node_iter);
    if (data->m_type != NodeDataType::numberData)
    {
        assert(false);
    }
    NumberData* param = static_cast<NumberData*>(data);
    Eigen::VectorX<Eigen::Vector3d> ders;
    curve->derivative_on_curve(param->m_data, 2, ders);
   
    Eigen::Vector3d tangent = ders[1];
    tangent.normalize();
    if (tangent.squaredNorm() < 0.1)
    {
        assert(false);
    }
    
    Eigen::Vector3d noraml = ders[2];
    noraml.normalize();
    noraml -= noraml.dot(tangent) * tangent;
    noraml.normalize();

    Eigen::Vector3d binormal = tangent.cross(noraml);

    if (noraml.squaredNorm() < 0.1)
    {
        assert(false);
    }
   
    PointData* pt = new PointData(ders[0]);
    PointData* normal_vec = new PointData(noraml);
    PointData* binormal_vec = new PointData(binormal);

    result[0].push_back(pt);
    result[1].push_back(normal_vec);
    result[2].push_back(binormal_vec);

    m_parent->m_node_data.push_back(pt);
    m_parent->m_node_data.push_back(normal_vec);
    m_parent->m_node_data.push_back(binormal_vec);

    // result[0].push_back(new PointData(ders[0]));
    // result[1].push_back(new PointData(noraml));
    // result[2].push_back(new PointData(binormal));
    return result;
}

static tnurbs::nunbs_surf3d* create_plane(const Eigen::Vector3d& origin, const Eigen::Vector3d& udir, const Eigen::Vector3d& vdir, double uLength, double vLength)
{
	std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(2);
	control_points[0].resize(3, 2);
	control_points[1].resize(3, 2);
    double uhalfLength = uLength / 2.0;
    double vhalfLength = vLength / 2.0;
	control_points[0].col(0) = origin - uhalfLength * udir - vhalfLength * vdir;
	control_points[0].col(1) = origin + uhalfLength * udir - vhalfLength * vdir;

	control_points[1].col(0) = origin - uhalfLength * udir + vhalfLength * vdir;
	control_points[1].col(1) = origin + uhalfLength * udir + vhalfLength * vdir;
	Eigen::VectorXd knots_vector(4);
	knots_vector[0] = knots_vector[1] = 0.0;
	knots_vector[2] = knots_vector[3] = 1.0;

	tnurbs::nunbs_surf3d* plane = new  tnurbs::nunbs_surf3d(knots_vector, knots_vector, control_points);
    return plane;

}

static tnurbs::nunbs3d* create_rectangle(const Eigen::Vector3d& origin, const Eigen::Vector3d& udir, const Eigen::Vector3d& vdir, double uLength, double vLength)
{
	Eigen::Matrix<double, 3, Eigen::Dynamic> control_points(3, 5);;
    double uhalfLength = uLength / 2.0;
    double vhalfLength = vLength / 2.0;
	control_points.col(0) = origin - uhalfLength * udir - vhalfLength * vdir;
	control_points.col(1) = origin + uhalfLength * udir - vhalfLength * vdir;
	control_points.col(2) = origin + uhalfLength * udir + vhalfLength * vdir;
	control_points.col(3) = origin - uhalfLength * udir + vhalfLength * vdir;
	control_points.col(4) = control_points.col(0);
	Eigen::VectorXd knots_vector(7);
	knots_vector[0] = knots_vector[1] = 0.0;
    knots_vector[2] = 0.25;
    knots_vector[3] = 0.5;
    knots_vector[4] = 0.75;
	knots_vector[5] = knots_vector[6] = 1.0;

	tnurbs::nunbs3d* rectangle = new  tnurbs::nunbs3d(knots_vector, control_points);
    return rectangle;
}
Rectangle::Rectangle(const std::string& id)
{
    m_id = id;
    size_t inputCount = 5;
    size_t outputCount = 1;
    m_inputs.resize(inputCount);
    m_outputs.resize(outputCount);
    for (int index = 0; index < inputCount; ++index)
    {
        auto& interface = m_inputs[index];
        interface.set_kind_data(KindData::AsItem);
    }
    for (int index = 0; index < outputCount; ++index)
    {
        auto& interface = m_outputs[index];
        interface.set_kind_data(KindData::AsItem);
    }
}

std::vector<std::vector<AbstractNodeData*>> Rectangle::calculate_inner(const NodeIterator& node_iter)
{
    std::vector<std::vector<AbstractNodeData*>> result{ {nullptr} };
    std::array<Eigen::Vector3d, 3> vecs;
    for (size_t index = 0; index < 3; ++index)
    {
        AbstractNodeData* data = get_data_item(index, node_iter);
        if (data->m_type != NodeDataType::pointData)
        {
            assert(false);
        }
        PointData* vec = static_cast<PointData*>(data);
        vecs[index] = vec->m_data;
    }
    std::array<double, 2> lengths;

    for (size_t index = 0; index < 2; ++index)
    {
        AbstractNodeData* data = get_data_item(index + 3, node_iter);
        if (data->m_type != NodeDataType::numberData)
        {
            assert(false);
        }
        NumberData* length = static_cast<NumberData*>(data);
        lengths[index] = length->m_data;
    }
    tnurbs::nunbs3d* rec = create_rectangle(vecs[0], vecs[1], vecs[2], lengths[0], lengths[1]);
	std::string gid = create_uuid();
	m_parent->m_curves[gid] = rec;
	tnurbs::curve_mesh_helper<tnurbs::nunbs3d> mh;
	tnurbs::disc_curve(rec, mh, tnurbs::TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
	tnurbs::mesh<1>* curve_mesh = new tnurbs::mesh<1>();
	tnurbs::mesh_help_to_mesh(mh, *curve_mesh);
	m_parent->m_curves_disc[gid] = curve_mesh;
	GeometryData* data = new GeometryData(gid);
    // std::unique_ptr<GeometryData> data = std::make_unique<GeometryData>(gid);
    data->m_geoemtryType = 1;
    result[0][0] = data;
	m_parent->m_nid_to_gid[m_id].push_back(gid);
    m_parent->m_node_data.push_back(data);
    return result;
}
