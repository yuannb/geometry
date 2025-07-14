#include "serialization.hpp"
#include <numberNode.hpp>
#include <InterpolateCurve.hpp>
#include <pointNode.hpp>
#include "MergeNode.hpp"

#include <charconv>
using json = nlohmann::ordered_json;
AbstractNodeData* serialize_utils::to_nodedata(json node_data)
{
	if (node_data.contains(ParamterType::type) == false)
	{
		return nullptr;
	}
	std::string data_type = node_data[ParamterType::type].get<std::string>();
	if (data_type == NodeDataType::numberData)
	{
		if (node_data.contains(ParamterType::data) == false)
		{
			return nullptr;
		}
		double value = node_data[ParamterType::data].get<double>();
		NumberData* data = new NumberData(value);
		return data;
	}
	else if (data_type == NodeDataType::pointData)
	{
		if (node_data.contains(ParamterType::data) == false)
		{
			return nullptr;
		}
		json vec_data = node_data[ParamterType::data];
		Eigen::Vector3d vec;
		for (int i = 0; i < 3; ++i)
		{
			vec[i] = vec_data[i].get<double>();
		}
		PointData* data = new PointData(vec);
		return data;
	}
	else if (data_type == NodeDataType::geometryData)
	{
		if (node_data.contains(ParamterType::data) == false)
		{
			return nullptr;
		}
		if (node_data.contains(ParamterType::geoemtryType) == false)
		{
			return nullptr;
		}
		std::string gid = node_data[ParamterType::data].get<std::string>();
		int geoType = node_data[ParamterType::geoemtryType].get<int>();
		GeometryData* data = new GeometryData(gid);
		data->m_geoemtryType = geoType;
		return data;
	}
	else
	{
		std::cout << "error" << std::endl;
	}
	return nullptr;
}

bool serialize_utils::to_path(const std::string& s, char delimiter, Path& path)
{
	std::string token;
	std::istringstream tokenStream(std::string(s.begin() + 1, s.end() - 1));

	// °´ delimiter ·Ö¸î×Ö·û´®
	path.m_path.clear();
	while (std::getline(tokenStream, token, delimiter))
	{
		size_t num;
		auto result = std::from_chars(token.data() + 1, token.data() + token.size() - 1, num);
		if (result.ec != std::errc{})
		{
			return false;
		}
		path.m_path.push_back(num);
	}

	return true;
}

bool serialize_utils::to_intface(json intf_json, NodeInterface& intf)
{
	if (intf_json.contains(ParamterType::data) == false)
	{
		return false;
	}
	json paramter_data = intf_json[ParamterType::data];

	for (auto& item : paramter_data.items())
	{
		Path path;
		to_path(item.key(), ',', path);
		std::vector<AbstractNodeData*> nodes;
		for (auto& value : item.value())
		{
			AbstractNodeData* data = to_nodedata(value);
			nodes.push_back(data);
		}
		intf.m_data.insert(std::make_pair(path, nodes));
	}

	return true;
}

