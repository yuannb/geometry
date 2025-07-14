#include "ReadJson.hpp"
#include "nodeFactory.hpp"
#include "databaseManager.hpp"
#include <functional>

using json = nlohmann::ordered_json;

static bool fillParamter(const json parameterJs, NodeInterface& parameter, std::map<std::string, NodeInterface*>& id2param, bool is_input)
{
	std::string parameter_id = parameterJs.at(ParamterType::id).template get<std::string>();
	parameter.m_id = parameter_id;
	id2param[parameter_id] = &parameter;
    std::string dataType = parameterJs.at(ParamterType::dataType).get<std::string>();
    if (is_input == true)
    {
		if (dataType == NodeDataType::numberData)
		{
			const json param_data = parameterJs.at(ParamterType::value);
			if (param_data.is_number() == false)
			{
				assert(false);
			}
			Path path({ 0 });
			NumberData* number = new NumberData(param_data.template get<double>());
			parameter.m_data[path] = { number };
		}
		else if (dataType == NodeDataType::pointData)
		{
			const json param_data = parameterJs.at(ParamterType::value);
			Path path({ 0 });
			if (param_data.is_array() == false || param_data.size() != 3)
			{
				assert(false);
			}
			Eigen::Vector3d point;
			for (size_t index = 0; index < 3; ++index)
			{
				point[index] = param_data[index].template get<double>();
			}
			PointData* point_data = new PointData(point);
			parameter.m_data[path] = { point_data };
		}
    }
    return true;
}

static bool fillNode(const json nodeJs, Node* node, const std::string& nodeId, std::map<std::string, NodeInterface*>& id2param, NodeGraph& ng)
{
    node->m_id = nodeId;
	const json node_inputs = nodeJs.at(ParamterType::inputs);

    size_t input_size = node_inputs.size();
    size_t actual_size = node->m_inputs.size();
    if (input_size > actual_size)
    {
        node->m_inputs.resize(input_size);
        for (size_t index = actual_size; index < input_size; ++index)
        {
            if (actual_size == 0)
            {
				node->m_inputs[index].m_kind = KindData::AsItem;
                continue;
            }
            node->m_inputs[index].m_kind = node->m_inputs[actual_size - 1].m_kind;
        }
    }

    int index = 0;
    for (const auto [name, value] : node_inputs.items())
    {
		fillParamter(value, node->m_inputs[index], id2param, true);
		index += 1;
    }
    
	const json node_outputs = nodeJs.at(ParamterType::outputs);
    size_t output_size = node_outputs.size();
    actual_size = node->m_outputs .size();
    if (output_size > actual_size)
    {
        node->m_outputs.resize(output_size);
        for (size_t index = actual_size; index < output_size; ++index)
        {
            if (actual_size == 0)
            {
				node->m_outputs[index].m_kind = KindData::AsItem;
                continue;
            }
            node->m_outputs[index].m_kind = node->m_outputs[actual_size - 1].m_kind;
        }
    }
    index = 0;
	for (const auto [name, value] : node_outputs.items())
	{
		fillParamter(value, node->m_outputs[index], id2param, false);
		++index;
	}

	ng.m_nodes[nodeId] = node;
	node->m_parent = &ng;
    return true;
}

#pragma optimize("", off)
bool CreateNodeGraph(const nlohmann::ordered_json& ng_data, NodeGraph& ng)
{
    if (ng_data.contains(ParamterType::nodes) == false)
    {
        return false;
    }
    const json ng_nodes = ng_data.at(ParamterType::nodes);

    if (ng_data.contains(ParamterType::order) == false)
    {
        return false;
    }
    const json ng_orders = ng_data.at(ParamterType::order);

    if (ng_data.contains(ParamterType::id) == false)
    {
        return false;
    }
    const json ng_id = ng_data.at(ParamterType::id);
    if (ng_id.is_string() == false)
    {
        return false;
    }
    ng.m_id = ng_id.get<std::string>();
    std::map<std::string, NodeInterface*> id2param;

    NodeGraphDBManager ngdb;
    bool is_exist = ngdb.is_ng_exits(ng.m_id);
    std::string ng_date = "";
    bool date_is_new{ false };
    if (is_exist == true)
    {
        ng_date = ngdb.get_ng_date(ng.m_id);
        // date_is_new = true;
        date_is_new = ng_date == ng_data.at("date").get<std::string>() ? true : false;
    }
    std::shared_ptr<spdlog::logger> logger = LoggerManager::getLogger();
	
    if (ng_data.contains(ParamterType::changed_nodes) == false)
	{
		logger->error("json of paramter: graphJson dont cotain the key of changed_nodes");
		return false;
	}
	// std::set<std::string> changed_nodes;
	for (const auto& element : ng_data.at(ParamterType::changed_nodes))
	{
		ng.changed_nodes.insert(element.get<std::string>());
    }

    if (ng_data.contains(ParamterType::depended_nodes) == false)
	{
		logger->error("json of paramter: graphJson dont cotain the key of depended_nodes");
		return false;
	}
	// std::set<std::string> depended_nodes;
	for (const auto& element : ng_data.at(ParamterType::depended_nodes))
	{
		ng.depended_nodes.insert(element.get<std::string>());
    }
    
    if (ng_data.contains(ParamterType::deleted_nodes) == false)
	{
		logger->error("json of paramter: graphJson dont cotain the key of status");
		return false;
	}
	// std::set<std::string> deleted_nodes;
	for (const auto& element : ng_data.at(ParamterType::deleted_nodes))
	{
		ng.deleted_nodes.insert(element.get<std::string>());
    }

    for (const json node_js : ng_nodes)
    {
        if (node_js.contains(ParamterType::id) == false)
        {
            return false;
        }
        const json node_id = node_js.at(ParamterType::id);
        
        if (node_id.is_string() == false)
        {
            return false;
        }

        bool idInDependedNodes = ng.depended_nodes.find(node_id.get<std::string>()) == ng.depended_nodes.end() ? false : true;
        bool idInChangedNodes = ng.changed_nodes.find(node_id.get<std::string>()) == ng.changed_nodes.end() ? false : true;
        ng.from_database = date_is_new;
        if (idInDependedNodes && date_is_new == true)
        {
            bool flag = ngdb.CreateNodeFromDB(node_js, ng, id2param);
            if (flag == false)
            {
                return false;
            }
        }
        else if (date_is_new && idInChangedNodes == false)
        {
            continue;
        }
        else
        {
            const json node_type = node_js.at(ParamterType::type);
            if (node_type.is_string() == false)
            {
                return false;
            }

            std::string type = node_type.template get<std::string>();
            std::string id = node_id.template get<std::string>();
            Node* node = NodeFactory::createObject(type);
            fillNode(node_js, node, id, id2param, ng);
            if (node == nullptr)
            {
                std::cout << "CreateNodeGraphFromJson : do not have this type node: " << type << std::endl;
                return false;
            }

        }

    }

    std::map<std::string, std::string> interface_relation;
    const json connects = ng_data.at("connections");
    for ( auto it = connects.begin(); it != connects.end(); ++it)
    {
        const json& current = *it;
        std::string from_id = current.at("from_id").template get<std::string>();
        std::string to_id = current.at("to_id").template get<std::string>();
        interface_relation[to_id] = from_id;
    }

    for (auto it = ng.m_nodes.begin(); it != ng.m_nodes.end(); ++it)
    {
        size_t inputs_count = it->second->m_inputs.size();
        for (size_t index = 0; index < inputs_count; ++index)
        {
            it->second->m_inputs[index].m_from_id = interface_relation[it->second->m_inputs[index].m_id];
            it->second->m_inputs[index].m_from_interface = id2param[interface_relation[it->second->m_inputs[index].m_id]];
        }
    }

    const json orders = ng_data.at("order");
    for (auto it = orders.begin(); it != orders.end(); ++it)
    {
        std::string id = it->template get<std::string>();
        ng.m_nodes_sort.push_back(id);
    }


    return true;
}
#pragma optimize("", on)
