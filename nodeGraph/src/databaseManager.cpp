#include "databaseManager.hpp"


using json = nlohmann::ordered_json;

void NodeGraphDBManager::check_error(int status)
{
	if (status != 0) {
		logger->error("NodeGraphDBManager: check_error: status != 0");
		throw DatabaseException(mysql_error(conn));
	}
}

void NodeGraphDBManager::escape_string(std::string& str)
{
	char* escaped = new char[str.length() * 2 + 1];
	mysql_real_escape_string(conn, escaped, str.c_str(), str.length());
	str = escaped;
	delete[] escaped;
}

NodeGraphDBManager::NodeGraphDBManager() {
	conn = mysql_init(nullptr);
	if (conn == nullptr)
	{
		logger->error("NodeGraphDBManager: construct function: Failed to initialize MySQL connection");
		throw DatabaseException("Failed to initialize MySQL connection");
	}

	if (nullptr == mysql_real_connect(conn, host, user, password, dbname, port, nullptr, 0))
	{
		logger->error("NodeGraphDBManager: construct function: Failed to connect to database, ERROR: {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to connect to database: ") + mysql_error(conn));
	}

	// 设置字符集
	if (mysql_set_character_set(conn, "utf8")) {
		logger->error("NodeGraphDBManager: construct function: Failed to set character set, ERROR: {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to set character set: ") + mysql_error(conn));
	}
}

NodeGraphDBManager::~NodeGraphDBManager() {
	if (conn) {
		std::cout << "Closing conn: " << mysql_ping(conn) << std::endl;
		mysql_close(conn);
		conn = nullptr;
	}
}

    // 创建必要的表
void NodeGraphDBManager::create_tables() {
	const char* queries[] = {
		"CREATE TABLE IF NOT EXISTS nodegraphs ("
		"  uuid VARCHAR(72) PRIMARY KEY,"
		"  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP"
		")",

		"CREATE TABLE IF NOT EXISTS nodes ("
		"  uuid VARCHAR(72) PRIMARY KEY,"
		"  type VARCHAR(72) NOT NULL,"
		"  input_count INT UNSIGNED NOT NULL,"
		"  output_count INT UNSIGNED NOT NULL,"
		"  nodegraph_uuid VARCHAR(72) NOT NULL,"
		"  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP"
		")",

		"CREATE TABLE IF NOT EXISTS parameters ("
		"  uuid VARCHAR(72) PRIMARY KEY,"
		"  is_input BOOLEAN DEFAULT false,"
		"  param_index INT UNSIGNED NOT NULL,"
		"  node_uuid VARCHAR(72) NOT NULL,"
		"  from_id VARCHAR(72) NOT NULL DEFAULT '',"
		"  to_id VARCHAR(72) NOT NULL DEFAULT '',"
		"  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP"
		")",

		"CREATE TABLE IF NOT EXISTS geometry ("
		"  uuid VARCHAR(72) PRIMARY KEY,"
		"  node_uuid VARCHAR(72) NOT NULL,"
		"  `value` TEXT,"
		"  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP"
		")",

		"CREATE TABLE IF NOT EXISTS parameter_values ("
		"  id INT AUTO_INCREMENT PRIMARY KEY,"
		"  parameter_uuid VARCHAR(72) NOT NULL,"
		"  `key` VARCHAR(255) NOT NULL,"
		"  `value` TEXT,"
		"  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP"
		")"
	};

        for (const auto& query : queries)
        {
            if (0 != mysql_query(conn, query))
            {
                logger->error("NodeGraphDBManager: create_tables function: Failed to create table : {}", mysql_error(conn));
                throw DatabaseException(std::string("Failed to create table: ") + mysql_error(conn));
            }
        }
    }

void NodeGraphDBManager::add_geometry(const std::string& node_uuid, const NodeGraph& ng)
{
	if (auto it = ng.m_nid_to_gid.find(node_uuid); it != ng.m_nid_to_gid.end())
	{
		auto& geo_id_vec = it->second;
		for (auto gid = geo_id_vec.begin(); gid != geo_id_vec.end(); ++gid)
		{
			std::string geo_string;
			if (auto pit = ng.m_points_disc.find(*gid); pit != ng.m_points_disc.end())
			{
				json point_json = point_to_json(pit->second);
				geo_string = point_json.dump();
			}
			else if (auto cit = ng.m_curves.find(*gid); cit != ng.m_curves.end())
			{
				const tnurbs::curve3d* curve = cit->second;
				json curve_json = curve_to_json(curve);
				geo_string = curve_json.dump();
			}
			else if (auto sit = ng.m_surfaces.find(*gid); sit != ng.m_surfaces.end())
			{
				const tnurbs::surf3d* surf = sit->second;
				json surf_json = surf_to_json(surf);
				geo_string = surf_json.dump();
			}
			else
			{
				logger->error("NodeGraphDBManager: add_geometry(): geometry is invalid");
				continue;
			}
			std::string query = "INSERT INTO geometry (uuid, node_uuid, value) VALUES ('"
				+ *gid + "', '" + it->first + "', '" + geo_string + "')";
			if (mysql_query(conn, query.c_str()))
			{
				logger->error("NodeGraphDBManager: add_geometry(): Failed to add geometry : {}", mysql_error(conn));
				throw DatabaseException(std::string("Failed to add geometry: ") + mysql_error(conn));
			}
		}
	}
	return;
}

    // 添加 NodeGraph
void NodeGraphDBManager::add_nodegraph(const NodeGraph& nodegraph) {
	std::string query = "INSERT INTO nodegraphs (uuid) VALUES ('" + nodegraph.m_id + "')";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: add_nodegraph function: Failed to add nodegraph : {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to add nodegraph: ") + mysql_error(conn));
	}

	// 保存所有节点
	for (const auto& node : nodegraph.m_nodes)
	{
		add_node(*(node.second), nodegraph.m_id);
	}
}

// 添加 Node
void NodeGraphDBManager::add_node(const Node& node, const std::string& nodegraph_uuid)
{
	std::string query = "INSERT INTO nodes (uuid, type, input_count, output_count, nodegraph_uuid) VALUES ('"
		+ node.m_id + "', '" + node.get_type() + "', '" + std::to_string(node.m_inputs.size()) + "', '" + std::to_string(node.m_outputs.size()) + "', '" + nodegraph_uuid + "')";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: add_node(): Failed to add node : {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to add node: ") + mysql_error(conn));
	}

	// 保存所有参数
	int param_index = 0;
	for (const auto& param : node.m_outputs)
	{
		add_parameter(param, false, param_index, node);
		++param_index;
	}
	add_geometry(node.m_id, *node.m_parent);
}
    // 添加 Node

    // 添加 Parameter
void NodeGraphDBManager::add_parameter(const NodeInterface& param, bool is_input, int index, const Node& node)
{
	std::string node_uuid = node.m_id;
	int str = is_input == true ? 1 : 0;

	json to_id = json(param.m_to_id);
	std::string to_id_str = to_id.dump();
	std::string query = "INSERT INTO parameters (uuid, is_input, param_index, node_uuid, from_id, to_id) VALUES ('"
		+ param.m_id + "', '" + std::to_string(str) + "', '" + std::to_string(index) + "', '" + node_uuid + "', '" + param.m_from_id + "', '" + to_id_str + "')";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: add_node(): Failed to add parameter : {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to add parameter: ") + mysql_error(conn));
	}

	// 保存所有键值对
	for (const auto& pair : param.m_data)
	{
		const Path& path = pair.first;
		const std::vector<AbstractNodeData*> node_dataes = pair.second;
		json valueVector = json::array();
		for (auto value : node_dataes)
		{
			valueVector.push_back(value->to_json());
		}
		std::string values = valueVector.dump();

		std::string pathstring = path.to_string();
		escape_string(pathstring);
		escape_string(values);

		std::string value_query = "INSERT INTO parameter_values (parameter_uuid, `key`, `value`) VALUES ('"
			+ param.m_id + "', '" + pathstring + "', '" + values + "')";
		if (mysql_query(conn, value_query.c_str()))
		{
			logger->error("NodeGraphDBManager: add_node(): Failed to add parameter_value : {}", mysql_error(conn));
			throw DatabaseException(std::string("Failed to add parameter value: ") + mysql_error(conn));
		}
	}
}

    // 查询 geometry
bool NodeGraphDBManager::get_geometry(const std::string& gid, NodeGraph& ng)
{
	// 查询所有键值对
	if (ng.m_curves.find(gid) != ng.m_curves.end())
	{
		return true;
	}
	if (ng.m_points_disc.find(gid) != ng.m_points_disc.end())
	{
		return true;
	}
	if (ng.m_surfaces.find(gid) != ng.m_surfaces.end())
	{
		return true;
	}
	std::string query = "SELECT `value` FROM geometry WHERE uuid = '" + gid + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: get_geometry(): Failed to query geometry {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to query geometry: ") + mysql_error(conn));
	}

	MYSQL_RES* result = mysql_store_result(conn);
	if (result)
	{
		MYSQL_ROW row;
		while ((row = mysql_fetch_row(result)))
		{
			std::string value = row[0] ? row[0] : "";
			json value_json = nlohmann::ordered_json::parse(value);
			if (value_json[ParamterType::geoemtryType].get<int>() == 0)
			{
				Eigen::Vector3d pt;
				json_to_point(value_json, pt);
				ng.m_points_disc.insert(std::make_pair(gid, pt));
			}
			else if (value_json[ParamterType::geoemtryType].get<int>() == 1)
			{
				tnurbs::curve3d* curve = json_to_curve(value_json);
				ng.m_curves.insert(std::make_pair(gid, curve));
			}
			else if (value_json[ParamterType::geoemtryType].get<int>() == 2)
			{
				tnurbs::surf3d* surf = json_to_surf(value_json);
				ng.m_surfaces.insert(std::make_pair(gid, surf));
			}
			else
			{
				logger->error("NodeGraphDBManager: get_geometry(): invalid geoemty data");
				continue;
			}
		}
		mysql_free_result(result);
	}

	return true;
}

bool NodeGraphDBManager::is_ng_exits(const std::string& ng_id)
{
	std::string query = "SELECT uuid FROM nodegraphs WHERE uuid = '" + ng_id + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: is_ng_exits(): Failed to query nodegraph : {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to query nodegraph: ") + mysql_error(conn));
	}
	MYSQL_RES* result = mysql_store_result(conn);
	if (result == nullptr)
	{
		logger->error("NodeGraphDBManager: is_ng_exits(): NodeGraph not found : {}", mysql_error(conn));
	}
	if (mysql_num_rows(result) == 0)
	{
		mysql_free_result(result);
		return false;
	}
	mysql_free_result(result);
	return true;
}

bool NodeGraphDBManager::update_ng_time(const std::string& ng_id)
{
	std::string query = "UPDATE nodegraphs SET created_at = CURRENT_TIMESTAMP WHERE uuid = '" + ng_id + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: is_ng_exits(): Failed to query nodegraph : {}", mysql_error(conn));
	}
	return true;
}
std::string NodeGraphDBManager::get_ng_date(const std::string& ng_id)
{
	std::string query = "SELECT created_at FROM nodegraphs WHERE uuid = '" + ng_id + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: get_ng_date(): Failed to query nodegraph : {}", mysql_error(conn));
		return "-1";
	}
	MYSQL_RES* result = mysql_store_result(conn);
	if (result == nullptr)
	{
		logger->error("NodeGraphDBManager: get_ng_date(): NodeGraph not found : {}", mysql_error(conn));
		return "-1";
	}
	if (mysql_num_rows(result) == 0)
	{
		mysql_free_result(result);
		return "-1";
	}
	MYSQL_ROW row = mysql_fetch_row(result);
	if (row[0] == nullptr)
	{
		mysql_free_result(result);
		logger->error("NodeGraphDBManager: get_ng_date(): NodeGraph not found : {}", mysql_error(conn));
		return "-1";
	}
	std::string reulst = row[0];
	mysql_free_result(result);
	return reulst;
}

    // 查询 NodeGraph
NodeGraph NodeGraphDBManager::get_nodegraph(const std::string& uuid)
{
	NodeGraph nodegraph(uuid);
	std::map<std::string, NodeInterface*> id2param;
	// 查询 NodeGraph 是否存在
	std::string query = "SELECT uuid FROM nodegraphs WHERE uuid = '" + uuid + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: get_nodegraph(): Failed to query nodegraph : {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to query nodegraph: ") + mysql_error(conn));
	}

	MYSQL_RES* result = mysql_store_result(conn);
	if (result == nullptr || mysql_num_rows(result) == 0)
	{
		mysql_free_result(result);
		logger->error("NodeGraphDBManager: get_nodegraph(): NodeGraph not found : {}", mysql_error(conn));
		throw DatabaseException("NodeGraph not found");
	}
	mysql_free_result(result);

	// 查询所有节点
	query = "SELECT uuid, type, input_count, output_count FROM nodes WHERE nodegraph_uuid = '" + uuid + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: get_nodegraph(): Failed to query nodes : {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to query nodes: ") + mysql_error(conn));
	}

	result = mysql_store_result(conn);
	if (result)
	{
		MYSQL_ROW row;
		while ((row = mysql_fetch_row(result)))
		{
			std::string node_uuid = row[0];
			int input_count = std::stoi(row[2]);
			int output_count = std::stoi(row[3]);
			Node* node = get_node(node_uuid, nodegraph, input_count, output_count, id2param);
			nodegraph.m_nodes.insert(std::make_pair(node->m_id, node));
		}
		mysql_free_result(result);
	}

	return nodegraph;
}

    // 查询 Node
Node* NodeGraphDBManager::get_node(const std::string& uuid, NodeGraph& ng, int input_count, int output_count, std::map<std::string, NodeInterface*>& id2param)
{

	// 查询所有参数
	std::string query = "SELECT uuid, param_index, from_id, to_id FROM parameters WHERE node_uuid = '" + uuid + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: get_node(): Failed to query parameters {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to query parameters: ") + mysql_error(conn));
	}

	MYSQL_RES* result = mysql_store_result(conn);
	Node* node = new Node(uuid);
	node->is_calculate = true;
	node->m_inputs.resize(input_count);
	node->m_outputs.resize(output_count);
	node->m_parent = &ng;
	if (result)
	{
		MYSQL_ROW row;
		auto rowsCount = result->row_count;
		node->m_outputs.reserve(rowsCount);
		while ((row = mysql_fetch_row(result)))
		{
			std::string param_uuid = row[0];
			if (row[1] == nullptr)
			{
				logger->error("NodeGraphDBManager: get_node(): data is nullptr");
				continue;
			}
			int param_index = std::stoi(row[1]);
			NodeInterface& intf = node->m_outputs[param_index];
			get_parameter(param_uuid, intf, node);
			intf.m_from_id = row[2];
			id2param[param_uuid] = &intf;

			json to_id_js = nlohmann::ordered_json::parse(row[3]);
			size_t to_id_count = to_id_js.size();
			intf.m_to_id.resize(to_id_count);
			for (size_t i = 0; i < to_id_count; ++i)
			{
				intf.m_to_id[i] = to_id_js[i].get<std::string>();
			}
		}
	}

	mysql_free_result(result);
	return node;
}

    // 查询 Parameter
bool NodeGraphDBManager::get_parameter(const std::string& uuid, NodeInterface& param, Node* node)
{
	param.m_id = uuid;

	// 查询所有键值对
	std::string query = "SELECT `key`, `value` FROM parameter_values WHERE parameter_uuid = '" + uuid + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: get_parameter(): Failed to query parameter values {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to query parameter values: ") + mysql_error(conn));
	}

	MYSQL_RES* result = mysql_store_result(conn);
	if (result)
	{
		MYSQL_ROW row;
		while ((row = mysql_fetch_row(result)))
		{
			std::string key = row[0] ? row[0] : "";
			std::string value = row[1] ? row[1] : "";

			Path path;
			serialize_utils::to_path(key, ',', path);
			json value_json = nlohmann::ordered_json::parse(value);
			std::vector<AbstractNodeData*> dataes;
			for (auto& data : value_json)
			{
				AbstractNodeData* node_data = serialize_utils::to_nodedata(data);
				dataes.push_back(node_data);
				node->m_parent->m_node_data.push_back(node_data);
				if (node_data->m_type == NodeDataType::geometryData)
				{
					GeometryData* geo_data = static_cast<GeometryData*>(node_data);
					get_geometry(geo_data->m_data, *node->m_parent);
					// node->m_parent->m_nid_to_gid[node->m_id].push_back(geo_data->m_data);
				}
			}
			param.m_data.insert(std::make_pair(path, dataes));
		}
		mysql_free_result(result);
	}

	return true;
}

    // 查询 geometry
bool NodeGraphDBManager::delete_geometry(const std::string& nid)
{
	// 查询所有键值对
	std::string query = "DELETE FROM geometry WHERE node_uuid = '" + nid + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: delete_geometry(): Failed to query geometry {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to query geometry: ") + mysql_error(conn));
	}
	int rowsAffected = mysql_affected_rows(conn);
	logger->info("NodeGraphDBManager: delete_geometry(): delete {} lines records", rowsAffected);
	return true;
}

bool NodeGraphDBManager::find_nodes_id(const std::string& ng_id, std::vector<std::string>& nodes_id)
{
	// 查询所有参数
	std::string query = "SELECT uuid FROM nodes WHERE nodegraph_uuid = '" + ng_id + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: find_nodes_id(): Failed to find_nodes_id {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to find_nodes_id: ") + mysql_error(conn));
	}

	MYSQL_RES* result = mysql_store_result(conn);
	if (result)
	{
		MYSQL_ROW row;
		auto rowsCount = result->row_count;
		nodes_id.clear();
		nodes_id.reserve(rowsCount);
		for (uint64_t index = 0; index < rowsCount; ++index)
		{
			row = mysql_fetch_row(result);
			std::string param_uuid = row[0];
			nodes_id.push_back(param_uuid);
		}
		mysql_free_result(result);
	}

	return true;
}


    // 查询 NodeGraph
bool NodeGraphDBManager::delete_nodegraph(const std::string& uuid) {
	// 查询 NodeGraph 是否存在
	std::vector<std::string> nodes_id;
	find_nodes_id(uuid, nodes_id);

	for (const std::string& nid : nodes_id)
	{
		delete_node(nid);
		delete_geometry(nid);
	}

	std::string query = "DELETE FROM nodegraphs WHERE uuid = '" + uuid + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: delete_nodegraph(): Failed to delete_nodegraph : {}", mysql_error(conn));
	}
	uint64_t rowsAffected = mysql_affected_rows(conn);
	logger->info("NodeGraphDBManager: delete_node_graph(): delete {} lines records", rowsAffected);

	return true;
}

bool NodeGraphDBManager::find_paramter_id(const std::string& node_id, std::vector<std::string>& paramters_id)
{
	// 查询所有参数
	std::string query = "SELECT uuid FROM parameters WHERE node_uuid = '" + node_id + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: find_paramter_id(): Failed to find paramter id {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to find paramter id: ") + mysql_error(conn));
	}

	MYSQL_RES* result = mysql_store_result(conn);
	if (result)
	{
		MYSQL_ROW row;
		auto rowsCount = result->row_count;
		paramters_id.clear();
		paramters_id.reserve(rowsCount);
		for (uint64_t index = 0; index < rowsCount; ++index)
		{
			row = mysql_fetch_row(result);
			std::string param_uuid = row[0];
			paramters_id.push_back(param_uuid);
		}
		mysql_free_result(result);
	}

	return true;
}


    // 查询 Node
bool NodeGraphDBManager::delete_node(const std::string& node_id)
{
	std::vector<std::string> paramters_id;
	find_paramter_id(node_id, paramters_id);
	for (const std::string& pid : paramters_id)
	{
		delete_parameter(pid);
	}
	std::string query = "DELETE FROM nodes WHERE uuid = '" + node_id + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: delete_node(): Failed to delete node {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to delete node: ") + mysql_error(conn));
	}
	uint64_t rowsAffected = mysql_affected_rows(conn);
	logger->info("NodeGraphDBManager: delete_node(): delete {} lines records", rowsAffected);

	return true;
}

    // 查询 Parameter
bool NodeGraphDBManager::delete_parameter(const std::string& uuid)
{
	// 查询所有键值对
	std::string query = "DELETE FROM parameter_values WHERE parameter_uuid = '" + uuid + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: delete_parameter(): Failed to delete parameter_values {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to delete parameter_values: ") + mysql_error(conn));
	}
	int rowsAffected = mysql_affected_rows(conn);
	logger->info("NodeGraphDBManager: delete_parameter_values(): delete {} lines records", rowsAffected);

	query = "DELETE FROM parameters WHERE uuid = '" + uuid + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: delete_parameter(): Failed to delete parameter {}", mysql_error(conn));
		throw DatabaseException(std::string("Failed to delete parameters: ") + mysql_error(conn));
	}

	rowsAffected = mysql_affected_rows(conn);
	logger->info("NodeGraphDBManager: delete_parameter(): delete {} lines records", rowsAffected);
	return true;
}


bool NodeGraphDBManager::CreateNodeFromDB(const json node_js, NodeGraph& ng, std::map<std::string, NodeInterface*>& id2param)
{
	if (node_js.contains(ParamterType::id) == false)
	{
		logger->error("the json of node dont contain the key of order");
		return false;
	}
	const json node_id = node_js.at(ParamterType::id);
	if (node_id.is_string() == false)
	{
		logger->error("the json of node_type is not string");
		return false;
	}
	std::string id = node_id.template get<std::string>();
	std::string query = "SELECT type, input_count, output_count FROM nodes WHERE uuid = '" + id + "'";
	if (mysql_query(conn, query.c_str()))
	{
		logger->error("NodeGraphDBManager: CreateNodeGraphFromDataBase(): Failed to query nodes : {}", mysql_error(conn));
		return false;
	}

	MYSQL_RES* result = mysql_store_result(conn);
	if (result)
	{
		MYSQL_ROW row = mysql_fetch_row(result);
		int input_count = std::stoi(row[1]);
		int output_count = std::stoi(row[2]);
		Node* node = get_node(id, ng, input_count, output_count, id2param);
		ng.m_nodes.insert(std::make_pair(node->m_id, node));
		mysql_free_result(result);
	}
	return true;
}