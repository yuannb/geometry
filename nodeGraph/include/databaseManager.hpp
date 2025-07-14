#pragma once
#include "NodeGraph.hpp"
#include "logger_manager.hpp"
#include "serialization.hpp"
#include "nodeFactory.hpp"

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#include <iostream>

#include "mysql.h" 

//TODO: delete
class DatabaseException : public std::runtime_error {
public:
    DatabaseException(const std::string& message) : std::runtime_error(message) {}
};

class NodeGraphDBManager 
{
private:
    MYSQL* conn;

    const char* host = "127.0.0.1";
    const char* user = "root";
    const char* password = "";
    const char* dbname = "database_test";
    const int port = 3306;

    std::shared_ptr<spdlog::logger> logger = LoggerManager::getLogger();

    // 私有工具方法
    void check_error(int status);

    void escape_string(std::string& str);

public:
    NodeGraphDBManager();

    ~NodeGraphDBManager();

    // 创建必要的表
    void create_tables();

    void add_geometry(const std::string& node_uuid, const NodeGraph& ng);

    // 添加 NodeGraph
    void add_nodegraph(const NodeGraph& nodegraph);

    // 添加 Node
    void add_node(const Node& node, const std::string& nodegraph_uuid);

    // 添加 Parameter
    void add_parameter(const NodeInterface& param, bool is_input, int index, const Node& node);

    // 查询 geometry
    bool get_geometry(const std::string& gid, NodeGraph& ng);

    bool is_ng_exits(const std::string& ng_id);

    bool update_ng_time(const std::string& ng_id);

    std::string get_ng_date(const std::string& ng_id);

    // 查询 NodeGraph
    NodeGraph get_nodegraph(const std::string& uuid);

    // 查询 Node
    Node* get_node(const std::string& uuid, NodeGraph& ng, int input_count, int output_count, std::map<std::string, NodeInterface*>& id2param);

    // 查询 Parameter
    bool get_parameter(const std::string& uuid, NodeInterface& param, Node* node);

    // 查询 geometry
    bool delete_geometry(const std::string& nid);

    bool find_nodes_id(const std::string& ng_id, std::vector<std::string>& nodes_id);

    // 查询 NodeGraph
    bool delete_nodegraph(const std::string& uuid);

    bool find_paramter_id(const std::string& node_id, std::vector<std::string>& paramters_id);

    // 查询 Node
    bool delete_node(const std::string& node_id);

    // 查询 Parameter
    bool delete_parameter(const std::string& uuid);

    bool CreateNodeFromDB(const nlohmann::ordered_json node_js, NodeGraph& ng, std::map<std::string, NodeInterface*>& id2param);
};
