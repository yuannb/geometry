#pragma once
#include <string>
#include <vector>
#include <type_traits>
#include "magic_enum.hpp"
#include "parameter.hpp"
#include "json.hpp"

class NodeGraph;



class NodeIterator
{
public:
	std::vector<InterfaceIter> m_node_iter;
    bool m_is_end{ false };
public:
    explicit NodeIterator(const std::vector<InterfaceIter>& node_iter);

    bool is_item_end() const;

    NodeIterator& operator++();

    bool is_end() const;
};


struct NodeType
{
    static constexpr const char* interpolateCurve = "InterpolateCurve";
    static constexpr const char* numberNode = "NumberNode";
    static constexpr const char* merge = "Merge";
    static constexpr const char* constructPoint = "ConstructPoint";
    static constexpr const char* extrude = "Extrude";
    static constexpr const char* ellipse = "Ellipse";
    static constexpr const char* rotate = "Rotate";
    static constexpr const char* ranges = "Ranges";
    static constexpr const char* scale = "Scale";
    static constexpr const char* curveSplit = "CurveSplit";
    static constexpr const char* loft = "Loft";
    static constexpr const char* flatten = "Flatten";
    static constexpr const char* listItem = "ListItem";
    static constexpr const char* line = "Line";
    static constexpr const char* curvesInt = "CurvesInt";
    static constexpr const char* distancePt = "DistancePt";
    static constexpr const char* EndPts = "EndPts";
    static constexpr const char* prepFrame = "PrepFrame";
    static constexpr const char* rectangle = "Rectangle";
    static constexpr const char* sweep = "Sweep";
};

std::string create_uuid();

class Node 
{
public:

    NodeGraph *m_parent;
    
    std::string m_id;
    std::vector<NodeInterface> m_inputs;
    std::vector<NodeInterface> m_outputs;

private:
    bool need_add_output{false};
public:
    bool is_calculate = false;
public:

    Node() = default;
    Node(const std::string &id): m_id(id) { }

    virtual ~Node() 
    {
    }
    
    std::vector<AbstractNodeData*> get_data_list(size_t index, const NodeIterator& node_iter)
    {
		auto iter = node_iter.m_node_iter[index].m_iter;
        return iter.m_path_iter->second;
    }
    
    AbstractNodeData* get_data_item(size_t index, const NodeIterator& node_iter)
    {
		auto iter = node_iter.m_node_iter[index].m_iter;
		return iter.m_path_iter->second[iter.m_index];
    }
    NodeIterator begin();

    bool set_parent_graph(NodeGraph *ng);

    bool add_output(const std::vector<Path>& node_iter, const std::vector<std::vector<AbstractNodeData*>>& datas);

    virtual std::string get_type() const 
    {
        return "unknown";
    }

    bool collect_data();

    size_t get_input_count() { return m_inputs.size(); }


     bool calculate_path(const NodeIterator& node_iter, std::vector<Path>& paths);


    virtual bool calculate();
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter)
    {
        return { {nullptr} };
    }

    nlohmann::ordered_json to_json() const;

};


