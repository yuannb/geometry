
#include "node.hpp"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <functional>
#include <random>


using json = nlohmann::ordered_json;

NodeIterator::NodeIterator(const std::vector<InterfaceIter>& node_iter)
{
	m_node_iter = node_iter;
}

bool NodeIterator::is_item_end() const
{
	for (auto& interfaceIter : m_node_iter)
	{
		if (interfaceIter.m_iter.m_kind == KindData::AsItem)
		{
			if (interfaceIter.is_item_end() == false)
			{
				return false;
			}
		}
	}
	return true;
}

NodeIterator& NodeIterator::operator++()
{
	bool need_plus{ true };
	for (auto& interfaceIter : m_node_iter)
	{
		if (interfaceIter.m_iter.m_kind == KindData::AsItem)
		{
			if (interfaceIter.is_item_end() == false)
			{
				need_plus = false;
				++interfaceIter;
			}
		}
	}

	if (need_plus == true)
	{
		for (auto& interfaceIter : m_node_iter)
		{
			if (interfaceIter.is_list_end() == false)
			{
				++interfaceIter;
			}
		}
	}
	return *this;
}

bool NodeIterator::is_end() const
{
	for (auto& interfaceIter : m_node_iter)
	{
		bool f1 = interfaceIter.is_list_end();
		if (f1 == false)
		{
			return false;
		}
		if (interfaceIter.m_iter.m_kind == KindData::AsItem)
		{
			bool f2 = interfaceIter.is_item_end();
			if (f2 == false)
			{
				return false;
			}
		}
	}
	return true;
}

NodeIterator Node::begin()
{
    size_t input_count = m_inputs.size();
    std::vector<InterfaceIter> iters(input_count);
    for (size_t index = 0; index != input_count; ++index)
    {
        iters[index] = m_inputs[index].begin();
    }
    return NodeIterator(iters);
}


bool Node::set_parent_graph(NodeGraph *ng)
{
    m_parent = ng;
    return true;
}

std::string create_uuid()
{
    std::stringstream stream;
    auto random_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 seed_engine(random_seed);
    std::uniform_int_distribution<std::size_t> random_gen ;
    std::size_t value = random_gen(seed_engine);
    stream << std::hex << value;
    
    return stream.str();
}
 bool Node::add_output(const std::vector<Path>& node_paths, const std::vector<std::vector<AbstractNodeData*>>& datas)
{
	 size_t output_count = m_outputs.size();
	 assert(node_paths.size() == output_count);
     for (size_t index = 0; index != node_paths.size(); ++index)
     {
         auto& output_param = m_outputs[index];
         auto it = output_param.m_data.find(node_paths[index]);
         if (it == output_param.m_data.end())
         {
             output_param[node_paths[index]] = datas[index];
         }
         else
         {
             auto& elements = output_param[node_paths[index]];
             elements.insert(elements.end(), datas[index].begin(), datas[index].end());
         }
     }
     return true;
}


bool Node::collect_data()
{
    size_t inputCount = m_inputs.size();
    for (size_t index = 0; index < inputCount; ++index)
    {
        NodeInterface& interface = m_inputs[index];
        interface.collect_data();
    }
    return true;
}

bool Node::calculate_path(const NodeIterator& node_iter, std::vector<Path>& paths)
{
    size_t count = m_outputs.size();
    paths.resize(count);
    bool f = node_iter.is_item_end();
    for (size_t index = 0; index != count; ++index)
    {
        Path& path = paths[index];
        auto& parameter = m_outputs[index];
        if (parameter.m_kind == KindData::AsItem)
        {
			// if (parameter.m_data.size() == 0)
			// {
			// 	path.m_path = std::vector<size_t>{ 0 };
            //     continue;
			// 	// return true;
			// }
			// path.m_path = parameter.m_data.rbegin()->first.m_path;
            if (need_add_output == true)
            {
                size_t param_count = parameter.m_data.size();
                path.m_path = std::vector<size_t>{ param_count };
            }
            else
            {
                if (parameter.m_data.size() == 0)
                {
					path.m_path = std::vector<size_t>{ 0 };
                    continue;
                }
                path.m_path = parameter.m_data.rbegin()->first.m_path;
            }
        }
        else
        {
			size_t param_count = parameter.m_data.size();
			path.m_path = std::vector<size_t>{ param_count };
        }
    }
    need_add_output = false;
    if (f == true)
    {
        need_add_output = true;
    }
    return true;
}

bool Node::calculate()
{
    if (is_calculate == true)
    {
        return true;
    }
    auto it = begin();
    

    while (it.is_end() == false)
    {
        std::vector<Path> paths;
        calculate_path(it, paths);
        std::vector<std::vector<AbstractNodeData*>> node_data = calculate_inner(it);
        add_output(paths, node_data);
        ++it;
    }
	std::vector<Path> paths;
	calculate_path(it, paths);
	std::vector<std::vector<AbstractNodeData*>> node_data = calculate_inner(it);
	add_output(paths, node_data);
    for (auto& input : m_inputs)
    {
        input.m_data.clear();
    }
    return false;
}


json Node::to_json() const
{
	json node_json;
	node_json["inputs"] = json::array();
	json& j = node_json["inputs"];
	for (auto& input : m_inputs)
	{
		j.push_back(input.to_json());
	}
	node_json["outputs"] = json::array();
	json& j_out = node_json["outputs"];
	for (auto& output : m_outputs)
	{
		j_out.push_back(output.to_json());
	}
	return node_json;
}
