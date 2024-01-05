#include "ReadJson.hpp"
#include "InterpolateCurve.hpp"
#include "MergeNode.hpp"
#include "numberNode.hpp"
#include "pointNode.hpp"

bool CreateNodeGraphFromJson(const json &js, NodeGraph &ng)
{
    json ng_data = js["graph"];
    json ng_nodes = ng_data["nodes"];
    json ng_orders = js["Order"];

    ng.m_id = ng_data["id"];
    std::map<std::string, NodeInterface*> id2param;

    for (auto it = ng_nodes.begin(); it != ng_nodes.end(); ++it)
    {
        std::string type = (*it)["type"].template get<std::string>();
        std::string id = (*it)["id"].template get<std::string>();
        
        if (type == "NumberNode")
        {
            json number_data = (*it)["outputs"];
            std::string number_id = number_data["output"]["id"].template get<std::string>();
            double value = number_data["output"]["value"];
            NumberNode *node = new NumberNode(id, value);
            node->m_outputs[0]->m_id = number_id;

            ng.m_nodes[id] = node;
            id2param[number_id] = node->m_outputs[0];

        }
        else if (type == "MergeNode")
        {
            json merge_data = (*it)["inputs"];
            size_t input_count = merge_data.size();
            MergeNode *node = new MergeNode(id, input_count);
            for (size_t index = 0; index < input_count; ++index)
            {
                std::string key = std::string("Data") + std::to_string(index + 1);
                std::string input_id = merge_data[key]["id"].template get<std::string>();
                node->m_inputs[index]->m_id = input_id;

                id2param[input_id] = node->m_inputs[index];
            }

            json merge_out_data = (*it)["outputs"]["output"];
            std::string output_id = merge_out_data["id"].template get<std::string>();
            node->m_outputs[0]->m_id = output_id;

            ng.m_nodes[id] = node;
            id2param[output_id] = node->m_outputs[0];
        }
        else if (type == "ConstructPoint")
        {
            json point_input_data = (*it)["inputs"];
            ConstructPoint *node = new ConstructPoint(id);
            for (size_t index = 0; index < 3; ++index)
            {
                std::string key = std::string("number") + std::to_string(index + 1);
                std::string input_id = point_input_data[key]["id"].template get<std::string>();
                node->m_inputs[index]->m_id = input_id;

                id2param[input_id] = node->m_inputs[index];
            }

            json point_out_data = (*it)["outputs"]["output"];
            std::string output_id = point_out_data["id"].template get<std::string>();
            node->m_outputs[0]->m_id = output_id;

            ng.m_nodes[id] = node;

            id2param[output_id] = node->m_outputs[0];
        }
        else if (type == "InterpolateCurve")
        {
            json interpolate_input_data = (*it)["inputs"]["points"];
            InterpolateCurve *node = new InterpolateCurve(id);
            std::string input_id = interpolate_input_data["id"].template get<std::string>();
            node->m_inputs[0]->m_id = input_id;

            id2param[input_id] = node->m_inputs[0];

            json interpolate_output_data = (*it)["outputs"]["output"];
            std::string output_id = interpolate_output_data["id"].template get<std::string>();
            node->m_outputs[0]->m_id = output_id;

            ng.m_nodes[id] = node;

            id2param[output_id] = node->m_outputs[0];
            node->m_parent = &ng;
        }
        else
        {
            std::cout << "CreateNodeGraphFromJson : do not have this type node: " << type << std::endl;
        }
    }

    std::map<std::string, std::string> interface_relation;
    
    json connects = ng_data["connections"];
    for (auto it = connects.begin(); it != connects.end(); ++it)
    {
        std::string to_id = (*it)["to"].template get<std::string>();
        std::string from_id = (*it)["from"].template get<std::string>();
        interface_relation[to_id] = from_id;
    }

    for (auto it = ng.m_nodes.begin(); it != ng.m_nodes.end(); ++it)
    {
        size_t inputs_count = it->second->m_inputs.size();
        for (size_t index = 0; index < inputs_count; ++index)
        {
            it->second->m_inputs[index]->m_from_id = interface_relation[it->second->m_inputs[index]->m_id];
            it->second->m_inputs[index]->m_from_interface = id2param[interface_relation[it->second->m_inputs[index]->m_id]];
        }
    }

    json orders = js["Order"];
    for (auto it = orders.begin(); it != orders.end(); ++it)
    {
        std::string id = it->template get<std::string>();
        ng.m_nodes_sort.push_back(id);
    }


    return true;
}
