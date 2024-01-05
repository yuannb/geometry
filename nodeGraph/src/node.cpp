
#include "node.hpp"



bool Node::set_parent_graph(NodeGraph *ng)
{
    m_parent = ng;
    return true;
}

bool Node::add_input(NodeInterface *input)
{
    m_inputs.push_back(input);
    return true;
}

bool Node::add_output(NodeInterface *output)
{
    m_outputs.push_back(output);
    return true;
}

bool Node::collect_data()
{
    size_t inputCount = m_inputs.size();
    for (size_t index = 0; index < inputCount; ++index)
    {

        NodeInterface *interface = m_inputs[index];
        interface->collect_data();
    }
    return true;
}

bool Node::calculate_path(std::vector<DataRelation> &relations) const
{
    size_t inputs_count = m_inputs.size();
    
    std::vector<std::vector<Path>> paths;
    paths.reserve(inputs_count);
    size_t max_path_count = 0; 
    std::vector<size_t> paths_count;
    paths_count.reserve(inputs_count);
    for (size_t index = 0; index < inputs_count; ++index)
    {
        paths.push_back(m_inputs[index]->get_all_path());
        size_t path_count = m_inputs[index]->get_path_count();
        if (path_count > max_path_count)
            max_path_count = path_count;
        paths_count.push_back(path_count);
    }
    relations.resize(max_path_count);
    for (size_t index = 0; index < max_path_count; ++index)
    {
        DataRelation &rel = relations[index];
        rel.m_paths.reserve(inputs_count);
        rel.m_indexs.reserve(inputs_count);
        for (size_t param_index = 0; param_index < inputs_count; ++param_index)
        {
            size_t current_index = index < paths_count[param_index] ? index : paths_count[param_index] - 1;
            rel.m_indexs.push_back(current_index);
            rel.m_paths.push_back(paths[param_index][current_index]);
        }
    }

    //第几个参数是否是asitem
    std::vector<unsigned> item_index;

    for (size_t index = 0; index < inputs_count; ++index)
    {
        bool is_item = m_inputs[index]->get_data_kind() == KindData::AsItem;
        if (is_item == true)
        {
            item_index.push_back(index);
        }
    }
    size_t item_count = item_index.size();
    if (item_count == 0)
        return true;
    for (size_t index = 0; index < max_path_count; ++index)
    {
        DataRelation &rel = relations[index];
        size_t max_values_count = 0;
        std::vector<size_t> values_count;
        values_count.resize(item_count);
        for (size_t i = 0; i < item_count; ++i)
        {
            size_t value_count = m_inputs[item_index[i]]->get_value_count(rel.m_paths[item_index[i]]);
            if (value_count > max_values_count)
                max_values_count = value_count;
            values_count[i] = value_count;
        }
        rel.m_list_indexs.resize(max_values_count);
        for (size_t j = 0; j < max_values_count; ++j)
        {
            rel.m_list_indexs[j].resize(inputs_count, NAN);
            for (size_t k = 0; k < item_count; ++k)
            {
                size_t current_index = j < values_count[k] ? j : values_count[k] - 1;
                rel.m_list_indexs[j][item_index[k]] = current_index;
            }
        }
    }

    return true;
}


