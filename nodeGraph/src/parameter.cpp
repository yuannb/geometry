#include "parameter.hpp"

Path::Path(const std::vector<std::vector<unsigned>> &path)
{
    size_t count = path.size();
    m_path.resize(count);
    for (size_t index = 0; index < count; ++index)
    {
        std::vector<unsigned> &current = m_path[index];
        current.insert(current.end(), path[index].begin(), path[index].end());
    }
}
Path::Path(const Path &path)
{
    size_t count = path.m_path.size();
    m_path.resize(count);
    for (size_t index = 0; index < count; ++index)
    {
        std::vector<unsigned> &current = m_path[index];
        current.insert(current.end(), path.m_path[index].begin(), path.m_path[index].end());
    }
}

bool Path::operator<(const Path &right) const
{
    size_t l1 = m_path.size();
    size_t l2 = right.m_path.size();
    if (l1 > l2)
        return false;
    if (l1 < l2)
        return true;
    for (size_t index = 0; index < l1; ++index)
    {
        size_t k1 = m_path[index].size();
        size_t k2 = m_path[index].size();
        if (k1 > k2)
            return false;
        if (k1 < k2)
            return true;
        for (size_t idx = 0; idx < k1; ++idx)
        {
            if (m_path[index][idx] > m_path[index][idx])
                return false;
            else if (m_path[index][idx] < m_path[index][idx])
                return true;
        }
    }
    return false;
}

std::vector<AbstractNodeData*> & NodeInterface::insert_path(const Path &path)
{
    return m_data[path];
}

bool NodeInterface::insert_element(const Path &path, AbstractNodeData* data)
{
    m_data[path].push_back(data);
    return true;
}

bool NodeInterface::collect_data()
{
    //暂时不判断类型
    m_data.clear();
    m_data.insert(m_from_interface->m_data.begin(), m_from_interface->m_data.end());
    return true;
}

size_t NodeInterface::get_path_count() const { return m_data.size(); }

std::vector<Path> NodeInterface::get_all_path() const
{
    std::vector<Path> paths;
    auto it = m_data.begin();
    auto end = m_data.end();
    for (; it != end; ++it)
    {
        paths.push_back(it->first);
    }
    return paths;
}

size_t NodeInterface::get_value_count(const Path &path)
{
    return m_data[path].size();
}

KindData NodeInterface::get_data_kind() const { return m_kind;  }

bool NodeInterface::set_kind_data(const KindData &kind) { m_kind = kind; return true; }

bool NodeInterface::set_from_id(const std::string &from_id)
{
    m_from_id = from_id;
    return true;
}
bool NodeInterface::set_from_interface(const NodeInterface *interface)
{
    m_from_interface = interface;
    return true;
}

bool NodeInterface::set_data(const DataRelation &rel, std::vector<AbstractNodeData *> datas)
{
    
    if (m_kind == KindData::AsList)
    {
        size_t count = rel.m_list_indexs.size();

        size_t path_count = rel.m_indexs.size();
        Path path;


        if (count == 0)
        {
            path.m_path.clear();
            path.m_path.resize(path_count);
            for (size_t index = 0; index < path_count; ++index)
            {
                path.m_path[index] = { rel.m_indexs[index] };
            }
            m_data[path].insert(m_data[path].end(), datas.begin(), datas.end());
        }
        else
        {
            if (datas.size() != count)
               return false;
            for (size_t index = 0; index < count; ++index)
            {
                path.m_path.clear();
                path.m_path.resize(path_count);
                for (size_t i = 0; i < path_count; ++i)
                {
                    path.m_path[i] = { rel.m_indexs[i] };
                    if (std::isnan(rel.m_list_indexs[index][i]) == false)
                    {
                        path.m_path[i].push_back(rel.m_list_indexs[index][i]);
                    }
                }
                m_data[path].push_back(datas[index]);
            }
        }

        // std::vector<Path> paths(count);

        // for (size_t index = 0; index < count; ++index)
        // {
        //     path.m_path.clear();
        //     path.m_path.resize(path_count);
        //     for (size_t i = 0; i < path_count; ++i)
        //     {
        //         path.m_path[i] = { rel.m_indexs[i] };
        //         if (std::isnan(rel.m_list_indexs[index][i]) == false)
        //         {
        //             path.m_path[i].push_back(rel.m_list_indexs[index][i]);
        //         }
        //     }
        //     m_data[path].push_back(datas[index]);
        // }
    }
    else
    {
        size_t path_count = rel.m_indexs.size();
        Path path;
        // path.m_path.clear();
        path.m_path.resize(path_count);
        for (size_t i = 0; i < path_count; ++i)
        {
            path.m_path[i] = { rel.m_indexs[i] };
            // if (std::isnan(rel.m_list_indexs[i][0]) == false)
            // {
            //     path.m_path[i].push_back(rel.m_list_indexs[i][0]);
            // }
        }
        m_data[path] = datas;
    }
    return true;
}

bool NodeInterface::set_data(const DataRelation &rel, std::vector<std::vector<AbstractNodeData *>> datas)
{
    
    if (m_kind == KindData::AsList)
    {
        size_t count = rel.m_list_indexs.size();
        if (datas.size() != count)
            return false;
        size_t path_count = rel.m_paths.size();
        Path path;
        
        for (size_t index = 0; index < count; ++index)
        {
            path.m_path.clear();
            path.m_path.resize(path_count);
            for (size_t i = 0; i < path_count; ++i)
            {
                path.m_path[i] = { rel.m_indexs[i] };
                if (std::isnan(rel.m_list_indexs[index][0]) == false)
                {
                    path.m_path[i].push_back(rel.m_list_indexs[index][i]);
                }
            }
            m_data[path] = datas[index];
        }
    }
    else
    {
        return false;
    }
    return true;
}


// const NodeInterface *NodeInterface::get_from_interface() { return m_from_interface; }

const NodeInterface *NodeInterface::get_from_interface() const { return m_from_interface; }


