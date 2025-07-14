#include "parameter.hpp"

using json = nlohmann::ordered_json;
bool Path::operator<(const Path &right) const
{
    size_t l1 = m_path.size();
    size_t l2 = right.m_path.size();
    if (l1 > l2)
        return false;
    if (l1 < l2)
        return true;
	for (size_t idx = 0; idx < l1; ++idx)
	{
		if (m_path[idx] < right.m_path[idx])
			return true;
		else if (m_path[idx] < m_path[idx])
			return false;
	}
    return false;
}

InterfaceIter NodeInterface::begin()
{
    InterfaceIter iter = InterfaceIter();
    iter.m_iter.m_index = 0;
    iter.m_iter.m_kind = m_kind;
    iter.m_iter.m_path_iter = iter.m_iter.m_path_begin = m_data.begin();
    iter.m_iter.m_path_end = m_data.end();
    return iter;
}


bool NodeInterface::collect_data()
{
    if (m_from_interface != nullptr)
    {
        for (auto it = m_data.begin(); it != m_data.end(); ++it)
        {
            size_t count = it->second.size();
            for (size_t index = 0; index < count; ++index)
            {
                AbstractNodeData* data = it->second[index];
                if (data != nullptr)
                {
                    delete data;
                }
            }
        }
		m_data.clear();
		m_data.insert(m_from_interface->m_data.begin(), m_from_interface->m_data.end());
    }
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

