#pragma once
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

enum KindData
{
    AsItem = 0,
    AsList = 1
};


struct Path
{
    std::vector<std::vector<unsigned>> m_path;
    Path() = default;
    Path(const std::vector<std::vector<unsigned>> &path);
    // {
    //     size_t count = path.size();
    //     m_path.resize(count);
    //     for (size_t index = 0; index < count; ++index)
    //     {
    //         std::vector<unsigned> &current = m_path[index];
    //         current.insert(current.end(), path[index].begin(), path[index].end());
    //     }
    // }
    Path(const Path &path);
    // {
    //     size_t count = path.m_path.size();
    //     m_path.resize(count);
    //     for (size_t index = 0; index < count; ++index)
    //     {
    //         std::vector<unsigned> &current = m_path[index];
    //         current.insert(current.end(), path.m_path[index].begin(), path.m_path[index].end());
    //     }
    // }
    bool operator<(const Path &right) const;
    // {
    //     size_t l1 = m_path.size();
    //     size_t l2 = right.m_path.size();
    //     if (l1 > l2)
    //         return false;
    //     if (l1 < l2)
    //         return true;
    //     for (size_t index = 0; index < l1; ++index)
    //     {
    //         size_t k1 = m_path[index].size();
    //         size_t k2 = m_path[index].size();
    //         if (k1 > k2)
    //             return false;
    //         if (k1 < k2)
    //             return true;
    //         for (size_t idx = 0; idx < k1; ++idx)
    //         {
    //             if (m_path[index][idx] > m_path[index][idx])
    //                 return false;
    //             else if (m_path[index][idx] < m_path[index][idx])
    //                 return true;
    //         }
    //     }
    //     return false;
    // }
};


struct DataRelation
{
    std::vector<Path> m_paths;
    std::vector<unsigned> m_indexs;
    std::vector<std::vector<unsigned>> m_list_indexs;
};


struct AbstractNodeData
{
    AbstractNodeData() = default;
    AbstractNodeData(const std::string &str): m_type(str) { }
    std::string m_type;
};


struct NumberData : public AbstractNodeData
{
    double m_data;
    NumberData(double num): m_data(num) { m_type = std::string("NumberData"); };
    ~NumberData() { }
};

struct IntegerData : public AbstractNodeData
{
    int m_data;
    IntegerData(int num): m_data(num) { m_type = std::string("IntegerData"); }
    ~IntegerData() { }
};

struct PointData : public AbstractNodeData
{
    Eigen::Vector3d m_data;
    PointData(double x, double y, double z): m_data(x, y, z) { m_type = std::string("PointData"); }
    ~PointData() { }
};

struct GeometryData : public AbstractNodeData
{
    std::string m_data;
    // std::string m_mesh_id;
    GeometryData(const std::string &geo_id): m_data(geo_id) { m_type = std::string("GeometryData"); }
    ~GeometryData() { }
};




class NodeInterface
{
public:

    //interface id
    std::string m_id;

    //from interface id
    std::string m_from_id;

    //from interface
    const NodeInterface *m_from_interface;


    //the below two var may be deleted?
    //to interface id
    std::vector<std::string> m_to_id;
    // to interface
    std::vector<NodeInterface *> m_end_interface;

    KindData m_kind;
    std::map<Path, std::vector<AbstractNodeData*>> m_data;
public:
    NodeInterface(const std::string &id) : m_id(id) { }
    NodeInterface() = default;
    ~NodeInterface( ) { }
    NodeInterface(const std::string &id, const std::string &from_id) : m_id(id), m_from_id(from_id) { }

    std::vector<AbstractNodeData*> &insert_path(const Path &path);
    // {
    //     return m_data[path];
    // }

    bool insert_element(const Path &path, AbstractNodeData* data);
    // {
    //     m_data[path].push_back(data);
    //     return true;
    // }

    bool collect_data();
    // {
    //     //暂时不判断类型
    //     m_data.clear();
    //     m_data.insert(m_data.end(), m_from_interface->m_data.begin(), m_from_interface->m_data.end());
    //     return true;
    // }

    size_t get_path_count() const;// { return m_data.size(); }

    std::vector<Path> get_all_path() const;
    // {
    //     std::vector<const Path &> paths;
    //     auto it = m_data.begin();
    //     auto end = m_data.end();
    //     for (; it != end; ++it)
    //     {
    //         paths.push_back(it->first);
    //     }
    //     return paths;
    // }

    size_t get_value_count(const Path &path);
    // {
    //     return m_data[path].size();
    // }

    KindData get_data_kind() const;// { return m_kind;  }

    bool set_kind_data(const KindData &kind);// { m_kind = kind; return true; }

    bool set_from_id(const std::string &from_id);
    // {
    //     m_from_id = from_id;
    //     return true;
    // }
    bool set_from_interface(const NodeInterface *interface);
    // {
    //     m_from_interface = interface;
    //     return true;
    // }

    bool set_data(const DataRelation &rel, std::vector<AbstractNodeData *> datas);
    // {
        
    //     if (m_kind == KindData::AsList)
    //     {
    //         size_t count = rel.m_list_indexs.size();
    //         if (datas.size() != count)
    //             return false;
    //         size_t path_count = rel.m_paths.size();
    //         Path path;
            
    //         for (size_t index = 0; index < count; ++index)
    //         {
    //             path.m_path.clear();
    //             path.m_path.resize(path_count);
    //             for (size_t i = 0; i < path_count; ++i)
    //             {
    //                 path.m_path[i] = { rel.m_indexs[i] };
    //                 if (isnan(rel.m_list_indexs[index]) == false)
    //                 {
    //                     path.m_path[i].push_back(rel.m_list_indexs[index]);
    //                 }
    //             }
    //             m_data[path].push_back(datas[index]);
    //         }
    //     }
    //     else
    //     {
    //         size_t path_count = rel.m_paths.size();
    //         Path path;
    //         // path.m_path.clear();
    //         path.m_path.resize(path_count);
    //         for (size_t i = 0; i < path_count; ++i)
    //         {
    //             path.m_path[i] = { rel.m_indexs[i] };
    //             if (isnan(rel.m_list_indexs[i]) == false)
    //             {
    //                 path.m_path[i].push_back(rel.m_list_indexs[i]);
    //             }
    //         }
    //         m_data[path] = datas;
    //     }
    //     return true;
    // }

    bool set_data(const DataRelation &rel, std::vector<std::vector<AbstractNodeData *>> datas);
    // {
        
    //     if (m_kind == KindData::AsList)
    //     {
    //         size_t count = rel.m_list_indexs.size();
    //         if (datas.size() != count)
    //             return false;
    //         size_t path_count = rel.m_paths.size();
    //         Path path;
            
    //         for (size_t index = 0; index < count; ++index)
    //         {
    //             path.m_path.clear();
    //             path.m_path.resize(path_count);
    //             for (size_t i = 0; i < path_count; ++i)
    //             {
    //                 path.m_path[i] = { rel.m_indexs[i] };
    //                 if (std::isnan(rel.m_list_indexs[index]) == false)
    //                 {
    //                     path.m_path[i].push_back(rel.m_list_indexs[index][i]);
    //                 }
    //             }
    //             m_data[path] = datas[index];
    //         }
    //     }
    //     else
    //     {
    //         return false;
    //     }
    //     return true;
    // }



    // NodeInterface *get_from_interface();// { return m_from_interface; }

    const NodeInterface *get_from_interface() const;// { return m_from_interface; };

};





