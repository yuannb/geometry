#pragma once
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "json.hpp"

enum KindData
{
    AsItem = 0,
    AsList = 1
};


struct Path
{
    std::vector<size_t> m_path;
    Path() = default;
    Path(const std::vector<size_t>& path) : m_path(path) {};
    Path(const Path& path): m_path(path.m_path){};
    bool operator<(const Path &right) const;
    std::string to_string() const
    {
        std::string path = "{ ";
        size_t path_count = m_path.size();
        for (size_t index = 0; index + 1 < path_count; ++index)
        {
            path += std::to_string(m_path[index]);
            path += ", ";
        }

        path += std::to_string(m_path.back());
        path += " }";
        return path;
    }
};


struct NodeDataType
{
    static constexpr const char*  unknown = "unknown";
    static constexpr const char*  numberData = "numberData";
    static constexpr const char*  integerData = "integerData";
    static constexpr const char*  pointData = "pointData";
    static constexpr const char*  geometryData = "geometryData";
    static constexpr const char*  booleanData = "booleanData";
};

struct ParamterType
{
    static constexpr const char*  curve = "curve";
    static constexpr const char*  id = "id";
    static constexpr const char*  dir = "dir";

    static constexpr const char*  value = "value";
    static constexpr const char*  surf = "surf";
    static constexpr const char*  inputs = "inputs";
    static constexpr const char*  outputs = "outputs";
    static constexpr const char*  output = "output";
    static constexpr const char*  xdir = "xdir";
    static constexpr const char*  ydir = "ydir";
    static constexpr const char*  normal = "normal";
    static constexpr const char*  origin = "origin";
    static constexpr const char*  major = "major";
    static constexpr const char*  minor = "minor";

    static constexpr const char*  geometry = "Geometry";
    static constexpr const char*  angle = "Angle";
    
    static constexpr const char*  low = "low";
    static constexpr const char*  high = "high";
    static constexpr const char*  count = "count";
    
    static constexpr const char*  series = "series";

    static constexpr const char*  factor = "factor";
    static constexpr const char*  tolerance = "tolerance";
    static constexpr const char*  loft = "Loft";
    static constexpr const char*  input = "input";
    static constexpr const char*  index = "index";
    static constexpr const char*  wrap = "wrap";
    static constexpr const char*  curve1 = "curve1";
    static constexpr const char*  curve2 = "curve2";
    static constexpr const char*  startPoint = "startPoint";
    static constexpr const char*  endPoint = "endPoint";
    static constexpr const char*  line = "line";
    static constexpr const char*  points = "points";
    static constexpr const char*  pt1 = "pt1";
    static constexpr const char*  pt2 = "pt2";
    static constexpr const char*  dist = "dist";
    static constexpr const char*  startParam = "startParam";
    static constexpr const char*  endParam = "endParam";
    static constexpr const char*  param = "param";
    static constexpr const char*  ulength = "ulength";
    static constexpr const char*  vlength = "vlength";
    static constexpr const char*  rail = "Rail";
    static constexpr const char*  section = "section";
    static constexpr const char*  type = "type";
    static constexpr const char*  data = "data";
    static constexpr const char*  geoemtryType = "geoemtryType";
    static constexpr const char*  curveType = "curveType";
    static constexpr const char*  surfType = "surfType";
    static constexpr const char*  knots = "knots";
    static constexpr const char*  uKnots = "uknots";
    static constexpr const char*  vKnots = "vknots";
    static constexpr const char*  controlPoints = "controlPoints";
    static constexpr const char*  degree = "degree";
    static constexpr const char*  uDegree = "uDegree";
    static constexpr const char*  vDegree = "vDegree";
    static constexpr const char*  point = "point";
    static constexpr const char*  geometryData = "geometryData";
    static constexpr const char*  nodes = "nodes";
    static constexpr const char*  order = "order";
    static constexpr const char*  dataType = "type";
    static constexpr const char*  status = "status";
    static constexpr const char*  original = "original";
    static constexpr const char*  changed_nodes = "changed_nodes";
    static constexpr const char*  deleted_nodes = "deleted_nodes";
    static constexpr const char*  depended_nodes = "depended_nodes";
    


};

struct AbstractNodeData
{
    AbstractNodeData() = default;
    AbstractNodeData(const std::string &str): m_type(str) { }
    std::string m_type;
    virtual ~AbstractNodeData() {};
    virtual nlohmann::ordered_json to_json() const = 0;
};


struct NumberData : public AbstractNodeData
{
    double m_data;
    NumberData(double num): m_data(num) { m_type = NodeDataType::numberData; };
    ~NumberData() { }
    virtual nlohmann::ordered_json to_json() const
    {
        return { {"type", m_type}, {"data", m_data} };
    }
};


struct BooleanData : public AbstractNodeData
{
    double m_data;
    BooleanData(double num): m_data(num) { m_type = std::string("BooleanData"); };
    ~BooleanData() { }
    virtual nlohmann::ordered_json to_json() const
    {
        return { {"type", "BooleanData"}, {"data", m_data} };
    }
};
struct IntegerData : public AbstractNodeData
{
    int m_data;
    IntegerData(int num): m_data(num) { m_type = std::string("IntegerData"); }
    ~IntegerData() { }
    virtual nlohmann::ordered_json to_json() const
    {
        return { {"type", "IntegerData"}, {"data", m_data} };
    }
};

struct PointData : public AbstractNodeData
{
    Eigen::Vector3d m_data;
    PointData(double x, double y, double z): m_data(x, y, z) { m_type = NodeDataType::pointData; }
    PointData(const Eigen::Vector3d& vec): m_data(vec) { m_type = NodeDataType::pointData; }
    ~PointData() { }
    virtual nlohmann::ordered_json to_json() const
    {
        return { {"type", NodeDataType::pointData}, {"data", {m_data[0], m_data[1], m_data[2]}}};
    }
};

struct GeometryData : public AbstractNodeData
{
    std::string m_data;
    int m_geoemtryType = -1;
    // std::string m_mesh_id;
    GeometryData(const std::string &geo_id): m_data(geo_id) { m_type = NodeDataType::geometryData; }
    ~GeometryData() { }
    virtual nlohmann::ordered_json to_json() const
    {
        return { {"type", NodeDataType::geometryData}, {"data", m_data}, {"geoemtryType", m_geoemtryType} };
    }
};


struct IterfaceIterInner
{
    std::map<Path, std::vector<AbstractNodeData*>>::iterator m_path_iter;
    std::map<Path, std::vector<AbstractNodeData*>>::iterator m_path_begin;
    std::map<Path, std::vector<AbstractNodeData*>>::iterator m_path_end;
    size_t m_index;
    KindData m_kind;
};


class InterfaceIter
{
public:
	IterfaceIterInner m_iter;
public:
    InterfaceIter() = default;
	explicit InterfaceIter(IterfaceIterInner iter)
	{
		m_iter.m_path_iter = iter.m_path_iter;
		m_iter.m_path_begin = iter.m_path_begin;
		m_iter.m_path_end = iter.m_path_end;
		m_iter.m_index = iter.m_index;
		m_iter.m_kind = iter.m_kind;
	}

	InterfaceIter& operator++()
	{
		if (m_iter.m_kind == KindData::AsItem)
		{
			size_t count = m_iter.m_path_iter->second.size();
			if (m_iter.m_index + 1 < count)
			{
				m_iter.m_index += 1;
			}
			else
			{
                ++m_iter.m_path_iter;
                m_iter.m_index = 0;
			}
		}
		else if (m_iter.m_kind == KindData::AsList)
		{
			if (m_iter.m_path_iter != m_iter.m_path_end)
			{
				++m_iter.m_path_iter;
			}
		}
		else
		{
			assert(false);
		}
        return *this;
	}
    bool is_list_end() const
    {
        auto current_it = m_iter.m_path_iter;
        ++current_it;
        if (current_it == m_iter.m_path_end)
        {
            return true;
        }
        return false;
    }

    bool is_item_end() const
    {
        if (m_iter.m_path_begin == m_iter.m_path_end)
        {
            return true;
        }
        return m_iter.m_index + 1 >= m_iter.m_path_iter->second.size();
    }

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
    ~NodeInterface( ) 
    {
        // for (auto data : m_data)
        // {
        //     auto vecs = data.second;
        //     for (auto vec : vecs)
        //     {
        //         delete vec;
        //     }
        // }
    }
    NodeInterface(const std::string &id, const std::string &from_id) : m_id(id), m_from_id(from_id) { }

    InterfaceIter begin();
    std::vector<AbstractNodeData*>& operator[](const Path& path) 
    {
        Path temp = path;
        return m_data[path]; 
    }
    std::vector<AbstractNodeData> &insert_path(const Path &path);

    bool insert_element(const Path &path, AbstractNodeData data);

    bool collect_data();

    size_t get_path_count() const;// { return m_data.size(); }

    std::vector<Path> get_all_path() const;

    size_t get_value_count(const Path &path);

    KindData get_data_kind() const;
    bool set_kind_data(const KindData &kind);
    bool set_from_id(const std::string &from_id);

    bool set_from_interface(const NodeInterface *interface);

    nlohmann::ordered_json to_json() const
    {
        nlohmann::ordered_json paramter;
        paramter[m_id] = {};
        auto& paramter_inner = paramter[m_id];
        paramter_inner["from_id"] = m_from_id;
        paramter_inner["to_id"] = {};
        for (auto& to_id : m_to_id)
        {
            paramter_inner["to_id"].push_back(to_id);
        }
        paramter_inner["data"] = {};
        auto& param_data = paramter_inner["data"];
        for (auto it = m_data.begin(); it != m_data.end(); ++it)
        {
            std::string str = it->first.to_string();
            param_data[str] = {};
            auto& data_inner = param_data[str];
            size_t data_count = it->second.size();
            for (size_t index = 0; index != data_count; ++index)
            {
                data_inner.push_back(it->second[index]->to_json());
            }
        }

        return paramter;
    }

};





