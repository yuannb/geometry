#pragma once
#include "node.hpp"

class ConstructPoint : public Node
{
public:
    ConstructPoint() = default;
    ConstructPoint(const std::string &id);
    virtual ~ConstructPoint() { }

    virtual std::string get_type() const override
    {
        return NodeType::constructPoint;
    }
    
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
};

class DistancePt : public Node
{
public:
    DistancePt() = default;
    DistancePt(const std::string &id);
    virtual ~DistancePt() { }

    virtual std::string get_type() const override
    {
        return NodeType::distancePt;
    }
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
};

class EndsPoints : public Node
{
public:
    EndsPoints() = default;
    EndsPoints(const std::string &id);
    virtual ~EndsPoints() 
    {
        for (auto& input : m_inputs)
        {
            input.m_data.clear();
        }
    }

    virtual std::string get_type() const override
    {
        return NodeType::EndPts;
    }
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
};

class PrepFrame : public Node
{
public:
    PrepFrame() = default;
    PrepFrame(const std::string &id);
    virtual ~PrepFrame() 
    {
        for (auto& input : m_inputs)
        {
            input.m_data.clear();
        }
    }

    virtual std::string get_type() const override
    {
        return NodeType::prepFrame;
    }
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
};

class Rectangle : public Node
{
public:
    Rectangle() = default;
    Rectangle(const std::string &id);
    virtual ~Rectangle() 
    {
        for (auto& input : m_inputs)
        {
            input.m_data.clear();
        }
    }

    virtual std::string get_type() const override
    {
        return NodeType::rectangle;
    }
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
};
