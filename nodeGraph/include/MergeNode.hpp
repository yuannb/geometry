#pragma once
#include "node.hpp"


class MergeNode : public Node
{
public:
    MergeNode() = default;

    MergeNode(const std::string &id, size_t count = 1);
 
    virtual ~MergeNode() 
    { 
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
		for (auto& output : m_outputs)
		{
			output.m_data.clear();
		}
    }

    
    virtual std::string get_type() const override
    {
        return NodeType::merge;
    }

    virtual bool calculate() override;
    
};

class ListItem : public Node
{
public:
    ListItem() = default;
    ListItem(const std::string &id, size_t count = 1);
 
    virtual ~ListItem() 
    { 
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
		for (auto& output : m_outputs)
		{
			output.m_data.clear();
		}
    }

    virtual std::string get_type() const override
    {
        return NodeType::listItem;
    }
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
    
};
class Ellipse : public Node
{
public:
    Ellipse() = default;
    Ellipse(const std::string &id);
 
    virtual ~Ellipse() 
    { 
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
    }

    virtual std::string get_type() const override
    {
        return NodeType::ellipse;
    }

    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
    
};

class Rotate : public Node
{
public:
    Rotate() = default;
    Rotate(const std::string &id);
 
    virtual ~Rotate() 
    { 
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
    }

    virtual std::string get_type() const override
    {
        return NodeType::rotate;
    }

    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
    
};

class Scale : public Node
{
public:
    Scale() = default;
    Scale(const std::string &id);
 
    virtual ~Scale() 
    { 
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
    }

    virtual std::string get_type() const override
    {
        return NodeType::scale;
    }

    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
    
};

class CurveSplit : public Node
{
public:
    CurveSplit() = default;
    CurveSplit(const std::string &id);
 
    virtual ~CurveSplit() 
    { 
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
    }

    virtual std::string get_type() const override
    {
        return NodeType::curveSplit;
    }

    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
    
};

class Loft : public Node
{
public:
    Loft() = default;
    Loft(const std::string &id);
 
    virtual ~Loft() 
    { 
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
    }

    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
    

    virtual std::string get_type() const override
    {
        return NodeType::loft;
    }
    // virtual bool calculate() override;
    
};

class Sweep : public Node
{
public:
    Sweep() = default;
    Sweep(const std::string &id);
 
    virtual ~Sweep() 
    { 
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
    }

    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
    

    virtual std::string get_type() const override
    {
        return NodeType::sweep;
    }
    // virtual bool calculate() override;
    
};
class Flatten : public Node
{
public:
    Flatten() = default;
    Flatten(const std::string &id);
 
    virtual ~Flatten() 
    { 
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
		
        for (auto& output : m_outputs)
		{
			output.m_data.clear();
		}
    }

    virtual std::string get_type() const override
    {
        return NodeType::flatten;
    }

    virtual bool calculate() override;
    
};

class Line : public Node
{
public:
    Line() = default;
    Line(const std::string& id);

    virtual ~Line()
    {
        for (auto& input : m_inputs)
        {
            input.m_data.clear();
        }

    }

    virtual std::string get_type() const override
    {
        return NodeType::line;
    }

    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);

};
