#pragma once
#include "node.hpp"


class NumberNode : public Node
{
public:
    NumberNode() = default ;
    NumberNode(const std::string &id);
    virtual ~NumberNode()
    {
		for (auto& input : m_inputs)
		{
            if (input.m_from_interface != nullptr)
            {
				input.m_data.clear();
            }
		}
    }

    virtual std::string get_type() const override
    {
        return NodeType::numberNode;
    }

    virtual bool collect_data();
    
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter)
    {
        return { {nullptr} };
    }
    virtual bool calculate() override;
    
};

class Ranges : public Node
{
public:
    Ranges() = default ;
    Ranges(const std::string &id);
    virtual ~Ranges() 
    {
		for (auto& input : m_inputs)
		{
            if (input.m_from_interface != nullptr)
            {
				input.m_data.clear();
            }
		}
    }

    virtual std::string get_type() const override
    {
        return NodeType::ranges;
    }
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
    
};

class Extrude : public Node
{
public:
    Extrude() = default ;
    Extrude(const std::string &id);
    virtual ~Extrude()
    {
		for (auto& input : m_inputs)
		{
			input.m_data.clear();
		}
    }

    virtual std::string get_type() const override
    {
        return NodeType::extrude;
    }
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);
    
};
