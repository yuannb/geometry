#pragma once
#include "node.hpp"


class NumberNode : public Node
{
public:
    NumberNode() = default ;
    NumberNode(const std::string &id, double num);
    virtual ~NumberNode() { }


    virtual bool collect_data();
    // {
        // return true;
    // }

    virtual bool calculate() override;
    // {
    //     return true;
    // }
    
};

