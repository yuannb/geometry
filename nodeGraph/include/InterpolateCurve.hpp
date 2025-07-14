#pragma once
#include "node.hpp"
#include "fit_nurbs.h"

class InterpolateCurve :  public Node
{
public:
    InterpolateCurve() = default;
    InterpolateCurve(const std::string &id);
    virtual ~InterpolateCurve() { }

    virtual std::string get_type() const override
    {
        return NodeType::interpolateCurve;
    }
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);

};

class CurveIntCurve :  public Node
{
public:
    CurveIntCurve() = default;
    CurveIntCurve(const std::string &id);
    virtual ~CurveIntCurve() { }

    virtual std::string get_type() const override
    {
        return NodeType::curvesInt;
    }
    virtual std::vector<std::vector<AbstractNodeData*>> calculate_inner(const NodeIterator& node_iter);

};
