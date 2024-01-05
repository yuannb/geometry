#pragma once
#include "nurbs_tool.h"
#include <Eigen/Dense>
#include "declare.h"

class old_curve
{
private:
    /* data */
    // int type;
public:
    old_curve(/* args */);
    ~old_curve();
    virtual Eigen::Vector3d get_start_point() = 0;
    virtual Eigen::Vector3d get_end_point() = 0;
};

namespace tnurbs
{
    template<typename T, int dim>
    class curve
    {
    public:
        using Type = T;
        static constexpr int dimension = dim;
    protected:

        //辅助计算使用, 一般不用
        Box<T, 1> m_interval;

    public:
        virtual  ~curve() { };

        constexpr virtual ENGEOMETRYTYPE get_type() const { return ENGEOMETRYTYPE::UNKOWN; }

        virtual ENUM_NURBS point_on_curve(T u, Eigen::Vector<T, dim> &point) const = 0;

    };

    using curve3d = curve<double, 3>;

}
