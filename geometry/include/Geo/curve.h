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

        virtual ENGEOMETRYTYPE get_type() const { return ENGEOMETRYTYPE::UNKOWN; }

        virtual ENUM_NURBS point_on_curve(T u, Eigen::Vector<T, dim> &point) const = 0;
        
        virtual ENUM_NURBS get_ends_point(std::array<Eigen::Vector<T, dim>, 2>& points) const = 0;

        virtual Box<T, 1> get_domain_box() const = 0;

        virtual ENUM_NURBS derivative_on_curve(T u, int n, Eigen::Vector<Eigen::Vector<T, dim>, Eigen::Dynamic> &result) const = 0;
        
        virtual bool is_closed() const
        {
            return false;
        }
    };

    using curve3d = curve<double, 3>;

}
