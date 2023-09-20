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


//是否需要加tag来让码农确定是什么类型的曲线, 或者加一个geometry基类, 在此类里面作处理?
template<typename curve_type>
class curve
{
public:
    // using derived_type = typename geo_traits<curve_type>::type;
    using point_number_type = typename geo_traits<curve_type>::point_number_type;
    using point_type = typename geo_traits<curve_type>::point_type;
private:

    Interval<point_number_type> m_interval;

public:

    // curve() = delete;

    ~curve() { };


    ENUM_NURBS point_on_curve(point_number_type u, point_type &point)
    {
        static_cast<curve_type*>(this)->point_on_curve(u, point);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    curve_type* get_derived()
    {
        return static_cast<curve_type*>(this);
    }

    const curve_type *get_derived() const
    {
        return static_cast<const curve_type*>(this);
    }

};
