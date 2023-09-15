#pragma once
#include "nurbs_tool.h"
#include <Eigen/Dense>


template<typename T> struct geo_traits;

template<typename T> struct geo_traits<const T> : geo_traits<T> {};

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


template<typename curve_type>
class curve
{
public:
    using type = typename geo_traits<curve_type>::type;
    using point_number_type = typename geo_traits<curve_type>::point_number_type;
    using point_type = typename geo_traits<curve_type>::point_type;
private:
    Interval<point_number_type> m_interval;

public:

    // curve() = delete;

    ~curve() { };


    ENUM_NURBS point_on_curve(point_number_type u, point_type &point)
    {
        static_cast<type*>(this)->point_on_curve(u, point);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
};


