
#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "declare.h"
#include "nurbs_tool.h"
class old_surface
{
private:
public:
    old_surface();
    ~old_surface();
    virtual Eigen::Vector3d get_normal() = 0;
};
namespace tnurbs
{
    using namespace tnurbs;
    // template<typename surface_type>
    template<typename T, int dim>
    class surface
    {
    // public:
        // using derived_type = typename geo_traits<surface_type>::type;
        // using point_number_type = typename geo_traits<surface_type>::point_number_type;
        // using point_type = typename geo_traits<surface_type>::point_type;
    public:
        // Interval<point_number_type> m_u_interval;
        // Interval<point_number_type> m_v_interval;
        Box<T, 2> m_interval;

    public:

        // curve() = delete;

        ~surface() { };

        Box<T, 2> get_interval() const { return m_interval; }

        //此函数应写在各个子类中, 做一些合法性检测, 目前先暂时凑合一下
        bool set_interval(T u_min, T u_max, T v_min, T v_max) const 
        {  
            m_interval.Min = Eigen::Vector2<T> { u_min, v_min };
            m_interval.Min = Eigen::Vector2<T> { u_max, v_max };
            return true;
        }

        bool set_interval(const Box<T, 2> &bx) const 
        {  
            m_interval = bx;
            return true;
        }


        // ENUM_NURBS point_on_surface(point_number_type u, point_number_type v,  point_type &point)
        // {
        //     return static_cast<surface_type*>(this)->point_on_surface(u, v, point);
        // }

        // surface_type* get_derived()
        // {
        //     return static_cast<surface_type*>(this);
        // }

        // const surface_type *get_derived() const
        // {
        //     return static_cast<const surface_type*>(this);
        // }
    };

}

