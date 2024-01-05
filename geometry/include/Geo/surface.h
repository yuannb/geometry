
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
    template<typename surface_type>
    class surface
    {
    public:
        // using derived_type = typename geo_traits<surface_type>::type;
        using point_number_type = typename geo_traits<surface_type>::point_number_type;
        using point_type = typename geo_traits<surface_type>::point_type;
    private:
        Interval<point_number_type> m_u_interval;
        Interval<point_number_type> m_v_interval;

    public:

        // curve() = delete;

        ~surface() { };


        ENUM_NURBS point_on_surface(point_number_type u, point_number_type v,  point_type &point)
        {
            return static_cast<surface_type*>(this)->point_on_surface(u, v, point);
        }

        surface_type* get_derived()
        {
            return static_cast<surface_type*>(this);
        }

        const surface_type *get_derived() const
        {
            return static_cast<const surface_type*>(this);
        }
    };

}

