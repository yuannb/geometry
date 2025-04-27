
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
    public:

        ~surface() { };

        virtual ENUM_NURBS derivative_on_surface(int n, T u, T v, Eigen::MatrixX<Eigen::Vector<T, dim>>& result) const = 0;

		virtual	ENUM_NURBS point_on_surface(T u, T v, Eigen::Vector<T, dim>& point) const = 0;
		
        virtual	Box<T, 2> get_domain() const = 0;

		virtual bool is_period(std::array<bool, 2>& periods) const = 0;
    };

}

