#pragma once
#include "nurbs_tool.h"

template<typename T, int dim, bool is_rational, int points_count, int degree>
class polynomial_curve
{

};

template<typename T, int dim, bool is_rational>
class polynomial_curve<T, dim, is_rational, -1, -1>
{
private:
    static constexpr int point_size = is_rational ? dim + 1 : dim;
    int m_degree;
    
    Eigen::Matrix<T, point_size, Eigen::Dynamic> m_points;
    Interval<T> m_interval;

public:
    using point_type = typename Eigen::Vector<T, point_size>;

public: 

    polynomial_curve(const polynomial_curve<T, dim, is_rational, -1, -1> &curve)
    {
        m_degree = curve.m_degree;
        m_points = curve.m_points;
        m_interval = curve.m_interval;
    }

    ENUM_NURBS point_on_curve(T u, Eigen::Vector<T, dim> &point)
    {
        Eigen::VectorX<T> coeff(m_degree + 1);
        for (int i = 0; i <= m_degree; ++i)
        {
            coeff[i] = std::pow(u, i);
        }
        Eigen::Vector<T, point_size> vec = m_points * coeff;
        point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(vec);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
};