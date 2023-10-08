#pragma once
#include "nurbs_tool.h"
#include "surface.h"
namespace tnurbs
{
    using namespace tnurbs;
    template<typename T, int dim, int u_degree, int v_degree, bool is_rational>
    class polynomial_surface : public surface<polynomial_surface<T, dim, u_degree, v_degree, is_rational>>
    {

    };

    template<typename T, int dim , bool is_rational>
    class polynomial_surface<T, dim, -1, -1, is_rational> : public surface<polynomial_surface<T, dim, -1, -1, is_rational>>
    {
        static constexpr int point_size = is_rational ? dim + 1 : dim;
        int m_u_degree;
        int m_v_degree;
        Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> m_points;

    public:
        polynomial_surface(const Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> &points)
        {
            m_v_degree = points.rows() - 1;
            m_u_degree = points[0].cols() - 1;
            m_points = points;
        }

        ~polynomial_surface() { }


        ENUM_NURBS point_on_surface(T u, T v, Eigen::Vector<T, dim> &point)
        {
            Eigen::VectorX<T> u_coeff(m_u_degree + 1);
            Eigen::VectorX<T> v_coeff(m_v_degree + 1);

            int min_degree = std::min(m_u_degree, m_v_degree);

            for (int index = 0; index <= min_degree; ++index)
            {
                u_coeff[index] = std::pow(u, index);
                v_coeff[index] = std::pow(v, index);
            }

            for (int index = min_degree + 1; index <= m_u_degree; ++index)
                u_coeff[index] = std::pow(u, index);
            for (int index = min_degree + 1; index <= m_v_degree; ++index)
                v_coeff[index] = std::pow(v, index);

            Eigen::Matrix<T, point_size, Eigen::Dynamic> temp_mat;
            temp_mat.resize(point_size, m_v_degree + 1);

            for (int index = 0; index <= m_v_degree; ++index)
            {
                temp_mat.col(index) = m_points[index] * u_coeff;
            }
            Eigen::Vector<T, point_size> vec = temp_mat * v_coeff;
            point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(vec);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        int get_u_degree() const
        {
            return m_u_degree;
        }
        int get_v_degree() const
        {
            return m_v_degree;
        }

        ENUM_NURBS get_control_points_row(int row, Eigen::Matrix<T, point_size, Eigen::Dynamic> &points) const
        {
            points = m_points[row];
            return ENUM_NURBS::NURBS_SUCCESS;
        }

    };


    template<typename T, int dim, int u_degree, int v_degree, bool is_rational>
    struct geo_traits<polynomial_surface<T, dim, u_degree, v_degree, is_rational> >
    {
        // static constexpr int point_size = is_rational ? dim + 1 : dim;
        // using type = nurbs_surface<T, dim, rows, cols, u_degree, v_degree, is_rational>;
        using point_number_type = T;
        using point_type = typename  Eigen::Vector<T, dim> ;
    };

}
