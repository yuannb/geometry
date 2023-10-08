#pragma once
#include "nurbs_tool.h"
namespace tnurbs
{
    using namespace tnurbs;
    template<typename T, int dim, bool is_rational, int points_count, int degree>
    class polynomial_curve : public curve<polynomial_curve<T, dim, is_rational, points_count, degree>>
    {

    };

    template<typename T, int dim, bool is_rational>
    class polynomial_curve<T, dim, is_rational, -1, -1> : public curve<polynomial_curve<T, dim, is_rational, -1, -1>>
    {
    private:
        static constexpr int point_size = is_rational ? dim + 1 : dim;
        int m_degree;
        
        Eigen::Matrix<T, point_size, Eigen::Dynamic> m_points;

    public: 

        polynomial_curve(const polynomial_curve<T, dim, is_rational, -1, -1> &curve)
        {
            m_degree = curve.m_degree;
            m_points = curve.m_points;
        }

        ~polynomial_curve() { }

        polynomial_curve(const Eigen::Matrix<T, point_size, Eigen::Dynamic> &points)
        {
            m_degree = points.cols() - 1;
            m_points = points;
        }

        ENUM_NURBS get_control_points(Eigen::Matrix<T, point_size, Eigen::Dynamic> &points) const
        {
            points = m_points;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS point_on_curve(T u, Eigen::Vector<T, dim> &point) const
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
        int get_degree() const
        {
            return m_degree;
        }
    };

    template<typename T, int dim, bool is_rational, int points_count, int degree>
    struct geo_traits<polynomial_curve<T, dim, is_rational, points_count, degree> >
    {
        static constexpr int point_size = is_rational ? dim + 1 : dim;
        using type = polynomial_curve<T, dim, is_rational, points_count, degree>;
        using point_number_type = T;
        using point_type = typename  Eigen::Vector<T, dim> ;
    };

    template<typename T, int dim, bool is_rational>
    ENUM_NURBS save_obj(T low, T high, const polynomial_curve<T, dim, is_rational, -1, -1> &polynomial_cur, const std::string &path)
    {
        T step = (high - low) / 100.0;
        std::vector<Eigen::Vector<T, dim>> points;
        for (int u_index = 0; u_index < 100; ++u_index)
        {
            Eigen::Vector<T, dim> point;
            polynomial_cur.point_on_curve(u_index * step + low, point);
            points.push_back(point);
        }
        std::ofstream outfile2(path);
        for (auto point : points)
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }



}
