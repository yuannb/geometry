#pragma once
#include "nurbs_tool.h"
template<typename T, int dim, bool is_rational, int points_count>
class bezier_curve
{
private:
    static constexpr int row_size = is_rational ? dim + 1 : dim;
    Eigen::Matrix<T, row_size, points_count> m_control_points;
public:
    bezier_curve() = default;
    bezier_curve(Eigen::Vector<Eigen::Vector<T, row_size>, points_count> const &points)
    {
        for (int col_index = 0; col_index < points_count; ++col_index)
        {
            m_control_points.col(col_index) = points[col_index];
        }
    }
    ENUM_NURBS point_on_curve(T u, Eigen::Vector<T, dim> &point) const
    {
        if (u < 0.0 || u > 1.0)
            return NURBS_PARAM_IS_OUT_OF_DOMAIN;
        Eigen::Vector<T, row_size> vec;
        DeCasteljaul<T, dim, points_count, is_rational>(m_control_points, u, vec);
        point = project_point<T, is_rational, row_size>::project_point_to_euclidean_space(vec);
        return NURBS_SUCCESS;
    }

    ENUM_NURBS point_on_curve_b(T u, Eigen::Vector<T, dim> &point) const
    {
        if (u < 0.0 || u > 1.0)
            return NURBS_PARAM_IS_OUT_OF_DOMAIN;
        Eigen::Vector<T, points_count> basis;
        AllBernstein<T, points_count>(u, basis);
        Eigen::Vector<T, row_size> vec = m_control_points * basis;
        point = project_point<T, is_rational, row_size>::project_point_to_euclidean_space(vec);
        return NURBS_SUCCESS;
    }

    ENUM_NURBS set_control_points(const Eigen::Matrix<T, row_size, points_count> &control_points)
    {
        m_control_points = control_points;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    ENUM_NURBS get_control_points(Eigen::Matrix<T, row_size, points_count> &control_points) const
    {
        control_points = m_control_points;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    
    ENUM_NURBS get_control_points(int index, Eigen::Vector<T, dim> &control_point) const
    {
        if (index >= points_count || index < 0)
            return ENUM_NURBS::NURBS_ERROR;
        if constexpr (is_rational == false)
        {
            control_point = m_control_points.col(index);
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        control_point = m_control_points.template block<dim, 1>(0, index) / m_control_points(dim, index);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS get_weight(int index, T &w) const
    {
        if (index >= points_count || index < 0)
            return ENUM_NURBS::NURBS_ERROR;        
        if constexpr (is_rational == false)
        {
            w = 1.0;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        w = m_control_points(dim, index);
        return ENUM_NURBS::NURBS_SUCCESS;
    }


};

template<typename T, int dim, bool is_rational>
class bezier_curve<T, dim, is_rational, -1>
{
private:
    static constexpr int row_size = is_rational ? dim + 1 : dim;
    Eigen::Matrix<T, row_size, Eigen::Dynamic> m_control_points;
public:
    bezier_curve() = default;
    bezier_curve(std::vector<Eigen::Vector<T, row_size>> const &points)
    {
        size_t col_size = points.size(); 
        m_control_points.resize(row_size, col_size);
        for (size_t col_index = 0; col_index < col_size; ++col_index)
        {
            m_control_points.col(col_index) = points[col_index];
        }
    }
    ENUM_NURBS point_on_curve(T u, Eigen::Vector<T, dim> &point) const
    {
        if (u < 0.0 || u > 1.0)
            return NURBS_PARAM_IS_OUT_OF_DOMAIN;
        Eigen::Vector<T, row_size> vec;
        DeCasteljaul<T, dim, is_rational>(m_control_points, m_control_points.cols(), u, vec);
        point = project_point<T, is_rational, row_size>::project_point_to_euclidean_space(vec);
        return NURBS_SUCCESS;
    }

    ENUM_NURBS point_on_curve_b(T u, Eigen::Vector<T, dim> &point) const
    {
        if (u < 0.0 || u > 1.0)
            return NURBS_PARAM_IS_OUT_OF_DOMAIN;
        Eigen::Vector<T, Eigen::Dynamic> basis;
        AllBernstein<T>(m_control_points.rows(), u, basis);
        Eigen::Vector<T, row_size> vec = m_control_points * basis;
        point = project_point<T, is_rational, row_size>::project_point_to_euclidean_space(vec);
        return NURBS_SUCCESS;
    }

    ENUM_NURBS set_control_points(const Eigen::Matrix<T, row_size, Eigen::Dynamic> &control_points)
    {
        m_control_points = control_points;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    ENUM_NURBS get_control_points(Eigen::Matrix<T, row_size, Eigen::Dynamic> &control_points) const
    {
        control_points = m_control_points;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    ENUM_NURBS get_control_points(int index, Eigen::Vector<T, dim> &control_point) const
    {
        if (index >= m_control_points.cols() || index < 0)
            return ENUM_NURBS::NURBS_ERROR;
        if constexpr (is_rational == false)
        {
            control_point = m_control_points.col(index);
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        control_point = m_control_points.template block<dim, 1>(0, index) / m_control_points(dim, index);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS get_weight(int index, T &w) const
    {
        if (index >= m_control_points.cols() || index < 0)
            return ENUM_NURBS::NURBS_ERROR;        
        if constexpr (is_rational == false)
        {
            w = 1.0;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        w = m_control_points(dim, index);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
};
