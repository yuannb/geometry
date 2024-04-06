#pragma once
#include "nurbs_tool.h"
#include "bezier_curve.h"

namespace tnurbs
{
    using namespace tnurbs;
    //TODO: U向控制点或者V向控制点个数不固定, 另一个方向固定

    /// @brief bezier curve class
    /// @tparam T : double float int
    /// @tparam dim : nurbs所在欧氏空间的维数
    /// @tparam row_count : U向控制点的个数(不严谨)
    /// @tparam col_count : V向控制点的个数(不严谨)
    /// @tparam is_rational : 是否时有理nurbs
    template<typename T, int dim, int row_count, int col_count, bool is_rational>
    class bezier_surface
    {
    private:
        constexpr static int row_size = is_rational ? dim + 1: dim;
        Eigen::Matrix<Eigen::Vector<T, row_size>, row_count, col_count> m_control_points; //p_ij(i : row, j : col)
    public:
        bezier_surface(Eigen::Matrix<Eigen::Vector<T, row_size>, row_count, col_count> &points)
        {
            m_control_points = points;
        }
        bezier_surface(const std::vector<Eigen::Vector<T, row_size>> &points)
        {
            int points_counts = points.size();
            assert(points_counts == row_count * col_count);
            int index = 0;
            for (int row_index = 0; row_index < row_count; ++row_index)
            {
                for (int col_index = 0; col_index < col_count; ++col_index)
                    m_control_points(row_index, col_index) = points[index++];
            }
        }
        ENUM_NURBS point_on_surface(T u,T v, Eigen::Vector<T, dim> &point) const
        {
            if (u < 0.0 || u > 1.0)
                return NURBS_PARAM_IS_OUT_OF_DOMAIN;
            if (v < 0.0 || v > 1.0)
                return NURBS_PARAM_IS_OUT_OF_DOMAIN;
            Eigen::Vector<T, row_size> vec;
            DeCasteljaul2<T, dim, is_rational, row_count, col_count>(m_control_points, u, v, vec);
            point = project_point<T, is_rational, row_size>::project_point_to_euclidean_space(vec);
            return NURBS_SUCCESS;
        }

        ENUM_NURBS get_control_points(Eigen::MatrixX<Eigen::Vector<T, row_size>>& control_points) const;
    
    };

    /// @brief bezier curve class
    /// @tparam T : double float int
    /// @tparam dim : nurbs所在欧氏空间的维数
    /// @tparam is_rational : 是否时有理nurbs
    template<typename T, int dim, bool is_rational>
    class bezier_surface<T, dim, -1, -1, is_rational>
    {
    private:
        constexpr static int row_size = is_rational ? dim + 1: dim;
        Eigen::Matrix<Eigen::Vector<T, row_size>, Eigen::Dynamic, Eigen::Dynamic> m_control_points; //p_ij(i : row, j : col)
    public:

        using iso_curve_type = bezier_curve<T, dim, is_rational, -1>;

        bezier_surface(int rows, int cols, std::vector<Eigen::Vector<T, row_size>> &points)
        {
            m_control_points.resize(rows, cols);
            int index = 0;
            for (int row_index = 0; row_index < rows; ++row_index)
            {
                for (int col_index = 0; col_index < cols; ++col_index)
                    m_control_points(row_index, col_index) = points[index++];
            }
        }
        bezier_surface(Eigen::MatrixX<Eigen::Vector<T, row_size>> &points)
        {
            m_control_points = points;
        }

        bezier_surface(const Eigen::VectorX<Eigen::Matrix<T, row_size, Eigen::Dynamic>>& points)
        {
            int rows = points.rows();
            int cols = points[0].cols();
            m_control_points.resize(cols, rows);
            for (int v_index = 0; v_index < rows; ++v_index)
            {
                for (int u_index = 0; u_index < cols; ++u_index)
                {
                    m_control_points(u_index, v_index) = points[v_index].col(u_index);
                }
            }
        }
        
        ENUM_NURBS point_on_surface(T u,T v, Eigen::Vector<T, dim> &point)
        {
            if (u < 0.0 || u > 1.0)
                return NURBS_PARAM_IS_OUT_OF_DOMAIN;
            if (v < 0.0 || v > 1.0)
                return NURBS_PARAM_IS_OUT_OF_DOMAIN;
            Eigen::Vector<T, row_size> vec;
            int rows = m_control_points.rows();
            int cols = m_control_points.cols();
            DeCasteljaul2<T, dim, is_rational>(m_control_points, rows, cols, u, v, vec);
            point = project_point<T, is_rational, row_size>::project_point_to_euclidean_space(vec);
            return NURBS_SUCCESS;
        }
    
        ENUM_NURBS get_control_points(Eigen::MatrixX<Eigen::Vector<T, row_size>>& control_points) const
        {
            control_points = m_control_points;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
    
    };


}

