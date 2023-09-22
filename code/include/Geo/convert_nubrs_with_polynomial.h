#pragma once

#include "nurbs_curve.h"
#include "polynomial_curve.h"
#include "nurbs_surface.h"
#include "polynomial_surface.h"
//将转换函数新写一个文件可以避免nurbs和polynomail的相互包含
template<typename T, int dim, bool is_rational, int points_count, int degree>
ENUM_NURBS convert_nubrs_to_polynomial(const nurbs_curve<T, dim, is_rational, points_count, degree> &old_curve,
    std::vector<polynomial_curve<T, dim, is_rational, points_count, degree>*> &polynomial_curves)
{
    std::vector<nurbs_curve<T, dim, is_rational, -1, -1> *> bezier_curves;
    old_curve.decompose_to_bezier(bezier_curves);
    int curves_count = bezier_curves.size();
    int curve_degree = old_curve.get_degree();
    polynomial_curves.resize(curves_count, nullptr);
    Eigen::MatrixX<T> M_p;
    bezier_to_power_matrix(curve_degree, M_p);
    for (int index = 0; index < curves_count; ++index)
    {
        Eigen::MatrixX<T> R_p;
        std::array<T, 2> ends_knots;
        bezier_curves[index]->get_ends_knots(ends_knots);
        T c = 1.0 / (ends_knots[1] - ends_knots[0]);
        T d = -ends_knots[0] * c;
        reparameter_matrix_of_p_degree(curve_degree, c, d, R_p);
        Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic> new_points = bezier_curves[index]->get_control_points() * M_p * R_p;
        polynomial_curve<T, dim, is_rational, -1, -1> *new_polynomial_curve = new polynomial_curve<T, dim, is_rational, -1, -1>(new_points);
        delete bezier_curves[index];
        polynomial_curves[index] = new_polynomial_curve;
    
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}

template<typename T, int dim, int rows, int cols, int u_degree, int v_degree, bool is_rational>
ENUM_NURBS convert_nubrs_to_polynomial(const nurbs_surface<T, dim, rows, cols, u_degree, v_degree, is_rational> &old_surface,
    Eigen::MatrixX<polynomial_surface<T, dim, u_degree, v_degree, is_rational>*> &polynomial_surfaces)
{
    Eigen::MatrixX<nurbs_surface<T, dim, rows, cols, u_degree, v_degree, is_rational> *> bezier_surfaces;
    old_surface.decompose_to_bezier(bezier_surfaces);
    int new_rows = bezier_surfaces.rows();
    int new_cols = bezier_surfaces.cols();
    polynomial_surfaces.resize(new_rows, new_cols);
    int surface_u_degree = old_surface.get_u_degree();
    int surface_v_degree = old_surface.get_v_degree();

    Eigen::MatrixX<T> M_p, M_q;
    bezier_to_power_matrix(surface_u_degree, M_p);
    bezier_to_power_matrix(surface_v_degree, M_q);
    M_q.transposeInPlace();
    std::vector<T> u_different_knots, v_different_knots;
    old_surface.get_u_different_knots(u_different_knots);
    old_surface.get_v_different_knots(v_different_knots);
    std::vector<Eigen::MatrixX<T>> R_ps(new_rows);

    for (int row_index = 1; row_index <= new_rows; ++row_index)
    {
        T c = 1.0 / (u_different_knots[row_index] - u_different_knots[row_index - 1]);
        T d = -u_different_knots[row_index - 1] * c;
        Eigen::MatrixX<T> R_p;
        reparameter_matrix_of_p_degree(surface_u_degree, c, d, R_p);
        R_ps[row_index - 1] = R_p;
    }
    for (int col_index = 1; col_index <= new_cols; ++col_index)
    {
        Eigen::MatrixX<T> R_q;
        T c = 1.0 / (v_different_knots[col_index] - v_different_knots[col_index - 1]);
        T d = -v_different_knots[col_index - 1] * c;
        reparameter_matrix_of_p_degree(surface_v_degree, c, d, R_q);
        R_q.transposeInPlace();
        std::vector<Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic>> temp(surface_v_degree + 1);
        Eigen::VectorX<Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic>> polynomial_points(surface_v_degree + 1);
        for (int row_index = 1; row_index <= new_rows; ++row_index)
        {
            //列优先存储, 所以先循环列, 再循环行
            auto* surf = bezier_surfaces(row_index - 1, col_index - 1);
            for (int index = 0; index <= surface_v_degree; ++index)
            {
                Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic> row_points;
                surf->get_control_points_row(index, row_points);
                temp[index] = row_points * M_p * R_ps[row_index - 1];
            }

            auto left_mat = R_q * M_q;

            for (int v_index = 0; v_index <= surface_v_degree; ++v_index)
            {
                polynomial_points[v_index].resize(is_rational ? dim + 1 : dim, surface_u_degree + 1);
                polynomial_points[v_index].setConstant(0.0);
                for (int u_index = 0; u_index <= surface_u_degree; ++u_index)
                {
                    for (int index = 0; index <= surface_v_degree; ++index)
                    {
                        polynomial_points[v_index].col(u_index) += left_mat(v_index, index) * temp[index].col(u_index);
                    }
                }
            }
            polynomial_surfaces(row_index - 1, col_index - 1) = new polynomial_surface<T, dim, -1, -1, is_rational>(polynomial_points);
            delete surf;
        }
    }

    return ENUM_NURBS::NURBS_SUCCESS;
}


//此处为单polynomial curve(surface) 转nurbs curve(surface). 多条polynomial的问题可以将每条polynomial转成nurbs
//然后调用nurbs的相关算法拼成一个, 然后再调用节点消去算法, 最终可以得到一条nurbs
template<typename T, int dim, bool is_rational, int points_count, int degree>
ENUM_NURBS convert_polynomial_to_nubrs(const polynomial_curve<T, dim, is_rational, points_count, degree> &old_curve, const Interval<T> &interval,
    nurbs_curve<T, dim, is_rational, points_count, degree> &nurbs)
{
    Eigen::MatrixX<T> M_p;
    int curve_degree = old_curve.get_degree();
    bezier_to_power_matrix(curve_degree, M_p);
    Eigen::MatrixX<T> MI_p;
    power_matrix_to_bezier(M_p, MI_p);

    T high = interval.get_high();
    T low = interval.get_low();
    T c = high - low;

    Eigen::MatrixX<T> R_p;
    reparameter_matrix_of_p_degree(curve_degree, c, low, R_p);

    Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic> control_points;
    old_curve.get_control_points(control_points);
    Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic> new_control_points = control_points * (R_p * MI_p);

    Eigen::VectorX<T> knots_vector(2 * curve_degree + 2);
    knots_vector.block(0, 0, curve_degree + 1, 1).setConstant(low);
    knots_vector.block(curve_degree + 1, 0, curve_degree + 1, 1).setConstant(high);
    nurbs.set_control_points(new_control_points);
    nurbs.set_knots_vector(knots_vector);
    nurbs.set_degree(curve_degree);

    return ENUM_NURBS::NURBS_SUCCESS;
}


template<typename T, int dim, int rows, int cols, int u_degree, int v_degree, bool is_rational>
ENUM_NURBS convert_polynomial_to_nubrs(const polynomial_surface<T, dim, u_degree, v_degree, is_rational> &old_surface,
    const Interval<T> &u_interval, const Interval<T> &v_interval,  nurbs_surface<T, dim, rows, cols, u_degree, v_degree, is_rational> &nurbs)
{
    int surface_u_degree = old_surface.get_u_degree();
    int surface_v_degree = old_surface.get_v_degree();

    Eigen::MatrixX<T> M_p, M_q;
    bezier_to_power_matrix(surface_u_degree, M_p);
    bezier_to_power_matrix(surface_v_degree, M_q);
    Eigen::MatrixX<T> MI_p, MI_q;
    power_matrix_to_bezier(M_p, MI_p);
    power_matrix_to_bezier(M_q, MI_q);
    MI_q.transposeInPlace();
    Eigen::MatrixX<T> R_p;
    T u_low = u_interval.get_low();
    T u_high = u_interval.get_high();
    reparameter_matrix_of_p_degree(surface_u_degree, u_high - u_low, u_low, R_p);

    Eigen::MatrixX<T> R_q;
    T v_low = v_interval.get_low();
    T v_high = v_interval.get_high();
    reparameter_matrix_of_p_degree(surface_v_degree, v_high - v_low, v_low, R_q);
    R_q.transposeInPlace();
    std::vector<Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic>> temp(surface_v_degree + 1);
    Eigen::VectorX<Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic>> nurbs_points(surface_v_degree + 1);

    for (int index = 0; index <= surface_v_degree; ++index)
    {
        Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic> row_points;
        old_surface.get_control_points_row(index, row_points);
        temp[index] = row_points * R_p * MI_p;
    }

    auto left_mat =  MI_q * R_q;

    for (int v_index = 0; v_index <= surface_v_degree; ++v_index)
    {
        nurbs_points[v_index].resize(is_rational ? dim + 1 : dim, surface_u_degree + 1);
        nurbs_points[v_index].setConstant(0.0);
        for (int u_index = 0; u_index <= surface_u_degree; ++u_index)
        {
            for (int index = 0; index <= surface_v_degree; ++index)
            {
                nurbs_points[v_index].col(u_index) += left_mat(v_index, index) * temp[index].col(u_index);
            }
        }
    }
    Eigen::VectorX<T> u_knots_vector(2 * surface_u_degree + 2);
    Eigen::VectorX<T> v_knots_vector(2 * surface_v_degree + 2);
    u_knots_vector.block(0, 0, surface_u_degree + 1, 1).setConstant(u_low);
    u_knots_vector.block(surface_u_degree + 1, 0, surface_u_degree + 1, 1).setConstant(u_high);
    v_knots_vector.block(0, 0, surface_v_degree + 1, 1).setConstant(v_low);
    v_knots_vector.block(surface_v_degree + 1, 0, surface_v_degree + 1, 1).setConstant(v_high);
    nurbs.set_control_points(nurbs_points);
    nurbs.set_uv_degree(surface_u_degree, surface_v_degree);
    nurbs.set_uv_knots(u_knots_vector, v_knots_vector);

    return ENUM_NURBS::NURBS_SUCCESS;
}


