#pragma once

#include <concepts>
#include "nurbs_curve.h"
#include "nurbs_surface.h"
#include "create_nurbs_arc.h"

template<typename T, int dim>
// requires (dim >= 2)
ENUM_NURBS create_bilinear_nurbs_surface(const Eigen::Vector<T, dim> &P00, const Eigen::Vector<T, dim> &P01, const Eigen::Vector<T, dim> &P10,
    const Eigen::Vector<T, dim> &P11, nurbs_surface<T, dim, -1, -1, -1, -1, false> &nurbs)
{
    Eigen::Vector<Eigen::Matrix<T, dim, Eigen::Dynamic>, Eigen::Dynamic> control_points(2);
    Eigen::Matrix<T, dim, Eigen::Dynamic> fisrt_row(dim, 2);
    fisrt_row << P00, P01;
    control_points[0] = fisrt_row;
    Eigen::Matrix<T, dim, Eigen::Dynamic> second_row(dim, 2);
    second_row << P10, P11;
    control_points[1] = second_row;
    Eigen::VectorX<T> knots_vector(4);
    knots_vector << 0, 0, 1, 1;
    nurbs.set_control_points(control_points);
    nurbs.set_uv_degree(1, 1);
    nurbs.set_uv_knots(knots_vector, knots_vector);
    return ENUM_NURBS::NURBS_SUCCESS; 
}


//dim >= 2
template<typename T, int dim = 3, bool is_rational>
// requires (dim >= 2)
ENUM_NURBS create_cylinder_surface(const nurbs_curve<T, dim, is_rational, -1, -1> &section_curve, T d, const Eigen::Vector<T, dim> &direction, nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &surf)
{
    constexpr int point_size = is_rational ? dim + 1 : dim;
    Eigen::VectorX<T> u_knots_vector = section_curve.get_knots_vector();
    Eigen::VectorX<T> v_knots_vector(4);
    v_knots_vector << 0, 0, 1, 1;
    Eigen::Matrix<T, point_size, Eigen::Dynamic> control_points = section_curve.get_control_points();
    Eigen::Matrix<T, point_size, Eigen::Dynamic> control_points2;
    if constexpr (is_rational == false)
    {
        control_points2 = control_points2.rowwise() + d * direction;
    }
    else
    {
        Eigen::Vector<T, dim> translate_vector = d * direction;
        int points_count = control_points.cols();
        control_points2.resize(point_size, points_count);
        for (int index = 0; index < points_count; ++index)
        {
            control_points2.template block<dim, 1>(0, index) = control_points.template block<dim, 1>(0, index) + translate_vector * control_points(dim, index);
            control_points2(dim, index) = control_points(dim, index);
        }
    }
    Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> surface_control_points(2);
    surface_control_points[0] = control_points;
    surface_control_points[1] = control_points2;
    surf.set_uv_degree(section_curve.get_degree(), 1);
    surf.set_control_points(surface_control_points);
    surf.set_uv_knots(u_knots_vector, v_knots_vector);
    return ENUM_NURBS::NURBS_SUCCESS;
}


//dim >= 2
//为了简化此函数, 要求两条曲线全是有理nurbs曲线, 不然内部需要处理一条曲线为有理和一条曲线为无理的繁琐工作
template<typename T, int dim = 3>
// requires (dim >= 2)
ENUM_NURBS create_ruled_surface(nurbs_curve<T, dim, true, -1, -1> section_curve1, nurbs_curve<T, dim, true, -1, -1> section_curve2, nurbs_surface<T, dim, -1, -1, -1, -1, true> &surf)
{
    
    int degree1 = section_curve1.get_degree();
    int degree2 = section_curve2.get_degree();
    int degree = std::max(degree1, degree2);
    if (degree > degree1)
    {
        section_curve1.degree_elevate(degree - degree1);

    }
    else if (degree > degree2)
    {
        section_curve2.degree_elevate(degree - degree2);
    }


    //将curve2的区间与curve1对齐
    Eigen::VectorX<T> knots_vector1 =  section_curve1.get_knots_vector();
    Eigen::VectorX<T> knots_vector2 =  section_curve2.get_knots_vector();
    int knots_vector1_size = knots_vector1.size();
    int knots_vector2_size = knots_vector2.size();
    T low1 = knots_vector1[0];
    T high1 = knots_vector1[knots_vector1_size - 1];
    T low2 = knots_vector2[0];
    T high2 = knots_vector2[knots_vector2_size - 1];
    T alpha = (high2 - low2) / (high1 - low1);
    T beta = low2 - alpha * low1;
    nurbs_curve<T, dim, true, -1, -1> new_section_curve2;
    section_curve2.curve_reparameter_with_linear_function(alpha, beta, new_section_curve2);

    knots_vector2 = new_section_curve2.get_knots_vector();
    std::vector<T> merge_knots_vector;
    std::vector<T> knots_vector1_add;
    std::vector<T> knots_vector2_add;

    merge_two_knots_vector<T>(degree, knots_vector1, knots_vector2, merge_knots_vector, knots_vector1_add, knots_vector2_add);
    
    int konts_vector1_add_size = knots_vector1_add.size();
    int konts_vector2_add_size = knots_vector2_add.size();
    int merge_knots_vector_vec_size = merge_knots_vector.size();
    Eigen::VectorX<T> knots_vector1_add_vec(konts_vector1_add_size);
    Eigen::VectorX<T> knots_vector2_add_vec(konts_vector2_add_size);
    Eigen::VectorX<T> merge_knots_vector_vec(merge_knots_vector_vec_size);
    for (int index = 0; index < konts_vector1_add_size; ++index)
    {
        knots_vector1_add_vec[index] = knots_vector1_add[index];
    }
    for (int index = 0; index < konts_vector2_add_size; ++index)
    {
        knots_vector2_add_vec[index] = knots_vector2_add[index];
    }
    for (int index = 0; index < merge_knots_vector_vec_size; ++index)
    {
        merge_knots_vector_vec[index] = merge_knots_vector[index];
    }

    
    section_curve1.refine_knots_vector(knots_vector1_add_vec);
    section_curve2.refine_knots_vector(knots_vector2_add_vec);

    Eigen::Matrix<T, dim + 1, Eigen::Dynamic> control_points1 = section_curve1.get_control_points();
    Eigen::Matrix<T, dim + 1, Eigen::Dynamic> control_points2 = section_curve2.get_control_points();
    Eigen::VectorX<Eigen::Matrix<T, dim + 1, Eigen::Dynamic>> surface_control_points(2);
    surface_control_points[0] = control_points1;
    surface_control_points[1] = control_points2;
    surf.set_uv_degree(degree, 1);
    surf.set_control_points(surface_control_points);
    Eigen::VectorX<T> v_knots_vector(4);
    v_knots_vector << 0, 0, 1, 1;
    surf.set_uv_knots(merge_knots_vector_vec, v_knots_vector);
    return ENUM_NURBS::NURBS_SUCCESS;
}

//为了简化此函数, 要求曲线是有理nurbs曲线, 不然内部需要处理曲线为无理的繁琐工作
template<typename T>
ENUM_NURBS create_revolved_surface(const nurbs_curve<T, 3, true, -1, -1> &v_curve, const Eigen::Vector<T, 3> &origin, const Eigen::Vector<T, 3> &aix_direction, T theta,
    nurbs_surface<T, 3, -1, -1, -1, -1, true> &surf)
{
    int narcs = std::ceil(theta / M_PI_2);
    T delta_theta = theta / static_cast<T>(narcs);
    std::vector<T> cosines(narcs), sines(narcs);
    T angle = delta_theta;
    for (int i = 0; i < narcs; ++i)
    {
        cosines[i] = std::cos(angle);
        sines[i] = std::sin(angle);
        angle += delta_theta;
    }

    int rows = v_curve.get_control_points_count();
    Eigen::Vector<T, 3> unit_direction = aix_direction.normalized();
    Eigen::VectorX<Eigen::Matrix<T, 4, Eigen::Dynamic>> control_points(rows);
    Eigen::VectorX<T> u_knots_vector;
    for (int row_index = 0; row_index < rows; ++row_index)
    {
        Eigen::Matrix<T, 4, Eigen::Dynamic> contro_points_row(4, 2 * narcs + 1);
        Eigen::Vector<T, 3> point;
        v_curve.get_control_point(row_index, point);
        T weight;
        v_curve.get_weight(row_index, weight);

        Eigen::Vector<T, 3> project_point = origin + (point - origin).dot(unit_direction) * unit_direction;
        Eigen::Vector<T, 3> X = point - project_point;
        T r = X.norm();
        X.normalize();
        Eigen::Vector3<T> Y = unit_direction.cross(X);
        nurbs_curve<T, 3, true, -1, -1> circle;
        create_nurbs_circle<T, 3>(project_point, X, Y, r, 0, theta, circle);
        control_points[row_index] = weight * circle.get_control_points();
        if (row_index == 0)
            u_knots_vector = circle.get_knots_vector();
    }
    surf.set_control_points(control_points);
    surf.set_uv_knots(u_knots_vector, v_curve.get_knots_vector());
    surf.set_uv_degree(2, v_curve.get_degree());

    return ENUM_NURBS::NURBS_SUCCESS;
}

