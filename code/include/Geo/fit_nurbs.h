#pragma once

#include <concepts>
#include "nurbs_curve.h"
#include "nurbs_surface.h"
#include "create_nurbs_arc.h"
#include <Eigen/Sparse>

enum ENPARAMETERIEDTYPE
{
    CHORD = 0,
    CENTRIPETAL
};
template<typename T, int dim>
ENUM_NURBS global_curve_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, int degree, 
    const std::vector<T> &params, nurbs_curve<T, dim, false, -1, -1> &nurbs)
{
    int params_size = params.size();

    Eigen::VectorX<T> knots_vector(params_size + degree + 1);
    knots_vector.block(0, 0, degree + 1, 1).setConstant(0.0);
    knots_vector.block(params_size, 0, degree + 1, 1).setConstant(1.0);
    int internel_knots_count = params_size - degree - 1;
    T inverse_degree = 1.0 / static_cast<T>(degree);
    auto it = params.begin() + 1;
    auto it2 = it + degree;
    
    for (int j = 1; j <= internel_knots_count; ++j)
    {
        //求和可以优化
        knots_vector[degree + j] = std::accumulate(it, it2, 0.0) * inverse_degree;
        ++it;
        ++it2;
    }

    
    int knots_size = knots_vector.size();
    if (params_size != knots_size - degree - 1)
        return ENUM_NURBS::NURBS_ERROR;
    Eigen::SparseMatrix<T> mat(params_size, params_size);
    mat.reserve(Eigen::VectorXi::Constant(params_size, degree + 1));
    mat.insert(0, 0) = 1.0;
    mat.insert(params_size - 1, params_size - 1) = 1.0;
    for (int i = 1; i < degree; ++i)
    {
        mat.insert(0, i) = 0.0;
        mat.insert(params_size - 1, params_size - 1 - i);
    }
    for (int row = 1; row < params_size - 1; ++row)
    {
        int index = -1;
        find_span<T>(params[row], degree, knots_vector, index);
        Eigen::VectorX<T> basis;
        basis_functions<T>(index, params[row], degree, knots_vector, basis);
        for (int i = 0; i <= degree; ++i)
        {
            mat.insert(row, index + i - degree) = basis[i];
        }
    }
    mat.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
    solver.compute(mat);
    if (solver.info() != Eigen::Success)
        return ENUM_NURBS::NURBS_ERROR;
    Eigen::Matrix<T, Eigen::Dynamic, dim> control_points = solver.solve(points.transpose());
    if (solver.info() != Eigen::Success)
        return ENUM_NURBS::NURBS_ERROR;
    Eigen::Matrix<T, dim, Eigen::Dynamic> transpose_points = control_points.transpose();
    nurbs.set_control_points(transpose_points);
    nurbs.set_knots_vector(knots_vector);
    nurbs.set_degree(degree);
    return ENUM_NURBS::NURBS_SUCCESS;
}

template<typename T, int dim, int parameteried_type>
ENUM_NURBS global_curve_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, int degree, 
        nurbs_curve<T, dim, false, -1, -1> &nurbs)
{
    int points_count = points.cols();
    std::vector<T> params(points_count);
    if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
    {
        T d = 0.0;
        std::vector<T> vec_norm(points_count - 1);
        for (int index = 1; index < points_count; ++index)
        {
            vec_norm[index - 1] = (points.col(index) - points.col(index - 1)).norm();
            d += vec_norm[index - 1];
        }
           
        params[0] = 0.0;
        //可以优化以下参数的取值范围?
        params[points_count - 1] = 1.0;
        for (int index = 1; index < points_count - 1; ++index)
        {
            params[index] = params[index - 1] + vec_norm[index - 1] / d;
        }
            
    }
    else if constexpr(parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
    {
        T d = 0.0;
        std::vector<T> vec_norm(points - 1);
        for (int index = 1; index < points_count; ++index)
        {
            vec_norm[index - 1] = std::sqrt((points.col(index) - points.col(index - 1)).norm());
            d += vec_norm[index - 1];
        }
           
        params[0] = 0.0;
        //可以优化以下参数的取值范围?
        params[points_count - 1] = 1.0;
        for (int index = 1; index < points_count - 1; ++index)
            params[index] = params[index - 1] + vec_norm[index - 1] / d;
    }
    else
    {
        return ENUM_NURBS::NURBS_ERROR;
    }

    return global_curve_interpolate<T, dim>(points, degree, params, nurbs);
}