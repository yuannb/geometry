#pragma once

#include <concepts>
#include "nurbs_curve.h"
#include "nurbs_surface.h"
#include "create_nurbs_arc.h"
#include <Eigen/Sparse>

namespace tnurbs
{
    using namespace tnurbs;
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
        
        for (int j = 1; j <= internel_knots_count; ++j)
        {
            //求和可以优化
            knots_vector[degree + j] = std::accumulate(it, it + degree, 0.0) * inverse_degree;
            ++it;
        }

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


    template<typename T, int dim>
    ENUM_NURBS global_curve_interpolate_with_ends_tangent(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Vector<T, dim> &D0,
       const Eigen::Vector<T, dim> &D1, int degree, const std::vector<T> &params, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        int params_size = params.size();

        Eigen::VectorX<T> knots_vector(params_size + degree + 3);
        knots_vector.block(0, 0, degree + 1, 1).setConstant(0.0);
        knots_vector.block(params_size + 2, 0, degree + 1, 1).setConstant(1.0);
        int internel_knots_count = params_size - degree;
        T inverse_degree = 1.0 / static_cast<T>(degree);
        auto it = params.begin();
        
        for (int j = 0; j <= internel_knots_count; ++j)
        {
            //求和可以优化
            knots_vector[degree + j + 1] = std::accumulate(it, it + degree, 0.0) * inverse_degree;
            ++it;
        }

        
        Eigen::SparseMatrix<T> mat(params_size + 2, params_size + 2);
        mat.reserve(Eigen::VectorXi::Constant(params_size + 2, degree + 1));
        mat.insert(0, 0) = 1.0;
        mat.insert(params_size + 1, params_size + 1) = 1.0;
        mat.insert(1, 0) = -1.0;
        mat.insert(1, 1) = 1.0;
        mat.insert(params_size, params_size) = -1.0;
        mat.insert(params_size, params_size + 1) = 1.0;
        for (int i = 1; i < degree; ++i)
        {
            mat.insert(0, i) = 0.0;
            mat.insert(params_size + 1, params_size + 1 - i);
        }
        for (int  i = 2; i < degree; ++i)
        {
            mat.insert(1, i) = 0.0;
            mat.insert(params_size, params_size - i + 1) = 0.0;
        }

        for (int row = 2; row < params_size; ++row)
        {
            int index = -1;
            find_span<T>(params[row - 1], degree, knots_vector, index);
            Eigen::VectorX<T> basis;
            basis_functions<T>(index, params[row - 1], degree, knots_vector, basis);
            for (int i = 0; i <= degree; ++i)
            {
                mat.insert(row, index + i - degree) = basis[i];
            }
        }

        Eigen::Matrix<T, Eigen::Dynamic, dim> new_points(params_size + 2, dim);
        new_points.template block<1, dim>(0, 0) = (points.col(0)).transpose();
        new_points.template block<1, dim>(1, 0) = D0.transpose() * (knots_vector[degree + 1] * inverse_degree);
        new_points.template block(2, 0, params_size - 2, dim) = (points.block(0, 1, dim, params_size - 2)).transpose();
        new_points.template block<1, dim>(params_size, 0) = D1.transpose() * ((1.0 - knots_vector[params_size + 1]) * inverse_degree);
        new_points.template block<1, dim>(params_size + 1, 0) = (points.col(params_size - 1)).transpose();

        mat.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
        solver.compute(mat);
        if (solver.info() != Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        Eigen::Matrix<T, Eigen::Dynamic, dim> control_points = solver.solve(new_points);
        if (solver.info() != Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        Eigen::Matrix<T, dim, Eigen::Dynamic> transpose_points = control_points.transpose();
        nurbs.set_control_points(transpose_points);
        nurbs.set_knots_vector(knots_vector);
        nurbs.set_degree(degree);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim, int parameteried_type>
    ENUM_NURBS global_curve_interpolate_with_ends_tangent(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Vector<T, dim> &D0, 
        const Eigen::Vector<T, dim> &D1, int degree,  nurbs_curve<T, dim, false, -1, -1> &nurbs)
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

        return global_curve_interpolate_with_ends_tangent<T, dim>(points, D0, D1, degree, params, nurbs);
    }


    template<typename T, int dim>
    ENUM_NURBS global_3degree_curve_interpolate_with_ends_tangent(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Vector<T, dim> &D0, 
        const Eigen::Vector<T, dim> &D1, const std::vector<T> &params, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        
        int points_count = points.cols();
        int params_count = params.size();
        if (points_count != params_count)
        {
            return ENUM_NURBS::NURBS_ERROR;
        }
        Eigen::VectorX<T> knots_vector(params_count + 6);
        knots_vector.template block<4, 1>(0, 0).setConstant(0.0);
        knots_vector.template block<4, 1>(params_count + 2, 0).setConstant(1.0);
        for (int i = 1; i < params_count - 1; ++i)
        {
            knots_vector[i + 3] = params[i];
        }
        Eigen::Matrix<T, dim, Eigen::Dynamic> control_points(dim, points_count + 2);
        control_points.col(0) = points.col(0);
        control_points.col(1) = (knots_vector[4] / 3.0) * D0 + control_points.col(0);
        control_points.col(points_count + 1) = points.col(points_count - 1);
        control_points.col(points_count) = control_points.col(points_count + 1) - ((1 - knots_vector[points_count + 1]) / 3.0) * D1;
        solve_tri_diagonal<T, dim>(points, knots_vector, control_points);

        nurbs.set_control_points(control_points);
        nurbs.set_knots_vector(knots_vector);
        nurbs.set_degree(3);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim, int parameteried_type>
    ENUM_NURBS global_3degree_curve_interpolate_with_ends_tangent(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Vector<T, dim> &D0, 
        const Eigen::Vector<T, dim> &D1, nurbs_curve<T, dim, false, -1, -1> &nurbs)
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

        return global_3degree_curve_interpolate_with_ends_tangent<T, dim>(points, D0, D1, params, nurbs);
    }




}




















