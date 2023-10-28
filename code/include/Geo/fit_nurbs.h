#pragma once

#include <concepts>
#include "nurbs_curve.h"
#include "nurbs_surface.h"
#include "create_nurbs_arc.h"
#include <Eigen/Sparse>
#include <array>
#include <set>

namespace tnurbs
{

    template<typename T, int dim>
    ENUM_NURBS make_params_by_chord(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, std::vector<T> &params)
    {
        params.clear();
        int points_count = points.cols();
        params.resize(points_count);

        T d = 0.0;
        std::vector<T> vec_norm(points_count - 1);
        for (int index = 1; index < points_count; ++index)
        {
            vec_norm[index - 1] = (points.col(index) - points.col(index - 1)).norm();
            d += vec_norm[index - 1];
        }
        if (d < DEFAULT_ERROR)
        {
            return ENUM_NURBS::NURBS_CHORD_IS_ZERO;
        }
        
        params[0] = 0.0;
        //可以优化以下参数的取值范围?
        params[points_count - 1] = 1.0;
        for (int index = 1; index < points_count - 1; ++index)
        {
            params[index] = params[index - 1] + vec_norm[index - 1] / d;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename T, int dim>
    ENUM_NURBS make_params_by_center(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, std::vector<T> &params)
    {
        params.clear();
        int points_count = points.cols();
        params.resize(points_count);
        T d = 0.0;
        std::vector<T> vec_norm(points - 1);
        for (int index = 1; index < points_count; ++index)
        {
            vec_norm[index - 1] = std::sqrt((points.col(index) - points.col(index - 1)).norm());
            d += vec_norm[index - 1];
        }
        if (d < DEFAULT_ERROR)
        {
            return ENUM_NURBS::NURBS_CHORD_IS_ZERO;
        }
        params[0] = 0.0;
        //可以优化以下参数的取值范围?
        params[points_count - 1] = 1.0;
        for (int index = 1; index < points_count - 1; ++index)
            params[index] = params[index - 1] + vec_norm[index - 1] / d;
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    enum ENPARAMETERIEDTYPE
    {
        CHORD = 0,
        CENTRIPETAL
    };


    template<typename T, int dim, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS make_surface_params(const Eigen::Vector<Eigen::Matrix<T, dim, Eigen::Dynamic>, Eigen::Dynamic> &points, 
        std::vector<T> &u_params, std::vector<T> &v_params)
    {
        int u_valid_num = points[0].cols();
        int u_count = u_valid_num;
        int v_valid_num = points.rows();
        int v_count = v_valid_num;
        
        u_params.resize(u_count, 0.0);
        v_params.resize(v_count, 0.0);

        for (int v_index = 0; v_index < v_count; ++v_index)
        {
            std::vector<T> u_temp_params;
            if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
            {
                ENUM_NURBS flag = make_params_by_chord(points[v_index], u_temp_params);
                if (flag == ENUM_NURBS::NURBS_CHORD_IS_ZERO)
                    v_valid_num -= 1;
                else
                {
                    if (u_temp_params.size() != static_cast<size_t> (u_count))
                        return ENUM_NURBS::NURBS_ERROR;
                }
            }
            else if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
            {
                ENUM_NURBS flag = make_params_by_center(points[v_index], u_temp_params);
                if (flag == ENUM_NURBS::NURBS_CHORD_IS_ZERO)
                    v_valid_num -= 1;
                else
                {
                    if (u_temp_params.size() != static_cast<size_t> (u_count))
                        return ENUM_NURBS::NURBS_ERROR;
                }
            }
            else
            {
                return ENUM_NURBS::NURBS_ERROR;
            }
            
            for (int index = 0; index < u_count; ++index)
                u_params[index] += u_temp_params[index];
        }
        if (v_valid_num == 0)
            return ENUM_NURBS::NURBS_ERROR;
        for (int index = 0; index < u_count; ++index)
            u_params[index] /= static_cast<T>(v_valid_num);

        //reverse matrix
        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> new_control_points;
        new_control_points.resize(u_count);
        for (int col_index = 0; col_index < u_count; ++col_index)
        {
            new_control_points[col_index].resize(dim, v_count);
            for (int row_index = 0; row_index < v_count; ++row_index)
            {
                new_control_points[col_index].col(row_index) = points[row_index].col(col_index);
            }
        }

        for (int u_index = 0; u_index < u_count; ++u_index)
        {
            std::vector<T> v_temp_params;
            if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
            {
                ENUM_NURBS flag = make_params_by_chord(new_control_points[u_index], v_temp_params);
                if (flag == ENUM_NURBS::NURBS_CHORD_IS_ZERO)
                    u_valid_num -= 1;
                else
                {
                    if (v_temp_params.size() != static_cast<size_t> (v_count))
                        return ENUM_NURBS::NURBS_ERROR;
                }
            }
            else if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
            {
                ENUM_NURBS flag = make_params_by_center(new_control_points[u_index], v_temp_params);
                if (flag == ENUM_NURBS::NURBS_CHORD_IS_ZERO)
                    u_valid_num -= 1;
                else
                {
                    if (v_temp_params.size() != static_cast<size_t> (v_count))
                        return ENUM_NURBS::NURBS_ERROR;
                }
            }
            else
            {
                return ENUM_NURBS::NURBS_ERROR;
            }
            for (int index = 0; index < v_count; ++index)
                v_params[index] += v_temp_params[index];
        }
        if (v_valid_num == 0)
            return ENUM_NURBS::NURBS_ERROR;
        for (int index = 0; index < v_count; ++index)
            v_params[index] /= static_cast<T>(u_valid_num);

        return ENUM_NURBS::NURBS_SUCCESS;
    }


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

    template<typename T, int dim, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS global_curve_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, int degree, 
            nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, params);
                
        }
        else if constexpr(parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
        {
            make_params_by_center<T, dim>(points, params);
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
        Eigen::Vector3d ve;
        ve.transpose();
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

    template<typename T, int dim, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS global_curve_interpolate_with_ends_tangent(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Vector<T, dim> &D0, 
        const Eigen::Vector<T, dim> &D1, int degree,  nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, params);
                
        }
        else if constexpr(parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
        {
            make_params_by_center<T, dim>(points, params);
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

    template<typename T, int dim, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS global_3degree_curve_interpolate_with_ends_tangent(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Vector<T, dim> &D0, 
        const Eigen::Vector<T, dim> &D1, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, params);           
        }
        else if constexpr(parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
        {
            make_params_by_center<T, dim>(points, params);
        }
        else
        {
            return ENUM_NURBS::NURBS_ERROR;
        }

        return global_3degree_curve_interpolate_with_ends_tangent<T, dim>(points, D0, D1, params, nurbs);
    }

    template<typename T, int dim, int degree>
    ENUM_NURBS global_2or3degree_hermite_curve(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Matrix<T, dim, Eigen::Dynamic> &ders, 
        const std::vector<T> &params, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        static_assert((degree == 2 || degree == 3), "the degree of interpolater curve is not 2 or three");
        int points_count = points.cols();
        int params_count = params.size();
        int ders_count = ders.cols();
        if (points_count != params_count || params_count != ders_count)
        {
            return ENUM_NURBS::NURBS_ERROR;
        }

        int knots_vector_size = 2 * points_count + degree + 1;
        Eigen::VectorX<T> knots_vector(knots_vector_size);

        if constexpr (degree == 2)
        {
            knots_vector.template block<2, 1>(0, 0).setConstant(0.0);
            knots_vector.template block<3, 1>(knots_vector_size - 3, 0).setConstant(1.0);
            int count = points_count - 1;
            for (int index = 0; index < count; ++index)
            {
                knots_vector[2 + index * 2] = params[index];
                knots_vector[3 + index * 2] = (params[index] + params[index + 1]) / 2.0;
            }
        }
        else if constexpr (degree == 3)
        {
            knots_vector.template block<4, 1>(0, 0).setConstant(0.0);
            knots_vector[4] = params[1] / 2.0;
            knots_vector.template block<4, 1>(knots_vector_size - 4, 0).setConstant(1.0);
            int count = points_count - 2;
            knots_vector[knots_vector_size - 5] = (params[count] + 1.0) / 2.0;
            
            for (int index = 1; index < count; ++index)
            {
                knots_vector[3 + index * 2] = (2 * params[index] + params[index + 1]) / 3.0;
                knots_vector[4 + index * 2] = (params[index] + 2 * params[index + 1]) / 3.0;
            }
        }

        int matrix_size = 2 * points_count;
        Eigen::SparseMatrix<T> mat(matrix_size, matrix_size);
        mat.reserve(Eigen::VectorXi::Constant(matrix_size, degree + 1));

        Eigen::Matrix<T, Eigen::Dynamic, dim> Q_and_D(matrix_size, dim);
        for (int row = 0; row < points_count; ++row)
        {
            int index = -1;
            find_span<T, degree>(params[row], knots_vector, index);
            Eigen::Matrix<T, degree + 1, 2> point_and_tangent;
            ders_basis_funs<T, degree, 1>(index, params[row], knots_vector, point_and_tangent);

            for (int i = 0; i <= degree; ++i)
            {
                mat.insert(row * 2, index + i - degree) = point_and_tangent(i, 0);
                mat.insert(row * 2 + 1, index + i - degree) = point_and_tangent(i, 1);
            }
            Q_and_D.row(row * 2) = points.col(row).transpose();
            Q_and_D.row(row * 2 + 1) = ders.col(row).transpose();
        }
        mat.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
        solver.compute(mat);
        if (solver.info() != Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        Eigen::Matrix<T, Eigen::Dynamic, dim> control_points = solver.solve(Q_and_D);
        if (solver.info() != Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        Eigen::Matrix<T, dim, Eigen::Dynamic> transpose_points = control_points.transpose();

        nurbs.set_control_points(transpose_points);
        nurbs.set_knots_vector(knots_vector);
        nurbs.set_degree(degree);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim, int degree, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS global_2or3degree_hermite_curve(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Matrix<T, dim, Eigen::Dynamic> &ders,
        nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        static_assert((degree == 2 || degree == 3), "the degree of interpolater curve is not 2 or three");
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, params);           
        }
        else if constexpr(parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
        {
            make_params_by_center<T, dim>(points, params);
        }
        else
        {
            return ENUM_NURBS::NURBS_ERROR;
        }

        return global_2or3degree_hermite_curve<T, dim, degree>(points, ders, params, nurbs);
    }


    template<typename T, int dim>
    ENUM_NURBS global_surface_interpolate(const Eigen::Vector<Eigen::Matrix<T, dim, Eigen::Dynamic>, Eigen::Dynamic> &points, int u_degree, int v_degree, 
        const std::vector<T> &u_params, const std::vector<T> &v_params, nurbs_surface<T, dim, -1, -1, -1, -1, false> &nurbs)
    {
        int u_params_size = static_cast<int> (u_params.size());
        int v_params_size = static_cast<int> (v_params.size());


        Eigen::VectorX<T> u_knots_vector(u_params_size + u_degree + 1);
        Eigen::VectorX<T> v_knots_vector(v_params_size + v_degree + 1);
        u_knots_vector.block(0, 0, u_degree + 1, 1).setConstant(0.0);
        v_knots_vector.block(0, 0, v_degree + 1, 1).setConstant(0.0);
        u_knots_vector.block(u_params_size, 0, u_degree + 1, 1).setConstant(1.0);
        v_knots_vector.block(v_params_size, 0, v_degree + 1, 1).setConstant(1.0);
        int u_internel_knots_count = u_params_size - u_degree - 1;
        int v_internel_knots_count = v_params_size - v_degree - 1;
        T u_inverse_degree = 1.0 / static_cast<T>(u_degree);
        T v_inverse_degree = 1.0 / static_cast<T>(v_degree);
        
        auto u_it = u_params.begin() + 1;
        for (int j = 1; j <= u_internel_knots_count; ++j)
        {
            //求和可以优化
            u_knots_vector[u_degree + j] = std::accumulate(u_it, u_it + u_degree, 0.0) * u_inverse_degree;
            ++u_it;
        }

        auto v_it = v_params.begin() + 1;
        for (int j = 1; j <= v_internel_knots_count; ++j)
        {
            //求和可以优化
            v_knots_vector[v_degree + j] = std::accumulate(v_it, v_it + v_degree, 0.0) * v_inverse_degree;
            ++v_it;
        }

        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> R(u_params_size);
        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> new_control_points(v_params_size);
        for (int index = 0; index < u_params_size; ++index)
        {
            R[index].resize(dim, v_params_size);
        }

        for (int index = 0; index < v_params_size; ++index)
        {
            new_control_points[index].resize(dim, u_params_size);
        }

        for (int v_index = 0; v_index < v_params_size; ++v_index)
        {
            nurbs_curve<T, dim, false, -1, -1> temp_nurbs;
            global_curve_interpolate<T, dim>(points[v_index], u_degree, u_params, temp_nurbs);
            Eigen::Matrix<T, dim, Eigen::Dynamic> control_points = temp_nurbs.get_control_points();
            for (int index = 0; index < u_params_size; ++index)
            {
                R[index].col(v_index) = control_points.col(index);
            }
        }
        for (int u_index = 0; u_index < u_params_size; ++u_index)
        {
            nurbs_curve<T, dim, false, -1, -1> temp_nurbs;
            global_curve_interpolate<T, dim>(R[u_index], v_degree, v_params, temp_nurbs);
            Eigen::Matrix<T, dim, Eigen::Dynamic> control_points = temp_nurbs.get_control_points();
            for (int index = 0; index < v_params_size; ++index)
            {
                new_control_points[index].col(u_index) = control_points.col(index);
            }
        }

        nurbs.set_control_points(new_control_points);
        nurbs.set_uv_knots(u_knots_vector, v_knots_vector);
        nurbs.set_uv_degree(u_degree, v_degree);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS global_surface_interpolate(const Eigen::Vector<Eigen::Matrix<T, dim, Eigen::Dynamic>, Eigen::Dynamic> &points, int u_degree, int v_degree, 
            nurbs_surface<T, dim, -1, -1, -1, -1, false> &nurbs)
    {
        std::vector<T> u_params, v_params;
        make_surface_params<T, dim, parameteried_type>(points, u_params, v_params);
        return global_surface_interpolate<T, dim>(points, u_degree, v_degree, u_params, v_params, nurbs);
    }


    //tangent里面的向量都三单位向量
    //flag = true表示保留角点; flag = false表示不保留角点
    template<typename T, int dim, bool flag>
    ENUM_NURBS make_tangent_by_5points(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, Eigen::Matrix<T, dim, Eigen::Dynamic> &tangent)
    {
        int points_count = points.cols();
        int coeff_count = points_count + 3;
        Eigen::Matrix<T, dim, Eigen::Dynamic> Q(dim, coeff_count);
        Q.block(0, 2, dim, points_count - 1) = points.block(0, 1, dim, points_count - 1) - points.block(0, 0, dim, points_count - 1);
        Q.col(1) = 2.0 * Q.col(2) - Q.col(3);
        Q.col(0) = 2.0 * Q.col(1) - Q.col(2);

        tangent.resize(dim, points_count);

        T coeff1 = (Q.col(0).cross(Q.col(1))).norm();
        for (int index = 1; index <= points_count; ++index)
        {
            T alpha;
            T coeff2 = Q.col(index + 1).cross(Q.col(index + 2)).squaredNorm();
            if (coeff1 < DEFAULT_ERROR && coeff2 < DEFAULT_ERROR)
            {
                if constexpr (flag == true)
                {
                    alpha = 1.0;
                }
                else
                {
                    alpha = 0.5;
                }
            }
            else
            {
                alpha = coeff1 / (coeff1 + coeff2);
            }

            tangent.col(index - 1) = (1.0 - alpha) * Q.col(index - 1) + alpha * Q.col(index);
            coeff1 = coeff2;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    //false 表示无理插值, true表示有理插值; 此函数仅在二次nurbs插值中使用
    template<typename T, int dim, bool is_rational = false>
    ENUM_NURBS evaluate_intersct_point(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Matrix<T, dim, Eigen::Dynamic> &tangents, 
        std::vector<Eigen::Vector<T, dim>> &new_interpolate_points, Eigen::Matrix<T, dim, Eigen::Dynamic> &control_points)
    {
        int points_count = points.cols();
        int tangents_count = tangents.cols();
        if (points_count != tangents_count)
        {
            return ENUM_NURBS::NURBS_ERROR;
        }

        //一次分配足够的空间(最差会浪费一半的内存)
        std::vector<Eigen::Vector<T, dim>> temp_control_points;
        temp_control_points.reserve(2 * points_count);
        new_interpolate_points.reserve(2 * points_count);
        temp_control_points.push_back(points.col(0));
        new_interpolate_points.push_back(points.col(0));
        for (int index = 0; index <= points_count - 2; ++index)
        {
            Eigen::Vector<T, 2> intersect_params;
            ENUM_NURBS flag = tow_line_intersect<T, dim>(points.col(index), tangents.col(index), points.col(index + 1), tangents.col(index + 1), intersect_params);
            if (flag == ENUM_NURBS::NURBS_SUCCESS)
            {
                if (intersect_params[0] > 0.0 && intersect_params[1] < 0.0)
                {
                    Eigen::Vector<T, dim> intersetc_point = (points.col(index) + intersect_params[0] * tangents.col(index) + points.col(index + 1) + \
                                                            intersect_params[1] * tangents.col(index + 1)) / 2.0;
                    temp_control_points.push_back(std::move(intersetc_point));
                    new_interpolate_points.push_back(points.col(index + 1));
                }
                else
                {
                    Eigen::Vector<T, dim> t1 = tangents.col(index);
                    Eigen::Vector<T, dim> t2 = tangents.col(index + 1);
                    Eigen::Vector<T, dim> q = points.col(index) - points.col(index + 1);
                    T len = q.norm();
                    t1.normalize();
                    t2.normalize();
                    q.normalize();
                    T costheta1 = std::abs(t1.dot(q));
                    T costheta2 = std::abs(t2.dot(q));
                    T num;
                    if constexpr (is_rational == false)
                        num = 0.75 * len;
                    else
                        num = 1.5 * len;
                    T gamma1 = num / (2 * costheta2 + costheta1);
                    T gamma2 = num / (2 * costheta1 + costheta2);

                    Eigen::Vector<T, dim> R1 = points.col(index) + gamma1 * t1;
                    Eigen::Vector<T, dim> R2 = points.col(index + 1) - gamma2 * t2;
                    Eigen::Vector<T, dim> Q = (gamma2 * R1 + gamma1 * R2) / (gamma1 + gamma2);
                    temp_control_points.push_back(std::move(R1));
                    temp_control_points.push_back(std::move(R2));

                    new_interpolate_points.push_back(std::move(Q));
                    new_interpolate_points.push_back(points.col(index + 1));
                }
            }
            else
            {
                Eigen::Vector<T, dim> t1 = tangents.col(index);
                Eigen::Vector<T, dim> t2 = tangents.col(index + 1);
                Eigen::Vector<T, dim> q = points.col(index) - points.col(index + 1);
                t1.normalize();
                t2.normalize();
                q.normalize();
                if (t1.cross(q).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR && t2.cross(q).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR)
                {
                    Eigen::Vector<T, dim> R = (points.col(index) + points.col(index + 1)) / 2.0;
                    temp_control_points.push_back(std::move(R));
                    new_interpolate_points.push_back(points.col(index + 1));
                }
                else if (t1.cross(t2).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR)
                {
                    T gamma = ((points.col(index) - points.col(index + 1)) / 2.0).norm();
                    if (gamma < DEFAULT_ERROR)
                        return ENUM_NURBS::NURBS_ERROR;

                    Eigen::Vector<T, dim> R1 = points.col(index) + gamma * t1;
                    Eigen::Vector<T, dim> R2 = points.col(index + 1) - gamma * t2;
                    Eigen::Vector<T, dim> Q = (gamma * R1 + gamma * R2) / (gamma + gamma);
                    temp_control_points.push_back(std::move(R1));
                    temp_control_points.push_back(std::move(R2));
                    
                    new_interpolate_points.push_back(std::move(Q));
                    new_interpolate_points.push_back(points.col(index + 1));
                }
                else
                {
                    return ENUM_NURBS::NURBS_ERROR;
                }
            }
        }
        temp_control_points.push_back(points.col(points_count - 1));
        int control_points_count = temp_control_points.size();
        control_points.resize(dim, control_points_count);
        for (int index = 0; index < control_points_count; ++index)
        {
            control_points.col(index) = temp_control_points[index];
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    //此函数仅在二次nurbs插值中使用
    template<typename T, int dim>
    ENUM_NURBS local_2degree_make_params(const std::vector<Eigen::Vector<T, dim>> &new_interpolate_points, const Eigen::Matrix<T, dim, Eigen::Dynamic> &control_points, 
        const Eigen::VectorX<T> &weight, Eigen::VectorX<T> &knots_vector)
    {
        int points_count = new_interpolate_points.size();
        int control_points_count = control_points.cols();
        if (points_count != control_points_count - 1)
            return ENUM_NURBS::NURBS_ERROR;
        knots_vector.resize(points_count + 4);
        knots_vector.template block<3, 1>(0, 0).setConstant(0.0);
        knots_vector.template block<3, 1>(points_count + 1, 0).setConstant(1.0);
        knots_vector[3] = 1.0;
        for (int index = 4; index <= points_count; ++index)
        {
            T len1 = (control_points.col(index - 2) - new_interpolate_points[index - 3]).norm();
            T len2 = (new_interpolate_points[index - 3] - control_points.col(index - 3)).norm();
            knots_vector[index] = knots_vector[index - 1] + (knots_vector[index - 1] - knots_vector[index - 2]) * len1 * weight[index - 2] / (len2 * weight[index - 3]) ;
        }

        T len1 = (control_points.col(points_count - 1) - new_interpolate_points[points_count - 2]).norm();
        T len2 = (new_interpolate_points[points_count - 2] - control_points.col(points_count - 2)).norm();
        T un = knots_vector[points_count] + (knots_vector[points_count] - knots_vector[points_count - 1]) * len1 * weight[points_count - 1] / (len2 * weight[points_count - 2]);
        knots_vector.block(3, 0, points_count - 2, 1) /= un;
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    //此函数仅在二次nurbs插值中使用
    template<typename T, int dim>
    ENUM_NURBS local_2degree_make_params(const std::vector<Eigen::Vector<T, dim>> &new_interpolate_points, const Eigen::Matrix<T, dim, Eigen::Dynamic> &control_points, Eigen::VectorX<T> &knots_vector)
    {
        int points_count = new_interpolate_points.size();
        int control_points_count = control_points.cols();
        if (points_count != control_points_count - 1)
            return ENUM_NURBS::NURBS_ERROR;
        knots_vector.resize(points_count + 4);
        knots_vector.template block<3, 1>(0, 0).setConstant(0.0);
        knots_vector.template block<3, 1>(points_count + 1, 0).setConstant(1.0);
        knots_vector[3] = 1.0;
        for (int index = 4; index <= points_count; ++index)
        {
            T len1 = (control_points.col(index - 2) - new_interpolate_points[index - 3]).norm();
            T len2 = (new_interpolate_points[index - 3] - control_points.col(index - 3)).norm();
            knots_vector[index] = knots_vector[index - 1] + (knots_vector[index - 1] - knots_vector[index - 2]) * len1 / len2;
        }

        T len1 = (control_points.col(points_count - 1) - new_interpolate_points[points_count - 2]).norm();
        T len2 = (new_interpolate_points[points_count - 2] - control_points.col(points_count - 2)).norm();
        T un = knots_vector[points_count] + (knots_vector[points_count] - knots_vector[points_count - 1]) * len1 / len2;
        knots_vector.block(3, 0, points_count - 2, 1) /= un;
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 平面上的点的非有理, 二次曲线的局部插值;推荐dim = 2时使用, 当dim > 2时, 需要保证所有的点全部在某个平面上
    /// @tparam T double float int...
    /// @tparam dim 点所在的欧式空间的维数
    /// @tparam flag 是否保留角点
    /// @param points 插值点
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim>
    ENUM_NURBS local_2degree_parabola_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Matrix<T, dim, Eigen::Dynamic> &tangents,
             nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        Eigen::Matrix<T, dim, Eigen::Dynamic> control_points;
        std::vector<Eigen::Vector<T, dim>> new_interpolate_points;
        ENUM_NURBS flag = evaluate_intersct_point<T, dim>(points, tangents, new_interpolate_points, control_points);

        if (flag != ENUM_NURBS::NURBS_SUCCESS)
            return flag;

         //计算节点
        Eigen::VectorX<T> knots_vector;
        local_2degree_make_params<T, dim>(new_interpolate_points, control_points, knots_vector);
        nurbs.set_control_points(control_points);
        nurbs.set_knots_vector(knots_vector);
        nurbs.set_degree(2);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 平面上的点的非有理, 二次曲线的局部插值;推荐dim = 2时使用, 当dim > 2时, 需要保证所有的点全部在某个平面上
    /// @tparam T double float int...
    /// @tparam dim 点所在的欧式空间的维数
    /// @tparam flag 是否保留角点
    /// @param points 插值点
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim, bool flag>
    ENUM_NURBS local_2degree_parabola_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        Eigen::Matrix<T, dim, Eigen::Dynamic> tangents;
        make_tangent_by_5points<T, dim, flag>(points, tangents);
        return local_2degree_parabola_interpolate<T, dim>(points, tangents, nurbs);
    }

    /// @brief 平面上的点的圆弧曲线的局部插值;推荐dim = 2时使用, 当dim > 2时, 需要保证所有的点全部在某个平面上
    /// @tparam T double float int...
    /// @tparam dim 点所在的欧式空间的维数
    /// @tparam flag 是否保留角点
    /// @param points 插值点
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim>
    ENUM_NURBS local_2degree_arc_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Matrix<T, dim, Eigen::Dynamic> &tangents,
             nurbs_curve<T, dim, true, -1, -1> &nurbs)
    {
        Eigen::Matrix<T, dim, Eigen::Dynamic> control_points;
        std::vector<Eigen::Vector<T, dim>> new_interpolate_points;
        ENUM_NURBS flag = evaluate_intersct_point<T, dim, true>(points, tangents, new_interpolate_points, control_points);

        if (flag != ENUM_NURBS::NURBS_SUCCESS)
            return flag;

        //计算权重
        int new_interpolate_points_count = new_interpolate_points.size();
        Eigen::Vector<T, Eigen::Dynamic> weight(new_interpolate_points_count + 1);
        weight[0] = 1.0;
        weight[new_interpolate_points_count] = 1.0;
        for (int index = 1; index < new_interpolate_points_count; ++index)
        {
            Eigen::Vector<T, dim> v1 = control_points.col(index) - new_interpolate_points[index - 1];
            Eigen::Vector<T, dim> v2 = control_points.col(index) - new_interpolate_points[index];
            T len1 = v1.squaredNorm();
            T len2 = v2.squaredNorm();
            if (v1.cross(v2).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR)
            {
                weight[index] = 1.0;
            }
            else if (std::abs(len1 - len2) < DEFAULT_ERROR * DEFAULT_ERROR)
            {
                Eigen::Vector<T, dim> v3 = new_interpolate_points[index] - new_interpolate_points[index - 1];
                v1.normalize();
                v3.normalize();
                T w = v1.dot(v3);
                if (w < 0.0)
                    return ENUM_NURBS::NURBS_ERROR;
                weight[index] = w;
            }
            else
            {
                Eigen::Vector<T, dim> v3 = new_interpolate_points[index] - new_interpolate_points[index - 1];
                Eigen::Vector<T, dim> M = 0.5 * (new_interpolate_points[index] + new_interpolate_points[index - 1]);
                Eigen::Vector<T, dim> mr = control_points.col(index) - M;
                v1.normalize();
                v2.normalize();
                v3.normalize();
                Eigen::Vector<T, dim> dir1 = (v1 + v3) / 2.0;
                Eigen::Vector<T, dim> dir2 = (v2 - v3) / 2.0;
                Eigen::Vector<T, 2> interset_params;
                ENUM_NURBS flag = tow_line_intersect<T, dim>(new_interpolate_points[index - 1], dir1, M, mr, interset_params);
                if (flag != ENUM_NURBS::NURBS_SUCCESS)
                {
                    return flag;
                }
                Eigen::Vector<T, dim> s1 = (interset_params[0] * dir1 + new_interpolate_points[index - 1] + interset_params[1] * mr + M) / 2.0;

                flag = tow_line_intersect<T, dim>(new_interpolate_points[index], dir2, M, mr, interset_params);
                if (flag != ENUM_NURBS::NURBS_SUCCESS)
                {
                    return flag;
                }
                Eigen::Vector<T, dim> s2 = (interset_params[0] * dir2 + new_interpolate_points[index] + interset_params[1] * mr + M) / 2.0;
                Eigen::Vector<T, dim> s = (s1 + s2) / 2.0;
                T sw = (s - M).norm() / mr.norm();
                if (sw < DEFAULT_ERROR || (1.0 - sw) < DEFAULT_ERROR)
                    return ENUM_NURBS::NURBS_ERROR;
                weight[index] = sw / (1.0 - sw);
            }
        }
        weight[new_interpolate_points_count] = 1.0;

         //计算节点
        Eigen::VectorX<T> knots_vector;
        local_2degree_make_params<T, dim>(new_interpolate_points, control_points, weight, knots_vector);

        nurbs.set_control_points(control_points, weight);
        nurbs.set_knots_vector(knots_vector);
        nurbs.set_degree(2);
        Eigen::Matrix<T, dim + 1, Eigen::Dynamic> tempppppp = nurbs.get_control_points();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 平面上的点的圆弧曲线的局部插值;推荐dim = 2时使用, 当dim > 2时, 需要保证所有的点全部在某个平面上
    /// @tparam T double float int...
    /// @tparam dim 点所在的欧式空间的维数
    /// @tparam flag 是否保留角点
    /// @param points 插值点
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim, bool flag = false>
    ENUM_NURBS local_2degree_arc_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, nurbs_curve<T, dim, true, -1, -1> &nurbs)
    {
        Eigen::Matrix<T, dim, Eigen::Dynamic> tangents;
        make_tangent_by_5points<T, dim, flag>(points, tangents);
        return local_2degree_arc_interpolate<T, dim>(points, tangents, nurbs);
    }

    /// @brief 局部三次杨条插值
    /// @tparam T double float int ...
    /// @tparam dim 点所在的欧式空间的维数
    /// @param points 插值点
    /// @param tangents_dir 插值点的切向(切方向, 不是切向量)
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim>
    ENUM_NURBS local_3degree_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Matrix<T, dim, Eigen::Dynamic> &tangents_dir, 
        nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        int points_count = points.cols();
        int tangents_count = points.cols();
        if (points_count != tangents_count)
            return ENUM_NURBS::NURBS_ERROR;
        
        //计算控制点
        int control_points_count = 2 * points_count;
        Eigen::Matrix<T, dim, Eigen::Dynamic> control_points(dim, control_points_count);
        control_points.col(0) = points.col(0);
        for (int index = 1; index < points_count; ++index)
        {
            //TODO: 判断单位化有没有成功
            Eigen::Vector<T, dim> T0 = tangents_dir.col(index - 1).normalized();
            Eigen::Vector<T, dim> T3 = tangents_dir.col(index).normalized();
            Eigen::Vector<T, dim> T03 = T0 + T3;

            Eigen::Vector<T, dim> P30 = points.col(index) - points.col(index - 1);
            T a = 16.0 - T03.squaredNorm();
            T b = 12.0 * T03.dot(P30);
            T c = -36.0 * P30.squaredNorm();

            T delta = b * b - 4.0 * a * c;
            if (delta < DEFAULT_ERROR * DEFAULT_ERROR)
                return ENUM_NURBS::NURBS_ERROR;
            T sqr_delta = std::sqrt(delta);
            T x1 = ((-1.0) * b - sqr_delta) / (2.0 * a);
            T x2 = ((-1.0) * b + sqr_delta) / (2.0 * a);
            T alpha = x1 > 0.0 ? x1 : x2;
            control_points.col(2 * index - 1) = points.col(index - 1) + (alpha / 3.0) * T0;
            control_points.col(2 * index) = points.col(index) - (alpha / 3.0) * T3;
        }
        control_points.col(control_points_count - 1) = points.col(points_count - 1);

        //计算节点矢量
        int konts_count = control_points_count + 4;
        Eigen::VectorX<T> knots_vector(konts_count);
        knots_vector.template block<6, 1>(0, 0).setConstant(0.0);
        knots_vector.template block<4, 1>(konts_count - 4, 0).setConstant(1.0);

        T u = 0.0;
        for (int index = 1; index <= points_count - 2; ++index)
        {
            u += 3.0 * (control_points.col(2 * index - 1) - points.col(index - 1)).norm();
            knots_vector[2 + index * 2] = knots_vector[3 + index * 2] = u;
        }
        u += 3.0 * (control_points.col(2 * (points_count - 1) - 1) - points.col(points_count - 2)).norm();
        knots_vector.block(4, 0, konts_count - 8, 1) /= u;
        
        nurbs.set_control_points(control_points);
        nurbs.set_knots_vector(knots_vector);
        nurbs.set_degree(3);
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 局部三次杨条插值
    /// @tparam T double float int ...
    /// @tparam dim 点所在的欧式空间的维数
    /// @param points 插值点
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim, bool flag = false>
    ENUM_NURBS local_3degree_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        Eigen::Matrix<T, dim, Eigen::Dynamic> tangents;
        make_tangent_by_5points<T, dim, flag>(points, tangents);
        return local_3degree_interpolate<T, dim>(points, tangents, nurbs);
    }


    /// @brief 局部三次杨条插值
    /// @tparam T double float int ...
    /// @tparam dim 点所在的欧式空间的维数
    /// @param points 插值点
    /// @param u_tangents_dir 插值点的u切向(切方向, 不是切向量)
    /// @param v_tangents_dir 插值点的v切向(切方向, 不是切向量)
    /// @param u_params 插值点的参数
    /// @param v_params 插值点的参数
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim>
    ENUM_NURBS local_bi3degree_interpolate(const Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> &points, const Eigen::MatrixX<Eigen::Vector<T, dim>> &u_tangents_dir,
        const Eigen::MatrixX<Eigen::Vector<T, dim>> &v_tangents_dir, const std::vector<T> &u_params, const std::vector<T> &v_params, nurbs_surface<T, dim, -1, -1, -1, -1, false> &nurbs)
    {
        int v_points_count = points.rows();
        int u_points_count = points[0].cols();

        int u_tangents_u_count = u_tangents_dir.cols();
        int u_tangents_v_count = u_tangents_dir.rows();
        int v_tangents_u_count = v_tangents_dir.cols();
        int v_tangents_v_count = v_tangents_dir.rows();
        int u_params_count = u_params.size();
        int v_params_count = v_params.size();

        if (v_points_count != u_tangents_v_count || v_points_count != v_tangents_v_count || v_points_count != v_params_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        if (u_points_count != u_tangents_u_count || u_points_count != v_tangents_u_count || u_points_count != u_params_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        
        //计算u向和v向的总弦长(下面的代码可以优化, 一次循环即可)
        std::vector<T> u_chord_lens(v_points_count);
        std::vector<T> v_chord_lens(u_points_count);
        for (int v_index = 0; v_index < v_points_count; ++v_index)
        {
            T len = 0.0;
            for (int u_index = 1; u_index < u_points_count; ++u_index)
            {
                len += (points[v_index].col(u_index) - points[v_index].col(u_index - 1)).norm();
            }
            u_chord_lens[v_index] = len;
        }

        for (int u_index = 0; u_index < u_points_count; ++u_index)
        {
            T len = 0.0;
            for (int v_index = 1; v_index < v_points_count; ++v_index)
            {
                len += (points[v_index].col(u_index) - points[v_index - 1].col(u_index)).norm();
            }
            v_chord_lens[u_index] = len;
        }

        std::vector<T> alpha_k(u_points_count - 1);
        std::vector<T> beta_l(v_points_count - 1);
        for (int index = 0; index < u_points_count - 2; ++index)
        {
            alpha_k[index] = (u_params[index + 1] - u_params[index]) / (u_params[index + 2] - u_params[index]);
        }
        T temp = u_params[u_points_count - 1] - u_params[u_points_count - 2];
        alpha_k[u_points_count - 2] = temp / (temp + (u_params[1] - u_params[0]));
        for (int index = 0; index < v_points_count - 2; ++index)
        {
            beta_l[index] = (v_params[index + 1] - v_params[index]) / (v_params[index + 2] - v_params[index]);
        }
        temp = v_params[v_points_count - 1] - v_params[v_points_count - 2];
        beta_l[v_points_count - 2] = temp / (temp + (v_params[1] - v_params[0]));

        Eigen::MatrixX<Eigen::Vector<T, dim>> Duv(v_points_count, u_points_count);
        for (int v_index = 1; v_index < v_points_count; ++v_index)
        {
            for (int u_index = 1; u_index < u_points_count; ++u_index)
            {
                int last_v_index = v_index + 1 == v_points_count ? 1 : v_index + 1;
                int last_u_index = u_index + 1 == u_points_count ? 1 : u_index + 1;
                Eigen::Vector<T, dim> first_u_dir = u_chord_lens[v_index - 1] * u_tangents_dir(v_index - 1, u_index).normalized();
                Eigen::Vector<T, dim> second_u_dir = u_chord_lens[v_index] * u_tangents_dir(v_index, u_index).normalized();
                Eigen::Vector<T, dim> third_u_dir = u_chord_lens[last_v_index] * u_tangents_dir(last_v_index, u_index).normalized();

                Eigen::Vector<T, dim> first_v_dir = v_chord_lens[u_index - 1] * v_tangents_dir(u_index - 1, v_index).normalized();
                Eigen::Vector<T, dim> second_v_dir = v_chord_lens[u_index] * v_tangents_dir(u_index, v_index).normalized();
                Eigen::Vector<T, dim> third_v_dir = v_chord_lens[last_u_index] * v_tangents_dir(last_u_index, v_index).normalized();

                Eigen::Vector<T, dim> duv = ((1.0 - beta_l[v_index - 1]) / (v_params[v_index] - v_params[v_index - 1])) * (second_u_dir - first_u_dir) + \
                                            (beta_l[v_index - 1] / (v_params[v_index + 1] - v_params[v_index])) * (third_u_dir - second_u_dir);
                Eigen::Vector<T, dim> dvu = ((1.0 - alpha_k[u_index - 1]) / (u_params[u_index] - u_params[u_index - 1])) * (second_v_dir - first_v_dir) + \
                                            (alpha_k[u_index - 1] / (u_params[u_index + 1] - u_params[u_index])) * (third_v_dir - second_v_dir);
                Duv(v_index, u_index) = (alpha_k[u_index - 1] * duv + beta_l[v_index - 1] * dvu) / (alpha_k[u_index - 1] + beta_l[v_index - 1]);
            }
        }

        //将一半边界处的Duv设置为另一半边界的Duv
        Duv(0, 0) = Duv(v_points_count - 1, u_points_count - 1);
        for (int index = 1; index < u_points_count; ++index)
        {
            Duv(0, index) = Duv(v_points_count - 1, index);
        }
        for (int index = 1; index < v_points_count; ++index)
        {
            Duv(index, 0) = Duv(index, u_points_count - 1);
        }

        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> control_points(2 * v_points_count);
        control_points[0].resize(dim, 2 * u_points_count);
        control_points[2 * v_points_count - 1].resize(dim, 2 * u_points_count);
        control_points[0].col(0) = points[0].col(0);
        control_points[0].col(2 * u_points_count - 1) = points[0].col(u_points_count - 1);
        control_points[2 * v_points_count - 1].col(0) = points[v_points_count - 1].col(0);
        control_points[2 * v_points_count - 1].col(2 * u_points_count - 1) = points[v_points_count - 1].col(u_points_count - 1);

        Eigen::MatrixX<Eigen::Vector<T, dim>> u_temp_control_points(v_points_count, 2 * u_points_count - 2);
        for (int v_index = 0; v_index < v_points_count; ++v_index)
        {
            for (int u_index = 0; u_index < u_points_count - 1; ++u_index)
            {
                T a = u_chord_lens[v_index] * (u_params[u_index + 1] - u_params[u_index]) / 3.0;
                u_temp_control_points(v_index, u_index * 2) = points[v_index].col(u_index) + a * u_tangents_dir(v_index, u_index).normalized();
                u_temp_control_points(v_index, u_index * 2 + 1) = points[v_index].col(u_index + 1) - a * u_tangents_dir(v_index, u_index + 1).normalized();
            }
        }

        Eigen::MatrixX<Eigen::Vector<T, dim>> v_temp_control_points(v_points_count * 2 - 2, u_points_count);
        for (int u_index = 0; u_index < u_points_count; ++u_index)
        {
            for (int v_index = 0; v_index < v_points_count - 1; ++v_index)
            {
                T a = v_chord_lens[u_index] * (v_params[v_index + 1] - v_params[v_index]) / 3.0;
                v_temp_control_points(v_index * 2, u_index) = points[v_index].col(u_index) + a * v_tangents_dir(v_index, u_index).normalized();
                v_temp_control_points(v_index * 2 + 1, u_index) = points[v_index + 1].col(u_index) - a * v_tangents_dir(v_index + 1, u_index).normalized();
            }
        }
        //中间的控制点
        for (int v_index = 1; v_index < v_points_count; ++v_index)
        {
            control_points[v_index * 2 - 1].resize(dim, 2 * u_points_count);
            control_points[v_index * 2].resize(dim, 2 * u_points_count);
            T delta_l = v_params[v_index] - v_params[v_index - 1];
            for (int u_index = 1; u_index < u_points_count; ++u_index)
            {
                T delta_k = u_params[u_index] - u_params[u_index - 1];
                T gamma = (delta_k * delta_l) / 9.0;
                control_points[v_index * 2 - 1].col(u_index * 2 - 1) = gamma * Duv(v_index - 1, u_index - 1) - points[v_index - 1].col(u_index - 1) + \
                                u_temp_control_points(v_index - 1, 2 * (u_index - 1)) + v_temp_control_points(2 * (v_index - 1), u_index - 1);
                control_points[v_index * 2 - 1].col(u_index * 2) = (-1 * gamma) * Duv(v_index - 1, u_index) - points[v_index - 1].col(u_index) + \
                                u_temp_control_points(v_index - 1, 2 * u_index - 1) + v_temp_control_points(2 * (v_index - 1), u_index);
                control_points[v_index * 2].col(u_index * 2 - 1) = (-1 * gamma) * Duv(v_index, u_index - 1) - points[v_index].col(u_index - 1) + \
                                u_temp_control_points(v_index , 2 * (u_index - 1)) + v_temp_control_points(2 * v_index - 1, u_index - 1);
                control_points[v_index * 2].col(u_index * 2) = gamma * Duv(v_index, u_index) - points[v_index].col(u_index) + \
                                u_temp_control_points(v_index, 2 * u_index - 1) + v_temp_control_points(2 * v_index - 1, u_index);

            }
        }
        //第一行和最后一行控制点
        for (int index = 1; index < u_points_count; ++index)
        {
            control_points[0].col(index * 2 - 1) = u_temp_control_points(0, 2 * (index - 1));
            control_points[0].col(index * 2) = u_temp_control_points(0, 2 * index - 1);
            control_points[2 * v_points_count - 1].col(index * 2 - 1) = u_temp_control_points(v_points_count - 1, 2 * (index - 1));
            control_points[2 * v_points_count - 1].col(index * 2) = u_temp_control_points(v_points_count - 1, 2 * index - 1);
        }
        //第一列和最后一列控制点
        for (int index = 1; index < v_points_count; ++index)
        {
            control_points[index * 2 - 1].col(0) = v_temp_control_points(2 * (index - 1), 0);
            control_points[index * 2].col(0) = v_temp_control_points(2 * index - 1, 0);
            control_points[index * 2 - 1].col(2 * v_points_count - 1) = v_temp_control_points(2 * (index - 1), u_points_count - 1);
            control_points[index * 2].col(2 * v_points_count - 1) = v_temp_control_points(2 * index - 1, u_points_count - 1);
        }

        Eigen::VectorX<T> u_knots(4 + u_params_count * 2), v_knots(4 + v_params_count * 2);
        u_knots.template block<4, 1>(0, 0).setConstant(u_params[0]);
        u_knots.template block<4, 1>(2 * u_params_count, 0).setConstant(u_params[u_params_count - 1]);
        for (int index = 1; index <= u_params_count - 2; ++index)
        {
            u_knots[2 * index + 2] = u_params[index];
            u_knots[2 * index + 3] = u_params[index];
        }

        v_knots.template block<4, 1>(0, 0).setConstant(v_params[0]);
        v_knots.template block<4, 1>(2 * v_params_count, 0).setConstant(v_params[v_params_count - 1]);
        for (int index = 1; index <= v_params_count - 2; ++index)
        {
            v_knots[2 * index + 2] = v_params[index];
            v_knots[2 * index + 3] = v_params[index];
        }
        nurbs.set_control_points(control_points);
        nurbs.set_uv_degree(3, 3);
        nurbs.set_uv_knots(u_knots, v_knots);

        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 局部三次杨条插值
    /// @tparam T double float int ...
    /// @tparam dim 点所在的欧式空间的维数
    /// @param points 插值点
    /// @param u_params 插值点的参数
    /// @param v_params 插值点的参数
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim>
    ENUM_NURBS local_bi3degree_interpolate(const Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> &points, const std::vector<T> &u_params, const std::vector<T> &v_params, nurbs_surface<T, dim, -1, -1, -1, -1, false> &nurbs)
    {
        int u_points_count = u_params.size();
        int v_points_count = v_params.size();
        if (u_points_count != points[0].cols())
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        if (v_points_count != points.rows())
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;

        Eigen::MatrixX<Eigen::Vector<T, dim>> u_tangents_dir(v_points_count, u_points_count);
        Eigen::MatrixX<Eigen::Vector<T, dim>> v_tangents_dir(v_points_count, u_points_count);

        for (int v_index = 0; v_index < v_points_count; ++v_index)
        {
            for (int u_index = 1; u_index < u_points_count - 1; ++u_index)
            {
                T delta_k = u_params[u_index] - u_params[u_index - 1];
                T next_delta_k = u_params[u_index + 1] - u_params[u_index];
                Eigen::Vector<T, dim> dk = (points[v_index].col(u_index) - points[v_index].col(u_index - 1)) / delta_k;
                Eigen::Vector<T, dim> dk1 = (points[v_index].col(u_index + 1) - points[v_index].col(u_index)) / next_delta_k;
                T alpha = delta_k / (next_delta_k + delta_k);
                u_tangents_dir(v_index, u_index) = (1.0 - alpha) * dk + alpha * dk1;
            }
            u_tangents_dir(v_index, 0) = (2.0 / (u_params[1] - u_params[0])) * (points[v_index].col(1) - points[v_index].col(0)) - u_tangents_dir(v_index, 1);
            u_tangents_dir(v_index, u_points_count - 1) = 2.0 / (u_params[u_points_count - 1] - u_params[u_points_count - 2]) * \
                    (points[v_index].col(u_points_count - 1) - points[v_index].col(u_points_count - 2)) - u_tangents_dir(v_index, u_points_count - 2);
        }

        for (int u_index = 0; u_index < u_points_count; ++u_index)
        {
            for (int v_index = 1; v_index < v_points_count - 1; ++v_index)
            {
                T delta_k = v_params[v_index] - v_params[v_index - 1];
                T next_delta_k = v_params[v_index + 1] - v_params[v_index];
                Eigen::Vector<T, dim> dk = (points[v_index].col(u_index) - points[v_index - 1].col(u_index)) / delta_k;
                Eigen::Vector<T, dim> dk1 = (points[v_index + 1].col(u_index) - points[v_index].col(u_index)) / next_delta_k;
                T alpha = delta_k / (next_delta_k + delta_k);
                v_tangents_dir(v_index, u_index) = (1.0 - alpha) * dk + alpha * dk1;
            }
            v_tangents_dir(0, u_index) = (2.0 / (v_params[1] - v_params[0])) * (points[1].col(u_index) - points[0].col(u_index)) - v_tangents_dir(1, u_index);
            v_tangents_dir(v_points_count - 1, u_index) = (2.0 / (v_params[v_points_count - 1] - v_params[v_points_count - 2]) ) * \
                (points[v_points_count - 1].col(u_index) - points[v_points_count - 2].col(u_index)) - v_tangents_dir(v_points_count - 2, u_index);
        }

        return local_bi3degree_interpolate<T, dim>(points, u_tangents_dir, v_tangents_dir, u_params, v_params, nurbs);
    }

    /// @brief 局部三次杨条插值
    /// @tparam T double float int ...
    /// @tparam dim 点所在的欧式空间的维数
    /// @param points 插值点
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim,  ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS local_bi3degree_interpolate(const Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> &points, nurbs_surface<T, dim, -1, -1, -1, -1, false> &nurbs)
    {
        std::vector<T> u_params, v_params;
        make_surface_params<T, dim, parameteried_type>(points, u_params, v_params);
        return local_bi3degree_interpolate<T, dim>(points, u_params, v_params, nurbs);
    }

    //内部使用
    template<typename T, int dim>
    ENUM_NURBS evaulate_matrix_for_approximation(const Eigen::VectorX<T> &knots, const std::vector<T> &params, int degree, int interpolate_points_count, Eigen::SparseMatrix<T> &N)
    {
        int control_points_count = knots.size() - degree - 1;
        int Nrows = std::min(degree + 1, control_points_count - 2);

        for (int col = 0; col < interpolate_points_count - 2; ++col)
        {
            int index = -1;
            find_span<T>(params[col + 1], degree, knots, index);
            Eigen::VectorX<T> basis;
            basis_functions<T>(index, params[col + 1], degree, knots, basis);

            if (index == degree && index == control_points_count -1)
            {
                for (int i = 0; i < Nrows; ++i)
                {
                    N.insert(index + i - degree, col) = basis[i + 1];
                }
            }
            else if (index == degree)
            {
                for (int i = 0; i < Nrows; ++i)
                {
                    if (i + 1 > degree)
                        N.insert(index - degree + i, col) = 0.0;
                    else
                        N.insert(index - degree + i, col) = basis[i + 1];
                }
            }
            else if (index == control_points_count - 1)
            {
                for (int i = 0; i < Nrows; ++i)
                {
                    if (index - degree + i - 1 == control_points_count - 2)
                        N.insert(control_points_count - 2 - Nrows, col) = 0.0;
                    else
                        N.insert(index - degree + i - 1, col) = basis[i];
                }

            }
            else
            {
                for (int i = 0; i < Nrows; ++i)
                {
                    N.insert(index - degree + i - 1, col) = basis[i];
                }
            }

        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    //内部使用
    template<typename T, int dim>
    ENUM_NURBS evaulate_R_matrix_for_approximation(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::VectorX<T> &knots, const std::vector<T> &params, int degree, int interpolate_points_count, Eigen::MatrixX<T> &R)
    {
        int control_points_count = knots.size() - degree - 1;
        R.resize(interpolate_points_count - 2, dim);
        for (int col = 0; col < interpolate_points_count - 2; ++col)
        {
            int index = -1;
            find_span<T>(params[col + 1], degree, knots, index);
            Eigen::VectorX<T> basis;
            basis_functions<T>(index, params[col + 1], degree, knots, basis);
            if (index == degree && index == control_points_count -1)
            {
                R.block(col, 0, 1, dim) = (points.col(col + 1) - basis[0] * points.col(0) -  basis[degree] * points.col(interpolate_points_count - 1)).transpose();
            }
            else if (index == degree)
            {
                R.block(col, 0, 1, dim) = (points.col(col + 1) - basis[0] * points.col(0)).transpose();
            }
            else if (index == control_points_count - 1)
            {
                R.block(col, 0, 1, dim) = (points.col(col + 1) - basis[degree] * points.col(interpolate_points_count - 1)).transpose();

            }
            else
            {
                R.block(col, 0, 1, dim) = (points.col(col + 1)).transpose();
            }
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    //内部使用
    template<typename T, int dim>
    ENUM_NURBS global_least_squares_curve_approximation(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const std::vector<T> &params, 
        const Eigen::VectorX<T> &knots, int degree, int control_points_count, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        int points_count = points.cols();
        Eigen::MatrixX<T> R;
        Eigen::SparseMatrix<T> N(control_points_count - 2, points_count - 2);
        int Nrows = std::min(degree + 1, control_points_count - 2);
        N.reserve(Eigen::VectorXi::Constant(points_count - 2, Nrows));
        evaulate_matrix_for_approximation<T, dim>(knots, params, degree, points_count, N);
        evaulate_R_matrix_for_approximation<T, dim>(points, knots, params, degree, points_count, R);

        Eigen::SparseMatrix<T> mat =  N * Eigen::SparseMatrix<T, Eigen::ColMajor>(N.transpose());
        Eigen::MatrixX<T> Rs = N * R;

        mat.makeCompressed();
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver;
        solver.compute(mat);
        if (solver.info() != Eigen::Success)
        {
            return ENUM_NURBS::NURBS_ERROR;
        }
        Eigen::Matrix<T, Eigen::Dynamic, dim> control_points = solver.solve(Rs);
        if (solver.info() != Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        Eigen::Matrix<T, dim, Eigen::Dynamic> new_control_points(dim, control_points_count);
        new_control_points.col(0) = points.col(0);
        new_control_points.col(control_points_count - 1) = points.col(points_count - 1);
        new_control_points.block(0, 1, dim, control_points_count - 2) = control_points.transpose();
        nurbs.set_control_points(new_control_points);
        nurbs.set_knots_vector(knots);
        nurbs.set_degree(degree);

        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim>
    ENUM_NURBS global_least_squares_curve_approximation(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const std::vector<T> &params, int degree, 
            int control_points_count, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        int params_count = params.size();
        int points_count = points.cols();
        if (params_count != points_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        if (points_count <= control_points_count || degree >= control_points_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;

        T d = static_cast<T>(points_count) / static_cast<T>(control_points_count - degree);
        Eigen::VectorX<T> knots(degree + control_points_count + 1);
        knots.block(0, 0, degree + 1, 1).setConstant(0.0);
        knots.block(control_points_count, 0, degree + 1, 1).setConstant(1.0);
        for (int index = 1; index <= control_points_count - 1 - degree; ++index)
        {
            int i = std::floor(index * d);
            T alpha = index * d - static_cast<T>(i);
            knots[degree + index] = (1.0 - alpha) * params[i - 1] + alpha * params[i];
        }
        return global_least_squares_curve_approximation<T, dim>(points, params, knots, degree, control_points_count, nurbs);
    }


    template<typename T, int dim, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS global_least_squares_curve_approximation(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, int degree, int control_points_count, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        int points_count = points.cols();
        if (points_count <= control_points_count || degree >= control_points_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
       
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, params);
        }
        else
        {
            make_params_by_center<T, dim>(points, params);
        }
        return global_least_squares_curve_approximation<T, dim>(points, params, degree, control_points_count, nurbs);
    }


    template<typename T, int dim, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS global_wc_least_squares_curve_approximation(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const std::vector<T> &W, const Eigen::Matrix<T, dim, Eigen::Dynamic> &Ders,
        const std::vector<int> &Index, const std::vector<T> &Wd, int degree, int control_points_count, nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        int points_count = points.cols();
        int w_count = W.size();
        if (points_count != w_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        int Wd_count = Wd.size();
        int index_count = Index.size();
        int ders_count = Ders.cols();
        if (Wd_count != index_count || Wd_count != ders_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        if (points_count < ders_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;

        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, params);
        }
        else
        {
            make_params_by_center<T, dim>(points, params);
        }

        Eigen::VectorX<T> knots(control_points_count + degree + 1);
        knots.block(0, 0, degree + 1, 1).setConstant(0.0);
        knots.block(control_points_count, 0, degree + 1, 1).setConstant(1.0);
        T d = static_cast<T>(points_count) / static_cast<T>(control_points_count - degree);
        for (int index = 1; index <= control_points_count - 1 - degree; ++index)
        {
            int i = std::floor(index * d);
            T alpha = index * d - static_cast<T>(i);
            knots[degree + index] = (1.0 - alpha) * params[i - 1] + alpha * params[i];
        }

        //简单起见, 此算法不使用稀疏矩阵; 有空优化一下
        std::vector<T> positive_w;
        std::set<int> constrain_index;
        std::set<int> unconstrain_index;
        positive_w.reserve(w_count);
        for (int index = 0; index < w_count; ++index)
        {
            if(W[index] < 0.0)
            {
                constrain_index.insert(index);
            }
            else
            {
                positive_w.push_back(W[index]);
                unconstrain_index.insert(index);
            }
        }
        
        std::set<int> ders_constrain_index;
        std::set<int> ders_unconstrain_index;
        std::vector<T> positive_Wd;
        positive_Wd.reserve(Wd_count);
        for (int index = 0; index < Wd_count; ++index)
        {
            if (Wd[index] < 0.0)
            {
                ders_constrain_index.insert(index);
            }
            else
            {
                positive_Wd.push_back(Wd[index]);
                ders_unconstrain_index.insert(index);
            }
        }

        int point_constrian_count = constrain_index.size();
        int ders_constrian_count = ders_constrain_index.size();
        int point_unconstrain_count = positive_w.size();
        int ders_unconstrain_count = positive_Wd.size();

        if (point_constrian_count + ders_constrian_count >= points_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        if (point_constrian_count + ders_constrian_count + points_count > point_unconstrain_count + ders_unconstrain_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;

        Eigen::MatrixX<T> WM(point_unconstrain_count + ders_unconstrain_count, point_unconstrain_count + ders_unconstrain_count);
        WM.setConstant(0.0);
        for (int index = 0; index < point_unconstrain_count; ++index)
            WM(index, index) = positive_w[index];
        for (int index = 0; index < ders_unconstrain_count; ++index)
            WM(index + point_unconstrain_count, index + point_unconstrain_count) = positive_Wd[index];

        //下面的代码可以优化
        Eigen::MatrixX<T> MD(point_constrian_count + ders_constrian_count, control_points_count);
        Eigen::Matrix<T, Eigen::Dynamic, dim> TD(point_constrian_count + ders_constrian_count, dim);
        MD.setConstant(0.0);
        int index = 0;
        for (auto it = constrain_index.begin(); it != constrain_index.end(); ++it)
        {
            int point_index = *it;
            int span = -1;
            find_span<T>(params[point_index], degree, knots, span);
            Eigen::VectorX<T> basis;
            basis_functions<T>(span, params[point_index], degree, knots, basis);
            MD.block(index, span - degree, 1, degree + 1) = basis.transpose();
            TD.row(index) = points.col(point_index).transpose();
            ++index;
        }
        index = 0;
        for (auto it = ders_constrain_index.begin(); it != ders_constrain_index.end(); ++it)
        {
            int point_index = *it;
            int span = -1;
            T u = params[Index[point_index]];
            find_span<T>(u, degree, knots, span);
            Eigen::Matrix<T, Eigen::Dynamic, 2> ders_basis(degree + 1, 2);
            ders_basis_funs<T, 1>(span, degree,  u, knots, ders_basis);
            MD.block(index + point_constrian_count, span - degree, 1, degree + 1) = ders_basis.col(1).transpose();
            TD.row(index + point_constrian_count) = Ders.col(point_index).transpose();
            ++index;
        }

        Eigen::MatrixX<T> ND(point_unconstrain_count + ders_unconstrain_count, control_points_count);
        Eigen::Matrix<T, Eigen::Dynamic, dim> SD(point_unconstrain_count + ders_unconstrain_count, dim);
        ND.setConstant(0.0);
        index = 0;
        for (auto it = unconstrain_index.begin(); it != unconstrain_index.end(); ++it)
        {
            int point_index = *it;
            int span = -1;
            find_span<T>(params[point_index], degree, knots, span);
            Eigen::VectorX<T> basis;
            basis_functions<T>(span, params[point_index], degree, knots, basis);
            ND.block(index, span - degree, 1, degree + 1) = basis.transpose();
            SD.row(index) = points.col(point_index).transpose();
            ++index;
        }

        index = 0;
        for (auto it = ders_unconstrain_index.begin(); it != ders_unconstrain_index.end(); ++it)
        {
            int point_index = *it;
            int span = -1;
            T u = params[Index[point_index]];
            find_span<T>(u, degree, knots, span);
            Eigen::Matrix<T, Eigen::Dynamic, 2> ders_basis(degree + 1, 2);
            ders_basis_funs<T, 1>(span, degree, u, knots, ders_basis);
            ND.block(index + point_unconstrain_count, span - degree, 1, degree + 1) = ders_basis.col(1).transpose();
            SD.row(index + point_unconstrain_count) = Ders.col(point_index).transpose();
            ++index;
        }


        int mat_count = control_points_count + ders_constrian_count + point_constrian_count;
        Eigen::MatrixX<T> mat(mat_count, mat_count);
        mat.setConstant(0.0);
        mat.block(0, 0, control_points_count, control_points_count) = ND.transpose() * WM * ND;
        mat.block(control_points_count, 0, point_constrian_count + ders_constrian_count, control_points_count) = MD;
        mat.block(0, control_points_count, control_points_count, point_constrian_count + ders_constrian_count) = MD.transpose();

        Eigen::Matrix<T, Eigen::Dynamic, dim> V(point_constrian_count + ders_constrian_count + control_points_count, dim);
        V.block(0, 0, control_points_count, dim) = ND.transpose() * WM * SD;
        V.block(control_points_count, 0, point_constrian_count + ders_constrian_count, dim) = TD;

        Eigen::BDCSVD<Eigen::MatrixX<T>,  Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
        Eigen::MatrixX<T> PA = matSvd.solve(V);
        if (matSvd.info() !=  Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;

        Eigen::MatrixX<T> check_mat = mat * PA - V;
        T max_elem = check_mat.maxCoeff();
        if (max_elem > DEFAULT_ERROR)
            return ENUM_NURBS::NURBS_ERROR;
        T min_elem = check_mat.minCoeff();
        if (min_elem < -1 * DEFAULT_ERROR)
            return ENUM_NURBS::NURBS_ERROR;
        Eigen::Matrix<T, dim, Eigen::Dynamic> control_points = PA.block(0, 0, control_points_count, dim).transpose();
        nurbs.set_control_points(control_points);
        nurbs.set_knots_vector(knots);
        nurbs.set_degree(degree);
        return ENUM_NURBS::NURBS_SUCCESS;

    }

    template<typename T, int dim, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS global_surface_approximation(const Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> &points, int u_degree, int v_degree, int u_control_points_count,
        int v_control_points_count, nurbs_surface<T, dim, -1, -1, -1, -1, false> &nurbs)
    {
        int v_points_count = points.rows();
        int u_points_count = points[0].cols();
        if (u_points_count <= u_control_points_count || u_degree >= u_control_points_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        if (v_points_count <= v_control_points_count || v_degree >= v_control_points_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        std::vector<T> u_params, v_params;
        make_surface_params<T, dim, parameteried_type>(points, u_params, v_params);
        
        //填充节点矢量
        T d = static_cast<T>(u_points_count) / static_cast<T>(u_control_points_count - u_degree);
        Eigen::VectorX<T> u_knots(u_degree + u_control_points_count + 1);
        u_knots.block(0, 0, u_degree + 1, 1).setConstant(0.0);
        u_knots.block(u_control_points_count, 0, u_degree + 1, 1).setConstant(1.0);
        for (int index = 1; index <= u_control_points_count - 1 - u_degree; ++index)
        {
            int i = std::floor(index * d);
            T alpha = index * d - static_cast<T>(i);
            u_knots[u_degree + index] = (1.0 - alpha) * u_params[i - 1] + alpha * u_params[i];
        }
        d = static_cast<T>(v_points_count) / static_cast<T>(v_control_points_count - v_degree);
        Eigen::VectorX<T> v_knots(v_degree + v_control_points_count + 1);
        v_knots.block(0, 0, v_degree + 1, 1).setConstant(0.0);
        v_knots.block(v_control_points_count, 0, v_degree + 1, 1).setConstant(1.0);
        for (int index = 1; index <= v_control_points_count - 1 - v_degree; ++index)
        {
            int i = std::floor(index * d);
            T alpha = index * d - static_cast<T>(i);
            v_knots[v_degree + index] = (1.0 - alpha) * v_params[i - 1] + alpha * v_params[i];
        }
        
        std::vector<Eigen::MatrixX<T>> temp_control_points;
        temp_control_points.reserve(v_points_count);

        Eigen::SparseMatrix<T> NU(u_control_points_count - 2, u_points_count - 2);
        int Nrows = std::min(u_degree + 1, u_control_points_count - 2);
        NU.reserve(Eigen::VectorXi::Constant(u_points_count - 2, Nrows));
        evaulate_matrix_for_approximation<T, dim>(u_knots, u_params, u_degree, u_points_count, NU);
        Eigen::SparseMatrix<T> mat =  NU * Eigen::SparseMatrix<T, Eigen::ColMajor>(NU.transpose());
        mat.makeCompressed();
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver;
        solver.compute(mat);
        if (solver.info() != Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        Eigen::MatrixX<T> R;   
        for (int index = 0; index < v_points_count; ++index)
        {
            evaulate_R_matrix_for_approximation<T, dim>(points[index], u_knots, u_params, u_degree, u_points_count, R);
            Eigen::MatrixX<T> Rs = NU * R;
            Eigen::Matrix<T, Eigen::Dynamic, dim> u_temp_control_points = solver.solve(Rs);
            if (solver.info() != Eigen::Success)
                return ENUM_NURBS::NURBS_ERROR;
            Eigen::Matrix<T, dim, Eigen::Dynamic> new_control_points(dim, u_control_points_count);
            new_control_points.col(0) = points[index].col(0);
            new_control_points.col(u_control_points_count - 1) = points[index].col(u_points_count - 1);
            new_control_points.block(0, 1, dim, u_control_points_count - 2) = u_temp_control_points.transpose();
            temp_control_points.push_back(std::move(new_control_points));
        }

        Eigen::SparseMatrix<T> NV(v_control_points_count - 2, v_points_count - 2);
        Nrows = std::min(v_degree + 1, v_control_points_count - 2);
        NV.reserve(Eigen::VectorXi::Constant(v_points_count - 2, Nrows));
        evaulate_matrix_for_approximation<T, dim>(v_knots, v_params, v_degree, v_points_count, NV);
        mat =  NV * Eigen::SparseMatrix<T, Eigen::ColMajor>(NV.transpose());
        mat.makeCompressed();
        
        solver.compute(mat);
        if (solver.info() != Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;

        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> control_points(v_control_points_count);
        for (int index = 0; index < v_control_points_count; ++index)
            control_points[index].resize(dim, u_control_points_count);
        
        for (int index = 0; index < u_control_points_count; ++index)
        {
            Eigen::Matrix<T, dim, Eigen::Dynamic> v_temp_control_points(dim, v_points_count);
            for (int i = 0; i < v_points_count; ++i)
                v_temp_control_points.col(i) = temp_control_points[i].col(index);


            evaulate_R_matrix_for_approximation<T, dim>(v_temp_control_points, v_knots, v_params, v_degree, v_points_count, R);
            Eigen::MatrixX<T> Rs = NV * R;
            Eigen::Matrix<T, Eigen::Dynamic, dim> cps = solver.solve(Rs);
            if (solver.info() != Eigen::Success)
                return ENUM_NURBS::NURBS_ERROR;
            control_points[0].col(index) = v_temp_control_points.col(0);
            control_points[v_control_points_count - 1].col(index) = v_temp_control_points.col(v_points_count - 1);

            for (int i = 1; i < v_control_points_count - 1; ++i)
            {
                control_points[i].col(index) = cps.row(i - 1).transpose();
            }
        }
        
        nurbs.set_control_points(control_points);
        nurbs.set_uv_degree(u_degree, v_degree);
        nurbs.set_uv_knots(u_knots, v_knots);

        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename T, int dim>
    ENUM_NURBS global_curve_approximation_err_bnd(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const std::vector<T> &params, int degree, nurbs_curve<T, dim, false, -1, -1> &nurbs, T E = DEFAULT_ERROR)
    {
        int params_count = params.size();
        int points_count = points.cols();
        if (params_count != points_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        if (degree < 1)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        Eigen::VectorX<T> knots(params_count + 2);
        knots[0] = knots[1] = 0.0;
        knots[params_count] = knots[params_count + 1] = 1.0;
        for (int index = 1; index < params_count - 1; ++index)
        {
            knots[index + 1] = params[index];
        }
        std::vector<T> ek(params_count, 0.0);
        std::vector<T> current_params = params;
        nurbs_curve<T, dim, false, -1, -1> temp_nurbs(knots, points);
        nurbs_curve<T, dim, false, -1, -1> simply_nurbs;
        for (int current_degree = 1; current_degree <= degree; ++current_degree)
        {
            temp_nurbs.remove_knots_bound_curve(current_params, ek, simply_nurbs, E);
            if (current_degree == degree)
                break;
            //将temp_nubrs的节点矢量的重复的提升1
            std::vector<T> new_knots;
            knots = simply_nurbs.get_knots_vector();
            new_knots.push_back(knots[0]);
            int knots_count = knots.size();
            for (int index = 1; index < knots_count; ++index)
            {
                if (knots[index] != new_knots.back())
                {
                    new_knots.push_back(new_knots.back());
                }
                new_knots.push_back(knots[index]);
            }
            new_knots.push_back(new_knots.back());
            knots = Eigen::Map<Eigen::VectorX<T>>(new_knots.data(), new_knots.size());
            int temp_points_count = knots.size() - current_degree - 1;
            //很可能不成功
            ENUM_NURBS flag = global_least_squares_curve_approximation<T, dim>(points, current_params, knots, current_degree + 1, temp_points_count, temp_nurbs);
            if (flag != ENUM_NURBS::NURBS_SUCCESS)
            {
                global_curve_interpolate<T, dim>(points, current_degree + 1, params, temp_nurbs);
                continue;
            }
            //寻找points的最近点(TODO: 根据反求的参数点重新排序拟合点(points))
            for (int index = 0; index < points_count; ++index)
            {
                T u;
                Eigen::Vector<T, dim> R;
                temp_nurbs.find_nearst_point_on_curve(points.col(index), u, R);
                ek[index] = (points.col(index) - R).norm();
                current_params[index] = u;
            }

        }
        Eigen::VectorX<T> temp_nurbs_knots = temp_nurbs.get_knots_vector();
        Eigen::VectorX<T> simply_nurbs_knots = simply_nurbs.get_knots_vector();
        if (temp_nurbs_knots.size() == simply_nurbs_knots.size())
        {
            nurbs.set_degree(degree);
            nurbs.set_control_points(simply_nurbs.get_control_points());
            nurbs.set_knots_vector(simply_nurbs.get_knots_vector());
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        
        ENUM_NURBS flag = global_least_squares_curve_approximation<T, dim>(points, current_params, knots, degree, knots.size() - degree - 1, temp_nurbs);
        if (flag != ENUM_NURBS::NURBS_SUCCESS)
        {
            nurbs.set_degree(degree);
            nurbs.set_control_points(temp_nurbs.get_control_points());
            nurbs.set_knots_vector(temp_nurbs.get_knots_vector());
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        //else
        nurbs.set_degree(degree);
        nurbs.set_control_points(simply_nurbs.get_control_points());
        nurbs.set_knots_vector(simply_nurbs.get_knots_vector());
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim, ENPARAMETERIEDTYPE parameteried_type>
    ENUM_NURBS global_curve_approximation_err_bnd(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, int degree, nurbs_curve<T, dim, false, -1, -1> &nurbs, T E = DEFAULT_ERROR)
    {
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, params);
        }
        else
        {
            make_params_by_center<T, dim>(points, params);
        }
        return global_curve_approximation_err_bnd<T, dim>(points, params, degree, nurbs, E);
    }

    template<typename T, int dim>
    ENUM_NURBS fit_with_conic(int ks, int ke, const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Vector<T, dim> &ts, 
        const Eigen::Vector<T, dim> &te, Eigen::Vector<T, dim> &R, T &w, T E)
    {
        if (ks >= ke)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        int points_count = points.cols();
        if (ks < 0 || ke >= points_count)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        
        if (ke - ks == 1)
        {
            nurbs_curve<T, dim, true, -1, -1> new_nurbs;
            Eigen::Matrix<T, dim, Eigen::Dynamic> tangents(dim, 2);
            tangents << ts, te;
            local_2degree_arc_interpolate<T, dim>(points.block(0, ks, dim, ke - ks + 1), tangents, new_nurbs);
            // Eigen::Vector<T, dim> control_point;
            new_nurbs.get_control_point(1, R);
            // R = control_point.template block<dim, 1>(0, 0) / control_point[dim];
            // w = control_point[dim];
            new_nurbs.get_weight(1, w);
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        Eigen::Vector2<T> params;
        ENUM_NURBS flag = tow_line_intersect<T, dim>(points.col(ks), ts, points.col(ke), te, params);
        if (flag != ENUM_NURBS::NURBS_SUCCESS)
        {
            Eigen::Vector<T, dim> v = points.col(ke) - points.col(ks);
            //有bug, 需要处理v很小的情况
            v.normalize();
            for (int index = ks + 1; index < ke; ++index)
            {
                Eigen::Vector<T, dim> v2 = points.col(index) - points.col(ks);
                v2.normalize();
                T cd = v.cross(v2).squaredNorm();
                if (cd > DEFAULT_ERROR * DEFAULT_ERROR)
                    return ENUM_NURBS::NURBS_ERROR;
            }
            R = (points.col(ks) + points.col(ke)) / 2.0;
            w = 1.0;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        //else
        if (params[0] <= 0.0 || params[1] >= 0.0)
            return ENUM_NURBS::NURBS_ERROR;
        T s = 0.0;
        Eigen::Vector<T, dim> v = points.col(ke) - points.col(ks);
        R = (points.col(ks) + params[0] * ts + points.col(ke) + params[1] * te) / 2.0;
        for (int index = ks + 1; index < ke; ++index)
        {
            Eigen::Vector<T, dim> v1 = points.col(index) - R;
            flag = tow_line_intersect<T, dim>(points.col(ks), v, R, v1, params);
            if (flag != ENUM_NURBS::NURBS_SUCCESS || params[0] <= 0.0 || params[0] >= 1.0 || params[1] <= 0.0)
                return ENUM_NURBS::NURBS_ERROR;
            // Eigen::Vector<T, dim> Q = ((points.col(ks) + params[0] * v + R + params[1] * v1) / 2.0);
            
            T a = std::sqrt(params[0] / (1.0 - params[0]));
            T u = a / (1.0 + a);
            T num = std::pow(1.0 - u, 2) * (points.col(index) - points.col(ks)).dot(R - points.col(index)) + std::pow(u, 2) * (points.col(index) - points.col(ke)).dot(R - points.col(index));
            T den = 2.0 * u * (1.0 - u) * (R - points.col(index)).squaredNorm();
            T wi = num / den;
            s += wi / (1.0 + wi);
        }
        s /= (ke - ks - 1);
        w = s / (1.0 - s);
        if (w < MIN_WEIGHT<T>::value || w > MAX_WEIGHT<T>::value)
            return ENUM_NURBS::NURBS_ERROR;
        Eigen::Vector<T, dim + 1> p0, p1, p2;
        p0.template block<dim, 1>(0, 0) = points.col(ks);
        p0[dim] = 1.0;// << points.col(ks) << 1.0;
        p1.template block<dim, 1>(0, 0) = R;
        p1[dim] = 1.0;;
        p1 *= w;
        p2.template block<dim, 1>(0, 0) = points.col(ke);
        p2[dim] = 1.0;
        // p2 << points.col(ke) << 1.0;
        Eigen::Matrix<T, dim + 1, Eigen::Dynamic> control_points(dim + 1, 3);
        control_points << p0, p1, p2;
        Eigen::VectorX<T> knots(6);
        knots << 0.0, 0.0, 0.0, 1.0, 1.0, 1.0;
        nurbs_curve<T, dim, true, -1, -1> new_nurbs(knots, control_points);
        for (int index = ks + 1; index < ke; ++index)
        {
            T u;
            Eigen::Vector<T, dim> nearst_ponts;
            flag = new_nurbs.find_nearst_point_on_curve(points.col(index), u, nearst_ponts);
            if (flag != ENUM_NURBS::NURBS_SUCCESS)
                return ENUM_NURBS::NURBS_ERROR;
            T dis = (nearst_ponts - points.col(index)).norm();
            if (dis > E)
                return ENUM_NURBS::NURBS_ERROR;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim>
    ENUM_NURBS fit_with_conic(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points,  const Eigen::Matrix<T, dim, Eigen::Dynamic> &tangents, 
        nurbs_curve<T, dim, true, -1, -1> &new_nurbs, T E)
    {
        int ks = 0;
        int points_count = points.cols();
        int ke = points_count - 1;
        std::vector<Eigen::Vector<T, dim>> control_points;
        control_points.push_back(points.col(0));
        std::vector<T> wi;
        wi.push_back(1.0);
        std::vector<T> temp_knots;
        temp_knots.push_back(0.0);
        temp_knots.push_back(0.0);
        temp_knots.push_back(0.0);
        while(true)
        {
            Eigen::Vector<T, dim> R;
            T w;
            ENUM_NURBS flag = fit_with_conic<T, dim>(ks, ke, points, tangents.col(ks), tangents.col(ke), R, w, E);
            if (flag != ENUM_NURBS::NURBS_SUCCESS)
            {
                ke = (ks + ke) / 2;
            }
            else
            {
                control_points.push_back(std::move(R));
                control_points.push_back(points.col(ke));
                wi.push_back(std::move(w));
                wi.push_back(1.0);
                T chord = (points.col(ke) - points.col(ks)).norm();
                temp_knots.push_back(temp_knots.back() + chord);
                temp_knots.push_back(temp_knots.back());
                if (ke == points_count - 1)
                    break;
                ks = ke;
                ke = points_count - 1;
            }
        }
        temp_knots.push_back(temp_knots.back());
        int control_points_count = control_points.size();
        Eigen::Matrix<T, dim + 1, Eigen::Dynamic> nurbs_control_points(dim + 1, control_points_count);
        for (int index = 0; index < control_points_count; ++index)
        {
            Eigen::Vector<T, dim + 1> p;
            p.template block<dim, 1>(0, 0) = control_points[index];
            p[dim] = 1.0;
            p *= wi[index];
            nurbs_control_points.col(index) = p;
        }
        Eigen::VectorX<T> knots = Eigen::Map<Eigen::VectorX<T>>(temp_knots.data(), temp_knots.size());
        new_nurbs.set_control_points(nurbs_control_points);
        new_nurbs.set_knots_vector(knots);
        new_nurbs.set_degree(2);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim, bool flag>
    ENUM_NURBS fit_with_conic(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, nurbs_curve<T, dim, true, -1, -1> &new_nurbs, T E)
    {
        Eigen::Matrix<T, dim, Eigen::Dynamic> tangents;
        make_tangent_by_5points<T, dim, flag>(points, tangents);
        return fit_with_conic(points, tangents, new_nurbs, E);
    }


}




















