#pragma once

#include <concepts>
#include "nurbs_curve.h"
#include "nurbs_surface.h"
#include "create_nurbs_arc.h"
#include <Eigen/Sparse>

namespace tnurbs
{

    template<typename T, int dim>
    ENUM_NURBS make_params_by_chord(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, int degree, std::vector<T> &params)
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
    ENUM_NURBS make_params_by_center(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, int degree, std::vector<T> &params)
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


    template<typename T, int dim, int parameteried_type>
    ENUM_NURBS make_surface_params(const Eigen::Vector<Eigen::Matrix<T, dim, Eigen::Dynamic>, Eigen::Dynamic> &points, int u_degree, int v_degree, 
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
                ENUM_NURBS flag = make_params_by_chord(points[v_index], u_degree, u_temp_params);
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
                ENUM_NURBS flag = make_params_by_center(points[v_index], u_degree, u_temp_params);
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
                ENUM_NURBS flag = make_params_by_chord(new_control_points[u_index], v_degree, v_temp_params);
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
                ENUM_NURBS flag = make_params_by_center(new_control_points[u_index], v_degree, v_temp_params);
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

    template<typename T, int dim, int parameteried_type>
    ENUM_NURBS global_curve_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, int degree, 
            nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, degree, params);
                
        }
        else if constexpr(parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
        {
            make_params_by_center<T, dim>(points, degree, params);
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
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, degree, params);
                
        }
        else if constexpr(parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
        {
            make_params_by_center<T, dim>(points, degree, params);
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
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, 3, params);           
        }
        else if constexpr(parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
        {
            make_params_by_center<T, dim>(points, 3, params);
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

    template<typename T, int dim, int degree, int parameteried_type>
    ENUM_NURBS global_2or3degree_hermite_curve(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::Matrix<T, dim, Eigen::Dynamic> &ders,
        nurbs_curve<T, dim, false, -1, -1> &nurbs)
    {
        static_assert((degree == 2 || degree == 3), "the degree of interpolater curve is not 2 or three");
        std::vector<T> params;
        if constexpr (parameteried_type == ENPARAMETERIEDTYPE::CHORD)
        {
            make_params_by_chord<T, dim>(points, degree, params);           
        }
        else if constexpr(parameteried_type == ENPARAMETERIEDTYPE::CENTRIPETAL)
        {
            make_params_by_center<T, dim>(points, degree, params);
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

    template<typename T, int dim, int parameteried_type>
    ENUM_NURBS global_surface_interpolate(const Eigen::Vector<Eigen::Matrix<T, dim, Eigen::Dynamic>, Eigen::Dynamic> &points, int u_degree, int v_degree, 
            nurbs_surface<T, dim, -1, -1, -1, -1, false> &nurbs)
    {
        std::vector<T> u_params, v_params;
        make_surface_params<T, dim, parameteried_type>(points, u_degree, v_degree, u_params, v_params);
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

         //计算节点
        Eigen::VectorX<T> knots_vector;
        


        //计算权重
        int new_interpolate_points_count = new_interpolate_points.size();
        // std::string dir("view2.obj");
        // std::ofstream outfile(dir);
        // for (int index = 0; index < new_interpolate_points_count; ++index)
        // {
        //     Eigen::Vector3d point = new_interpolate_points[index];
        //     outfile << "v " << point[0] << " " <<
        //     point[1] << " " << point[2] << std::endl;    
        // }
        // outfile.close();

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


        local_2degree_make_params<T, dim>(new_interpolate_points, control_points, weight, knots_vector);

        std::cout << knots_vector << std::endl;
        std::cout << control_points << std::endl;
        std::cout << weight << std::endl;
        nurbs.set_control_points(control_points, weight);
        //重新参数化, 修改权重
        nurbs.set_knots_vector(knots_vector);
        nurbs.set_degree(2);
        Eigen::Matrix<T, dim + 1, Eigen::Dynamic> tempppppp = nurbs.get_control_points();
        std::cout << tempppppp << std::endl;
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 平面上的点的圆弧曲线的局部插值;推荐dim = 2时使用, 当dim > 2时, 需要保证所有的点全部在某个平面上
    /// @tparam T double float int...
    /// @tparam dim 点所在的欧式空间的维数
    /// @tparam flag 是否保留角点
    /// @param points 插值点
    /// @param nurbs 插值生成的nurbs曲线
    /// @return 错误码
    template<typename T, int dim, bool flag>
    ENUM_NURBS local_2degree_arc_interpolate(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, nurbs_curve<T, dim, true, -1, -1> &nurbs)
    {
        Eigen::Matrix<T, dim, Eigen::Dynamic> tangents;
        make_tangent_by_5points<T, dim, flag>(points, tangents);
        return local_2degree_arc_interpolate<T, dim>(points, tangents, nurbs);
    }

}




















