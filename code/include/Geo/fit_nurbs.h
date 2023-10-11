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




}




















