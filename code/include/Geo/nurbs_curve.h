#pragma once
#include "nurbs_tool.h"
#include "bezier_curve.h"

/// @brief nurbs curve class
/// @tparam T : double float int ...
/// @tparam dim : nurbs所在欧氏空间的维数
/// @tparam is_rational : 是否是有理的
/// @tparam points_count : 控制点个数
/// @tparam degree : nurbs的阶数
template<typename T, int dim, bool is_rational, int points_count, int degree>
class nurbs_curve
{
    static int constexpr order = degree + 1;
    static int constexpr rows = is_rational ? dim + 1 : dim;
    Eigen::Vector<T, degree + 1 + points_count> m_knots_vector;
    Eigen::Matrix<T, rows, points_count> m_control_points;
public:
    nurbs_curve() = default;
    nurbs_curve(const Eigen::Vector<T, degree + 1 + points_count> &knots_vector,
                const Eigen::Matrix<T, rows, points_count> &control_points) : m_knots_vector(knots_vector),
                m_control_points(control_points) { };
    
    ENUM_NURBS set_control_points(const Eigen::Matrix<T, rows, points_count> &control_points)
    {
        m_control_points = control_points;
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS set_knots_vector(const Eigen::Vector<T, degree + 1 + points_count> &knots_vector)
    {
        m_knots_vector = knots_vector;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    
    /// @brief 计算nurbs曲线上的点
    /// @param u 曲线参数
    /// @param point out_put_param nurbs曲线参数u对应的像
    /// @return ENUM_NURBS错误码
    ENUM_NURBS point_on_curve(T u, Eigen::Vector<T, dim> &point) const
    {
        int index = -1;
        find_span<T, points_count, degree>(u, m_knots_vector, index);
        Eigen::Vector<T, degree + 1> coeff;
        basis_functions<T, degree, points_count>(index, u, m_knots_vector, coeff);
        Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, order) * coeff;
        point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
    /// @param u 曲线参数
    /// @param n 求(0, 1, 2 ... n)右导数
    /// @param result out_put_pram, result[i]为nurbs曲线的第n次导数
    /// @return ENUM_NURBS错误码
    ENUM_NURBS derivative_on_curve(T u, int n, Eigen::Vector<Eigen::Vector<T, dim>, Eigen::Dynamic> &result) const
    {
        int index = -1;
        find_span<T, points_count, degree>(u, m_knots_vector, index);
        Eigen::Matrix<T, degree + 1, Eigen::Dynamic, Eigen::ColMajor, degree + 1, degree + 1> ders(degree + 1, n + 1);
        ders_basis_funs<T, degree, points_count>(index, n, u, m_knots_vector, ders);
        for (int idx = 0; idx <= n; ++idx)
        {
            Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, degree + 1) * ders.col(idx);
            result[idx] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
    /// @tparam n 求(0, 1, 2 ... n)右导数
    /// @param u 曲线参数
    /// @param result out_put_pram, result[i]为nurbs曲线的第n次导数
    /// @return ENUM_NURBS错误码
    template<int n>
    ENUM_NURBS derivative_on_curve(T u, Eigen::Vector<Eigen::Vector<T, dim>, n + 1> &result) const
    {
        int index = -1;
        find_span<T, points_count, degree>(u, m_knots_vector, index);
        Eigen::Matrix<T, degree + 1, n + 1> ders(degree + 1, n + 1);
        ders_basis_funs<T, degree, points_count, n>(index, u, m_knots_vector, ders);
        for (int idx = 0; idx <= n; ++idx)
        {
            Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, degree + 1) * ders.col(idx);
            result[idx] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 计算nurbs曲线在参数u处的0-n阶导数
    /// @tparam n 求导阶数
    /// @param u 参数u
    /// @param ders out_put_params ders(j)表示nurbs的j阶导数在参数u处的值
    /// @return ENUM_NURBS错误码
    template<int n>
    ENUM_NURBS curve_derivs_alg2(T u, Eigen::Vector<Eigen::Vector<T, dim>, n + 1> &ders) const
    {
        int span = -1;
        find_span<T, points_count, degree>(u, m_knots_vector, span);
        Eigen::Matrix<T, degree + 1, degree + 1> basis_values;
        all_basis_functions<T, degree, points_count>(span, u, m_knots_vector, basis_values);
        Eigen::Vector<Eigen::Matrix<T, rows, Eigen::Dynamic>, n + 1> PK;
        curve_deriv_cpts<T, points_count, degree, n, rows>(span - degree, span, m_knots_vector, m_control_points, PK);
        for (int k = 0; k <= n; ++k)
        {
            ders[k] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(PK[k] * basis_values.col(k));
        }

        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 计算nurbs曲线在参数u处的0-n阶导数
    /// @param n 求导阶数
    /// @param u 参数u
    /// @param ders out_put_params ders(j)表示nurbs的j阶导数在参数u处的值
    /// @return ENUM_NURBS错误码
    ENUM_NURBS curve_derivs_alg2(int n, T u, Eigen::VectorX<Eigen::Vector<T, dim>> &ders) const
    {
        int span = -1;
        find_span<T, points_count, degree>(u, m_knots_vector, span);
        Eigen::Matrix<T, degree + 1, degree + 1> basis_values;
        all_basis_functions<T, degree, points_count>(span, u, m_knots_vector, basis_values);
        Eigen::VectorX<Eigen::Matrix<T, rows, Eigen::Dynamic>> PK(n + 1);
        curve_deriv_cpts<T, points_count, degree, rows>(n, span - degree, span, m_knots_vector, m_control_points, PK);
        for (int k = 0; k <= n; ++k)
        {
            ders[k] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(PK[k] * basis_values.col(k));
        }

        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief nurbs曲线插入节点r次
    /// @tparam r 插入节点次数
    /// @param u 插入节点
    /// @return ENUM_NURBS错误码
    template<int r>
    ENUM_NURBS insert_knots(T u, nurbs_curve<T, dim, is_rational, points_count + r, degree> &new_nurbs)
    {
        constexpr int knots_size = points_count + degree + 1;
        // if (std::abs(u - m_knots_vector[knots_size - 1]) < KNOTS_VECTOR_EPS)
        if (u == m_knots_vector[knots_size - 1])
            return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
        // if (std::abs(u - m_knots_vector[0]) < KNOTS_VECTOR_EPS)
        if (u == m_knots_vector[0])
            return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
        int span = -1;
        if(find_span<T,points_count,  degree>(u, m_knots_vector,span) != ENUM_NURBS::NURBS_SUCCESS)
            return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
        int current_index = span;
        int repeat = 0;
        while (current_index >= 0)
        {
            // if (std::abs(m_knots_vector[current_index--] - u) < KNOTS_VECTOR_EPS)
            if (m_knots_vector[current_index--] == u)
            {
                ++repeat;
                continue;
            }
            break;
        }
        Eigen::Matrix<T, rows, points_count + r> new_control_points;
        curve_knots_insert<T, r, rows, degree, points_count>(m_knots_vector, m_control_points, u, span, repeat, new_control_points);
        Eigen::Vector<T, knots_size + r> new_knots_vector;
        new_knots_vector.block(0, 0, span + 1, 1) = m_knots_vector.block(0, 0, span + 1, 1);
        for (int i = 1; i <= r; ++i)
            new_knots_vector[span + i] = u;
        new_knots_vector.block(span + 1 + r, 0, knots_size - span - 1, 1) = m_knots_vector.block(span + 1, 0, knots_size - span - 1, 1);
        new_nurbs.set_control_points(new_control_points);
        new_nurbs.set_knots_vector(new_knots_vector);

        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 利用插入节点的方法求曲线上的点
    /// @param u 参数u
    /// @param point 曲面上的点
    /// @return ENUM_NURBS错误码
    ENUM_NURBS pnt_by_corner_cut(T u, Eigen::Vector<T, dim> &point)
    {
        Eigen::Vector<T, rows> ration_point;
        curve_pnt_by_corner_cut<T, rows, degree, points_count>(u, m_knots_vector, m_control_points, ration_point);
        point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(ration_point);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 节点细化
    /// @param insert_knots 插入节点数组 
    /// @return ENUM_NURBS错误码
    ENUM_NURBS refine_knots_vector(const Eigen::VectorX<T> &insert_knots, nurbs_curve<T, dim, is_rational, -1, degree> &new_nurbs)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        refine_knots_vector_curve<T, rows, degree + points_count + 1, degree>(m_knots_vector, m_control_points, insert_knots, new_knots_vector, new_control_points);
        new_nurbs.set_knots_vector(new_knots_vector);
        new_nurbs.set_control_points(new_control_points);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 将nurbs曲线分解成bezier曲线
    /// @param bezier_curves 分解的bezier曲线, 内存用户释放
    /// @return ENUM_NURBS错误码
    ENUM_NURBS decompose_to_bezier(std::vector<bezier_curve<T, dim, is_rational, degree + 1> *> &bezier_curves)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::VectorX<Eigen::Matrix<T, rows, degree + 1>> new_control_points;
        int interval_count = -1;
        find_interval_segment_count<T, points_count, degree>(m_knots_vector, interval_count);
        new_control_points.resize(interval_count);
        decompose_curve_to_bezier<T, rows, points_count,  degree>(interval_count, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
        bezier_curves.resize(interval_count, nullptr);
        for (int index = 0; index < interval_count; ++index)
        {
            bezier_curve<T, dim, is_rational, degree + 1> *cur = new bezier_curve<T, dim, is_rational, degree + 1>();
            cur->set_control_points(new_control_points[index]);
            bezier_curves[index] = cur;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 节点去除
    /// @param u 去除的节点
    /// @param count 期望去除的节点次数
    /// @param time 实际去除是节点次数
    /// @param new_nurbs 去除节点新生成的nurbs
    /// @param error 去点去除的误差
    /// @return ENUM_NURBS错误码
    ENUM_NURBS remove_knots(T u, int count, int &time, nurbs_curve<T, dim, is_rational, -1, degree> &new_nurbs, T error = DEFAULT_ERROR)
    {
        constexpr int knots_size = degree + 1 + points_count;
        if (u == m_knots_vector[0] || u == m_knots_vector[knots_size - 1])
        {
            time = 0;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        int span = -1;
        find_span<T, points_count, degree>(u, m_knots_vector, span);
        if (m_knots_vector[span] != u)
            return ENUM_NURBS::NURBS_ERROR;
        
        int repeat = 1;
        int current = span - 1;
        while (current >= 0 && m_knots_vector[current] == u)
        {
            repeat += 1;
            current -= 1;
        }
        count = std::min(repeat, count);
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        remove_curve_knots<T, degree, points_count, rows, is_rational>(span, count, repeat, m_knots_vector, m_control_points,
                   new_knots_vector, new_control_points, time, error);
        new_nurbs.set_control_points(new_control_points);
        new_nurbs.set_knots_vector(new_knots_vector);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 升阶
    /// @param t 提升的阶数
    /// @return ENUM_NURBS错误码
    ENUM_NURBS degree_elevate(int t, nurbs_curve<T, dim, is_rational, -1, -1> &curve)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        degree_elevate_curve<T, rows>(t, degree, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
        
        curve.set_control_points(new_control_points);
        curve.set_knots_vector(new_knots_vector);
        curve.set_degree(degree + t);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

};

template<typename T, int dim, bool is_rational, int degree>
class nurbs_curve<T, dim, is_rational, -1, degree>
{
    static int constexpr order = degree + 1;
    static int constexpr rows = is_rational ? dim + 1 : dim;
    Eigen::Vector<T, Eigen::Dynamic> m_knots_vector; //degree + 1 + points_count
    Eigen::Matrix<T, rows, Eigen::Dynamic> m_control_points; //(rows, points_count)
public:
    nurbs_curve() = default;
    nurbs_curve(const Eigen::Vector<T, Eigen::Dynamic> &knots_vector,
                const Eigen::Matrix<T, rows, Eigen::Dynamic> &control_points) : m_knots_vector(knots_vector),
                m_control_points(control_points) { };
    
    ENUM_NURBS set_control_points(const Eigen::Matrix<T, rows, Eigen::Dynamic> &control_points)
    {
        m_control_points = control_points;
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS set_knots_vector(const Eigen::VectorX<T> &knots_vector)
    {
        m_knots_vector = knots_vector;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    
    /// @brief 计算nurbs曲线上的点
    /// @param u 曲线参数
    /// @param point out_put_param nurbs曲线参数u对应的像
    /// @return ENUM_NURBS错误码
    ENUM_NURBS point_on_curve(T u, Eigen::Vector<T, dim> &point) const
    {
        int index = -1;
        find_span<T, degree>(u, m_knots_vector, index);
        Eigen::Vector<T, degree + 1> coeff;
        basis_functions<T, degree>(index, u, m_knots_vector, coeff);
        Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, order) * coeff;
        point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    
    /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
    /// @param u 曲线参数
    /// @param n 求(0, 1, 2 ... n)右导数
    /// @param result out_put_pram, result[i]为nurbs曲线的第n次导数
    /// @return ENUM_NURBS错误码
    ENUM_NURBS derivative_on_curve(T u, int n, Eigen::Vector<Eigen::Vector<T, dim>, Eigen::Dynamic> &result) const
    {
        int index = -1;
        find_span<T, degree>(u, m_knots_vector, index);
        Eigen::Matrix<T, degree + 1, Eigen::Dynamic, Eigen::ColMajor, degree + 1, degree + 1> ders(degree + 1, n + 1);
        ders_basis_funs<T, degree>(index, n, u, m_knots_vector, ders);
        for (int idx = 0; idx <= n; ++idx)
        {
            Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, degree + 1) * ders.col(idx);
            result[idx] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
    /// @tparam n 求(0, 1, 2 ... n)右导数
    /// @param u 曲线参数
    /// @param result out_put_pram, result[i]为nurbs曲线的第n次导数
    /// @return ENUM_NURBS错误码
    template<int n>
    ENUM_NURBS derivative_on_curve(T u, Eigen::Vector<Eigen::Vector<T, dim>, Eigen::Dynamic> &result) const
    {
        int index = -1;
        find_span<T, degree>(u, m_knots_vector, index);
        Eigen::Matrix<T, degree + 1, n + 1> ders(degree + 1, n + 1);
        ders_basis_funs<T, degree, n>(index, u, m_knots_vector, ders);
        for (int idx = 0; idx <= n; ++idx)
        {
            Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, degree + 1) * ders.col(idx);
            result[idx] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 计算nurbs曲线在参数u处的0-n阶导数
    /// @tparam n 求导阶数
    /// @param u 参数u
    /// @param ders out_put_params ders(j)表示nurbs的j阶导数在参数u处的值
    /// @return ENUM_NURBS错误码
    template<int n>
    ENUM_NURBS curve_derivs_alg2(T u, Eigen::Vector<Eigen::Vector<T, dim>, n + 1> &ders) const
    {
        int span = -1;
        find_span<T, degree>(u, m_knots_vector, span);
        Eigen::Matrix<T, degree + 1, degree + 1> basis_values;
        all_basis_functions<T, degree>(span, u, m_knots_vector, basis_values);
        Eigen::Vector<Eigen::Matrix<T, rows, Eigen::Dynamic>, n + 1> PK;
        curve_deriv_cpts<T, degree, n, rows>(span - degree, span, m_knots_vector, m_control_points, PK);
        for (int k = 0; k <= n; ++k)
        {
            ders[k] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(PK[k] * basis_values.col(k));
        }

        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 计算nurbs曲线在参数u处的0-n阶导数
    /// @param n 求导阶数
    /// @param u 参数u
    /// @param ders out_put_params ders(j)表示nurbs的j阶导数在参数u处的值
    /// @return ENUM_NURBS错误码
    ENUM_NURBS curve_derivs_alg2(int n, T u, Eigen::VectorX<Eigen::Vector<T, dim>> &ders) const
    {
        int span = -1;
        find_span<T, degree>(u, m_knots_vector, span);
        Eigen::Matrix<T, degree + 1, degree + 1> basis_values;
        all_basis_functions<T, degree>(span, u, m_knots_vector, basis_values);
        Eigen::VectorX<Eigen::Matrix<T, rows, Eigen::Dynamic>> PK(n + 1);
        curve_deriv_cpts<T, degree, rows>(n, span - degree, span, m_knots_vector, m_control_points, PK);
        for (int k = 0; k <= n; ++k)
        {
            ders[k] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(PK[k] * basis_values.col(k));
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief nurbs曲线插入节点r次
    /// @tparam r 插入节点次数
    /// @param u 插入节点
    /// @return ENUM_NURBS错误码
    ENUM_NURBS insert_knots(T u, int r)
    {
        int span = -1;
        int knots_size = m_knots_vector.size();
        if(find_span<T>(u, degree, m_knots_vector,span) != ENUM_NURBS::NURBS_SUCCESS)
            return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
        int current_index = span;
        // if (std::abs(u - m_knots_vector[knots_size - 1]) < KNOTS_VECTOR_EPS)
        if (u == m_knots_vector[knots_size - 1])
            return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
        // if (std::abs(u - m_knots_vector[0]) < KNOTS_VECTOR_EPS)
        if (u == m_knots_vector[0])
            return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
        int repeat = 0;
        while (current_index >= 0)
        {
            if (m_knots_vector[current_index--] == u)
            // if (std::abs(m_knots_vector[current_index--] - u) < KNOTS_VECTOR_EPS)
            {
                ++repeat;
                continue;
            }
            break;
        }
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        curve_knots_insert<T, rows>(r, m_control_points.cols(), degree, m_knots_vector, m_control_points, u, span, repeat, new_control_points);
        m_control_points = new_control_points;
        Eigen::Vector<T, Eigen::Dynamic> new_knots_vector(knots_size + r);
        new_knots_vector.block(0, 0, span + 1, 1) = m_knots_vector.block(0, 0, span + 1, 1);
        for (int i = 1; i <= r; ++i)
            new_knots_vector[span + i] = u;
        new_knots_vector.block(span + 1 + r, 0, knots_size - span - 1, 1) = m_knots_vector.block(span + 1, 0, knots_size - span - 1, 1);
        m_knots_vector = new_knots_vector;
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 利用插入节点的方法求曲线上的点
    /// @param u 参数u
    /// @param point 曲面上的点
    /// @return ENUM_NURBS错误码
    ENUM_NURBS pnt_by_corner_cut(T u, Eigen::Vector<T, dim> &point)
    {
        int control_points_count = m_control_points.cols();
        Eigen::Vector<T, rows> ration_point;
        curve_pnt_by_corner_cut<T, rows, degree>(u, control_points_count, m_knots_vector, m_control_points, ration_point);
        point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(ration_point);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 节点细化
    /// @param insert_knots 插入节点数组 
    /// @return ENUM_NURBS错误码
    ENUM_NURBS refine_knots_vector(const Eigen::VectorX<T> &insert_knots)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        int knots_size = m_knots_vector.size();
        refine_knots_vector_curve<T, rows>(knots_size, degree, m_knots_vector, m_control_points, insert_knots, new_knots_vector, new_control_points);
        m_knots_vector = new_knots_vector;
        m_control_points = new_control_points;
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 将nurbs曲线分解成bezier曲线
    /// @param bezier_curves 分解的bezier曲线, 内存用户释放
    /// @return ENUM_NURBS错误码
    ENUM_NURBS decompose_to_bezier(std::vector<bezier_curve<T, dim, is_rational, degree + 1> *> &bezier_curves)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::VectorX<Eigen::Matrix<T, rows, degree + 1>> new_control_points;
        int interval_count = -1;
        find_interval_segment_count<T>(degree, m_knots_vector, interval_count);
        new_control_points.resize(interval_count);
        decompose_curve_to_bezier<T, rows, degree>(interval_count, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
        // int curve_count = new_control_points.size();
        bezier_curves.resize(interval_count, nullptr);
        for (int index = 0; index < interval_count; ++index)
        {
            bezier_curve<T, dim, is_rational, degree + 1> *cur = new bezier_curve<T, dim, is_rational, degree + 1>();
            cur->set_control_points(new_control_points[index]);
            bezier_curves[index] = cur;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 节点去除
    /// @param u 去除的节点
    /// @param count 期望去除的节点次数
    /// @param time 实际去除是节点次数
    /// @param error 去点去除的误差
    /// @return ENUM_NURBS错误码
    ENUM_NURBS remove_knots(T u, int count, int &time, T error = DEFAULT_ERROR)
    {
        int knots_size = m_knots_vector.size();
        if (u == m_knots_vector[0] || u == m_knots_vector[knots_size - 1])
        {
            time = 0;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        int span = -1;
        find_span<T, degree>(u, m_knots_vector, span);
        if (m_knots_vector[span] != u)
            return ENUM_NURBS::NURBS_ERROR;
        
        int repeat = 1;
        int current = span - 1;
        while (current >= 0 && m_knots_vector[current] == u)
        {
            repeat += 1;
            current -= 1;
        }
        count = std::min(repeat, count);
        remove_curve_knots<T, rows, is_rational>(span, count, degree, repeat, m_knots_vector, m_control_points, time, error);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 升阶
    /// @param t 提升的阶数
    /// @return ENUM_NURBS错误码
    ENUM_NURBS degree_elevate(int t, nurbs_curve<T, dim, is_rational, -1, -1> &curve)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        degree_elevate_curve<T, rows>(t, degree, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
        
        curve.set_control_points(new_control_points);
        curve.set_knots_vector(new_knots_vector);
        curve.set_degree(degree + t);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

};

/// @brief 
/// @tparam T 
/// @tparam dim 
/// @tparam is_rational 
template<typename T, int dim, bool is_rational>
class nurbs_curve<T, dim, is_rational, -1, -1>
{
    static int constexpr rows = is_rational ? dim + 1 : dim;
    int m_degree;
    Eigen::Vector<T, Eigen::Dynamic> m_knots_vector; //degree + 1 + points_count
    Eigen::Matrix<T, rows, Eigen::Dynamic> m_control_points; //(rows, points_count)
public:
    nurbs_curve() = default;
    nurbs_curve(const Eigen::Vector<T, Eigen::Dynamic> &knots_vector,
                const Eigen::Matrix<T, rows, Eigen::Dynamic> &control_points) : m_knots_vector(knots_vector),
                m_control_points(control_points) { m_degree = m_knots_vector.size() - m_control_points.cols() - 1; }
    

    ENUM_NURBS set_degree(int degree) { m_degree = degree; return ENUM_NURBS::NURBS_SUCCESS; }

    ENUM_NURBS set_control_points(const Eigen::Matrix<T, rows, Eigen::Dynamic> &points) { m_control_points = points; return ENUM_NURBS::NURBS_SUCCESS; }

    ENUM_NURBS set_knots_vector(const Eigen::VectorX<T> &knots_vector) { m_knots_vector = knots_vector; return ENUM_NURBS::NURBS_SUCCESS; }


    /// @brief 计算nurbs曲线上的点
    /// @param u 曲线参数
    /// @param point out_put_param nurbs曲线参数u对应的像
    /// @return ENUM_NURBS错误码
    ENUM_NURBS point_on_curve(T u, Eigen::Vector<T, dim> &point) const
    {
        int index = -1;
        find_span<T>(u,m_degree, m_knots_vector, index);
        Eigen::Vector<T, Eigen::Dynamic> coeff(m_degree + 1);
        basis_functions<T>(index, u, m_degree, m_knots_vector, coeff);
        Eigen::Vector<T, rows> vec = m_control_points.block(0, index - m_degree, rows, m_degree + 1) * coeff;
        point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

        /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
    /// @param u 曲线参数
    /// @param n 求(0, 1, 2 ... n)右导数
    /// @param result out_put_pram, result[i]为nurbs曲线的第n次导数
    /// @return ENUM_NURBS错误码
    ENUM_NURBS derivative_on_curve(T u, int n, Eigen::Vector<Eigen::Vector<T, dim>, Eigen::Dynamic> &result) const
    {
        int index = -1;
        find_span<T>(u, m_degree, m_knots_vector, index);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ders(m_degree + 1, n + 1);
        ders_basis_funs<T>(index, n, m_degree, u, m_knots_vector, ders);
        for (int idx = 0; idx <= n; ++idx)
        {
            Eigen::Vector<T, rows> vec = m_control_points.block(0, index - m_degree, rows,m_degree + 1) * ders.col(idx);
            result[idx] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
    /// @tparam n 求(0, 1, 2 ... n)右导数
    /// @param u 曲线参数
    /// @param result out_put_pram, result[i]为nurbs曲线的第n次导数
    /// @return ENUM_NURBS错误码
    template<int n>
    ENUM_NURBS derivative_on_curve(T u, Eigen::Vector<Eigen::Vector<T, dim>, Eigen::Dynamic> &result) const
    {
        int index = -1;
        find_span<T>(u, m_degree, m_knots_vector, index);
        Eigen::Matrix<T, Eigen::Dynamic, n + 1> ders(m_degree + 1, n + 1);
        ders_basis_funs<T, n>(index, m_degree, u, m_knots_vector, ders);
        for (int idx = 0; idx <= n; ++idx)
        {
            Eigen::Vector<T, rows> vec = m_control_points.block(0, index - m_degree, rows, m_degree + 1) * ders.col(idx);
            result[idx] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 计算nurbs曲线在参数u处的0-n阶导数
    /// @tparam n 求导阶数
    /// @param u 参数u
    /// @param ders out_put_params ders(j)表示nurbs的j阶导数在参数u处的值
    /// @return ENUM_NURBS错误码
    template<int n>
    ENUM_NURBS curve_derivs_alg2(T u, Eigen::Vector<Eigen::Vector<T, dim>, n + 1> &ders) const
    {
        int span = -1;
        find_span<T>(u, m_degree, m_knots_vector, span);
        Eigen::MatrixX<T> basis_values;
        all_basis_functions<T>(span, m_degree, u, m_knots_vector, basis_values);
        Eigen::Vector<Eigen::Matrix<T, rows, Eigen::Dynamic>, n + 1> PK;
        curve_deriv_cpts<T, n, rows>(m_degree, span - m_degree, span, m_knots_vector, m_control_points, PK);
        for (int k = 0; k <= n; ++k)
        {
            ders[k] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(PK[k] * basis_values.col(k));
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 计算nurbs曲线在参数u处的0-n阶导数
    /// @param n 求导阶数
    /// @param u 参数u
    /// @param ders out_put_params ders(j)表示nurbs的j阶导数在参数u处的值
    /// @return ENUM_NURBS错误码
    ENUM_NURBS curve_derivs_alg2(int n, T u, Eigen::VectorX<Eigen::Vector<T, dim>> &ders) const
    {
        int span = -1;
        find_span<T>(u,m_degree,  m_knots_vector, span);
        Eigen::MatrixX<T> basis_values;
        all_basis_functions<T>(span, m_degree, u, m_knots_vector, basis_values);
        Eigen::VectorX<Eigen::Matrix<T, rows, Eigen::Dynamic>> PK(n + 1);
        curve_deriv_cpts<T, rows>(m_degree, n, span - m_degree, span, m_knots_vector, m_control_points, PK);
        for (int k = 0; k <= n; ++k)
        {
            ders[k] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(PK[k] * basis_values.col(k));
        }        
        
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    // /// @brief nurbs曲线插入节点r次
    // /// @tparam r 插入节点次数
    // /// @param u 插入节点
    // /// @return ENUM_NURBS错误码
    // template<int r>
    // ENUM_NURBS insert_knots(T u)
    // {
    //     int span = -1;
    //     int knots_size = m_knots_vector.size();
    //     if(find_span<T>(u, m_degree, m_knots_vector,span) != ENUM_NURBS::NURBS_SUCCESS)
    //         return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
    //     int current_index = span;
    //     if (std::abs(u - m_knots_vector[knots_size - 1]) < KNOTS_VECTOR_EPS)
    //         return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
    //     if (std::abs(u - m_knots_vector[0]) < KNOTS_VECTOR_EPS)
    //         return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
    //     int repeat = 0;
    //     while (current_index >= 0)
    //     {
    //         if (std::abs(m_knots_vector[current_index--] - u) < KNOTS_VECTOR_EPS)
    //         {
    //             ++repeat;
    //             continue;
    //         }
    //         break;
    //     }
    //     Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
    //     curve_knots_insert<T, r, rows>(m_control_points.cols(), m_degree, m_knots_vector, m_control_points, u, span, repeat, new_control_points);
    //     m_control_points = new_control_points;
    //     Eigen::Vector<T, Eigen::Dynamic> new_knots_vector(knots_size + r);
    //     new_knots_vector.block(0, 0, span + 1, 1) = m_knots_vector.block(0, 0, span + 1, 1);
    //     for (int i = 1; i <= r; ++i)
    //         new_knots_vector[span + i] = u;
    //     new_knots_vector.block(span + 1 + r, 0, knots_size - span - 1, 1) = m_knots_vector.block(span + 1, 0, knots_size - span - 1, 1);
    //     m_knots_vector = new_knots_vector;
    //     return ENUM_NURBS::NURBS_SUCCESS;
    // }

    /// @brief nurbs曲线插入节点r次
    /// @tparam r 插入节点次数
    /// @param u 插入节点
    /// @return ENUM_NURBS错误码
    ENUM_NURBS insert_knots(T u, int r)
    {
        int span = -1;
        int knots_size = m_knots_vector.size();
        if(find_span<T>(u, m_degree, m_knots_vector,span) != ENUM_NURBS::NURBS_SUCCESS)
            return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
        int current_index = span;
        // if (std::abs(u - m_knots_vector[knots_size - 1]) < KNOTS_VECTOR_EPS)
        if (u == m_knots_vector[knots_size - 1])
            return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
        // if (std::abs(u - m_knots_vector[0]) < KNOTS_VECTOR_EPS)
        if (u == m_knots_vector[0])
            return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
        int repeat = 0;
        while (current_index >= 0)
        {
            // if (std::abs(m_knots_vector[current_index--] - u) < KNOTS_VECTOR_EPS)
            if (m_knots_vector[current_index--] == u)
            {
                ++repeat;
                continue;
            }
            break;
        }
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        curve_knots_insert<T, rows>(r, m_control_points.cols(), m_degree, m_knots_vector, m_control_points, u, span, repeat, new_control_points);
        m_control_points = new_control_points;
        Eigen::Vector<T, Eigen::Dynamic> new_knots_vector(knots_size + r);
        new_knots_vector.block(0, 0, span + 1, 1) = m_knots_vector.block(0, 0, span + 1, 1);
        for (int i = 1; i <= r; ++i)
            new_knots_vector[span + i] = u;
        new_knots_vector.block(span + 1 + r, 0, knots_size - span - 1, 1) = m_knots_vector.block(span + 1, 0, knots_size - span - 1, 1);
        m_knots_vector = new_knots_vector;
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 利用插入节点的方法求曲线上的点
    /// @param u 参数u
    /// @param point 曲面上的点
    /// @return ENUM_NURBS错误码
    ENUM_NURBS pnt_by_corner_cut(T u, Eigen::Vector<T, dim> &point) const
    {
        int control_points_count = m_control_points.cols();
        Eigen::Vector<T, rows> ration_point;
        curve_pnt_by_corner_cut<T, rows>(u, control_points_count, m_degree, m_knots_vector, m_control_points, ration_point);
        point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(ration_point);
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 节点细化
    /// @param insert_knots 插入节点数组 
    /// @return ENUM_NURBS错误码
    ENUM_NURBS refine_knots_vector(const Eigen::VectorX<T> &insert_knots)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        int knots_size = m_knots_vector.size();
        refine_knots_vector_curve<T, rows>(knots_size, m_degree, m_knots_vector, m_control_points, insert_knots, new_knots_vector, new_control_points);
        m_knots_vector = new_knots_vector;
        m_control_points = new_control_points;
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 将nurbs曲线分解成bezier曲线
    /// @param bezier_curves 分解的bezier曲线, 内存用户释放
    /// @return ENUM_NURBS错误码
    ENUM_NURBS decompose_to_bezier(std::vector<bezier_curve<T, dim, is_rational, -1> *> &bezier_curves)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::VectorX<Eigen::Matrix<T, rows, Eigen::Dynamic>> new_control_points;
        int interval_count = -1;
        find_interval_segment_count<T>(m_degree, m_knots_vector, interval_count);
        new_control_points.resize(interval_count);
        decompose_curve_to_bezier<T, rows>(m_degree, interval_count, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
        // int interval_count = new_control_points.size();
        bezier_curves.resize(interval_count, nullptr);
        for (int index = 0; index < interval_count; ++index)
        {
            bezier_curve<T, dim, is_rational, -1> *cur = new bezier_curve<T, dim, is_rational, -1>();
            cur->set_control_points(new_control_points[index]);
            bezier_curves[index] = cur;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 节点去除
    /// @param u 去除的节点
    /// @param count 期望去除的节点次数
    /// @param time 实际去除是节点次数
    /// @param error 去点去除的误差
    /// @return ENUM_NURBS错误码
    ENUM_NURBS remove_knots(T u, int count, int &time, T error = DEFAULT_ERROR)
    {
        int knots_size = m_knots_vector.size();
        if (u == m_knots_vector[0] || u == m_knots_vector[knots_size - 1])
        {
            time = 0;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        int span = -1;
        find_span<T>(u, m_degree, m_knots_vector, span);
        if (m_knots_vector[span] != u)
            return ENUM_NURBS::NURBS_ERROR;
        
        int repeat = 1;
        int current = span - 1;
        while (current >= 0 && m_knots_vector[current] == u)
        {
            repeat += 1;
            current -= 1;
        }
        count = std::min(repeat, count);
        remove_curve_knots<T, rows, is_rational>(span, count, m_degree, repeat, m_knots_vector, m_control_points, time, error);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 升阶
    /// @param t 提升的阶数
    /// @return ENUM_NURBS错误码
    ENUM_NURBS degree_elevate(int t)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        degree_elevate_curve<T, rows>(t, m_degree, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
        m_degree += t;
        m_knots_vector = new_knots_vector;
        m_control_points = new_control_points;
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS degree_reduce(T error = DEFAULT_ERROR)
    {
        Eigen::VectorX<T> new_knots_vector(m_knots_vector.size() - 1);
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        ENUM_NURBS flag = degree_reduce_curve<T, rows, is_rational>(m_degree, m_knots_vector, m_control_points, new_knots_vector, new_control_points, error);
        if (flag != ENUM_NURBS::NURBS_SUCCESS)
            return flag;
        m_control_points = new_control_points;
        m_knots_vector = new_knots_vector;
        m_degree -= 1;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
};


