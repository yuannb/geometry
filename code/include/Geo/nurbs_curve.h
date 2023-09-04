#pragma once
#include "nurbs_tool.h"
#include "bezier_curve.h"
#include <array>
// #include "ThreadPool.h"

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

    int is_include_konts(T u, T eps = 0.0)
    {
        int knots_size = m_knots_vector.size();
        for (int index = 0; index < knots_size; ++index)
        {
            if (std::abs(u - m_knots_vector[index]) <= eps)
                return index;
        }
        return INDEX_IS_OUTSIDE_OF_KNOTS_VECTOR;
    }

    int get_knots_count() const
    {
        return m_knots_vector.size();
    }

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
    

    int get_knots_multiplicity(int index) const
    {
        int konts_count = m_knots_vector.size();
        T knot = m_knots_vector[index];
        if (index >= konts_count || index < 0)
            return 0;
        int count = 1;
        int current_index = index - 1;
        while (current_index)
        {
            if (m_knots_vector[current_index] == knot)
            {
                ++count;
                --current_index;
            }
            else
                break;
        }
        current_index = index + 1;
        while (current_index < konts_count)
        {
            if (m_knots_vector[current_index] == knot)
            {
                ++count;
                ++current_index;
            }
            else
                break;
        }
        return count;
    }

    ENUM_NURBS get_ends_point(std::array<Eigen::Vector<T, dim>, 2> &points) const
    {
        points[0] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
        // int points_count = m_control_points.cols();
        points[1] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(points_count - 1));
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
            // std::cout << new_control_points[index] << std::endl;
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

    /// @brief 升阶
    /// @tparam t 提升的阶数
    /// @param curve 提升后的曲线
    /// @return ENUM_NURBS错误码
    template<int t>
    ENUM_NURBS degree_elevate(nurbs_curve<T, dim, is_rational, -1, degree + t> &curve)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        degree_elevate_curve<T, rows>(t, degree, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
        
        curve.set_control_points(new_control_points);
        curve.set_knots_vector(new_knots_vector);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 将nurbs曲线降低一阶
    /// @param error 误差
    /// @param curve 降阶后的曲线
    /// @return ENUM_NURBS错误码
    ENUM_NURBS degree_reduce(nurbs_curve<T, dim, is_rational, -1, degree - 1> &curve, T error = DEFAULT_ERROR)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        ENUM_NURBS flag = degree_reduce_curve<T, rows, is_rational>(degree, m_knots_vector, m_control_points, new_knots_vector, new_control_points, error);
        if (flag != ENUM_NURBS::NURBS_SUCCESS)
            return flag;
        curve.set_control_points(new_control_points);
        curve.set_knots_vector(new_knots_vector);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    bool is_closed() const
    {
        Eigen::Vector<T, dim> start = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
        Eigen::Vector<T, dim> end = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(points_count - 1));
        return (start - end).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR;
    }

    //TODO: 写个多线程, 多线程每次结果可能都会不同
    ENUM_NURBS find_nearst_point_on_curve(const Eigen::Vector<T, dim> &point, T &u, Eigen::Vector<T, dim> &nearst_point, T eps = DEFAULT_ERROR) const
    {
        T min = m_knots_vector[0];
        int knots_size = m_knots_vector.size();
        T max = m_knots_vector(knots_size - 1);

        //将节点分成(max - min + 1) * 100份；以后需要优化
        int step_count = (max - min + 1) * 10;
        T step = (max - min) / step_count;
        T min_length = (point - m_control_points.col(0)).squaredNorm();
        nearst_point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
        u = m_knots_vector[0];
        Eigen::Vector<Eigen::Vector<T, dim>, 3> ders_vec;
        for (int index = 0; index <= step_count; ++index)
        {
            T distance;
            T current_u = index * step + min;
            for (int loop_index = 0; loop_index < MAX_ITERATE_DEEP; ++loop_index)
            {
                curve_derivs_alg2<2>(current_u, ders_vec);
                Eigen::Vector<T, dim> dist_vec = ders_vec[0] - point;
                distance = dist_vec.squaredNorm();
                if (distance < eps * eps)
                {
                    u = current_u;
                    nearst_point = ders_vec[0];
                    return ENUM_NURBS::NURBS_SUCCESS;
                }

                if (min_length > distance)
                {
                    min_length = distance;
                    u = current_u;
                    nearst_point = ders_vec[0];
                }

                T cos_angle = ders_vec[1].dot(dist_vec);
                T tanget_vec_square_len = ders_vec[1].squaredNorm();
                if (cos_angle * cos_angle < tanget_vec_square_len * distance * eps * eps)
                {
                    break;
                }
            
                T next_u = current_u - cos_angle / (ders_vec[2].dot(dist_vec) + tanget_vec_square_len);
                bool is_closed_flag = is_closed();
                if (is_closed_flag)
                {
                    if (next_u < min)
                        next_u = max - (min - next_u);
                    else if (next_u > max)
                        next_u = min + (next_u - max);
                }
                else
                {
                    if (next_u < min)
                        next_u = min;
                    else if (next_u > max)
                        next_u = max;
                }
                if ((next_u - current_u) * (next_u - current_u) * tanget_vec_square_len < eps * eps)
                    break;
                current_u = next_u;
            }

        }
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
    
    T get_knot(int index) const
    {
        return m_knots_vector[index];
    }

    int get_knots_count() const
    {
        return m_knots_vector.size();
    }
    
    ENUM_NURBS get_ends_point(std::array<Eigen::Vector<T, dim>, 2> &points) const
    {
        points[0] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
        int points_count = m_control_points.cols();
        points[1] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(points_count - 1));
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    int get_knots_multiplicity(int index) const
    {
        int konts_count = m_knots_vector.size();
        T knot = m_knots_vector[index];
        if (index >= konts_count || index < 0)
            return 0;
        int count = 1;
        int current_index = index - 1;
        while (current_index)
        {
            if (m_knots_vector[current_index] == knot)
            {
                ++count;
                --current_index;
            }
            else
                break;
        }
        current_index = index + 1;
        while (current_index < konts_count)
        {
            if (m_knots_vector[current_index] == knot)
            {
                ++count;
                ++current_index;
            }
            else
                break;
        }
        return count;
    }

    int is_include_konts(T u, T eps = 0.0)
    {
        int knots_size = m_knots_vector.size();
        for (int index = 0; index < knots_size; ++index)
        {
            if (std::abs(u - m_knots_vector[index]) <= eps)
                return index;
        }
        return INDEX_IS_OUTSIDE_OF_KNOTS_VECTOR;
    }


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
            // std::cout << new_control_points[index] << std::endl;
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
    /// @param curve 升阶后的曲线
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

    /// @brief 升阶
    /// @param t 提升的阶数
    /// @param curve 升阶后的曲线
    /// @return ENUM_NURBS错误码
    template<int t>
    ENUM_NURBS degree_elevate(nurbs_curve<T, dim, is_rational, -1, degree + t> &curve)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        degree_elevate_curve<T, rows>(t, degree, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
        curve.set_control_points(new_control_points);
        curve.set_knots_vector(new_knots_vector);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 将nurbs曲线降低一阶
    /// @param error 误差
    /// @param curve 降阶后的曲线
    /// @return ENUM_NURBS错误码
    ENUM_NURBS degree_reduce(nurbs_curve<T, dim, is_rational, -1, degree - 1> &curve, T error = DEFAULT_ERROR)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        ENUM_NURBS flag = degree_reduce_curve<T, rows, is_rational>(degree, m_knots_vector, m_control_points, new_knots_vector, new_control_points, error);
        if (flag != ENUM_NURBS::NURBS_SUCCESS)
            return flag;
        curve.set_control_points(new_control_points);
        curve.set_knots_vector(new_knots_vector);
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    bool is_closed() const
    {
        int points_count = m_control_points.cols();
        Eigen::Vector<T, dim> start = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
        Eigen::Vector<T, dim> end = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(points_count - 1));
        return (start - end).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR;
    }


    //TODO: 写个多线程, 多线程每次结果可能都会不同
    ENUM_NURBS find_nearst_point_on_curve(const Eigen::Vector<T, rows> &point, T &u, Eigen::Vector<T, rows> &nearst_point, T eps = DEFAULT_ERROR) const
    {
        T min = m_knots_vector[0];
        int knots_size = m_knots_vector.size();
        T max = m_knots_vector(knots_size - 1);

        //将节点分成(max - min + 1) * 100份；以后需要优化
        int step_count = (max - min + 1) * 10;
        T step = (max - min) / step_count;
        T min_length = (point - m_control_points.col(0)).squaredNorm();
        nearst_point = nearst_point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
        u = m_knots_vector[0];
        Eigen::Vector<Eigen::Vector<T, dim>, 3> ders_vec;
        for (int index = 0; index <= step_count; ++index)
        {
            T distance;
            T current_u = index * step + min;
            for (int loop_index = 0; loop_index < MAX_ITERATE_DEEP; ++loop_index)
            {
                curve_derivs_alg2<2>(current_u, ders_vec);
                Eigen::Vector<T, dim> dist_vec = ders_vec[0] - point;
                distance = dist_vec.squaredNorm();
                if (distance < eps * eps)
                {
                    u = current_u;
                    nearst_point = ders_vec[0];
                    return ENUM_NURBS::NURBS_SUCCESS;
                }
                if (min_length > distance)
                {
                    min_length = distance;
                    u = current_u;
                    nearst_point = ders_vec[0];
                }

                T cos_angle = ders_vec[1].dot(dist_vec);
                T tanget_vec_square_len = ders_vec[1].squaredNorm();
                if (cos_angle * cos_angle < tanget_vec_square_len * distance * eps * eps)
                {
                    break;
                }
            
                T next_u = current_u - cos_angle / (ders_vec[2].dot(dist_vec) + tanget_vec_square_len);
                bool is_closed_flag = is_closed();
                if (is_closed_flag)
                {
                    if (next_u < min)
                        next_u = max - (min - next_u);
                    else if (next_u > max)
                        next_u = min + (next_u - max);
                }
                else
                {
                    if (next_u < min)
                        next_u = min;
                    else if (next_u > max)
                        next_u = max;
                }
                if ((next_u - current_u) * (next_u - current_u) * tanget_vec_square_len < eps * eps)
                    break;
                current_u = next_u;
            }
        
            // if (min_length > distance)
            // {
            //     min_length = distance;
            //     u = current_u;
            //     nearst_point = ders_vec[0];
            // }
        }
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
                m_control_points(control_points) 
                { 
                    // int knots_size = m_knots_vector.size();
                    // int control_points_count = m_control_points.cols();
                    m_degree = m_knots_vector.size() - m_control_points.cols() - 1; 
                }
    nurbs_curve(const nurbs_curve<T, dim, is_rational, -1, -1> &nurbs_curve_to_copy)
    {
        m_knots_vector = nurbs_curve_to_copy.get_knots_vector();
        m_degree = nurbs_curve_to_copy.get_degree();
        m_control_points = nurbs_curve_to_copy.get_control_points();
    }

    ENUM_NURBS get_ends_point(std::array<Eigen::Vector<T, dim>, 2> &points) const
    {
        points[0] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
        int points_count = m_control_points.cols();
        // std::cout << m_control_points << std::endl;
        points[1] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(points_count - 1));
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    int get_knots_count() const
    {
        return m_knots_vector.size();
    }

    int get_next_different_index(int index) const
    {
        T knot = m_knots_vector[index];
        int knots_count = m_knots_vector.size();
        for (int i = index + 1; i < knots_count; ++i)
        {
            if (m_knots_vector[i] != knot)
                return i;
        }
        return INDEX_IS_OUTSIDE_OF_KNOTS_VECTOR;
    }

    int get_knots_multiplicity(int index) const
    {
        int konts_count = m_knots_vector.size();
        T knot = m_knots_vector[index];
        if (index >= konts_count || index < 0)
            return 0;
        int count = 1;
        int current_index = index - 1;
        while (current_index)
        {
            if (m_knots_vector[current_index] == knot)
            {
                ++count;
                --current_index;
            }
            else
                break;
        }
        current_index = index + 1;
        while (current_index < konts_count)
        {
            if (m_knots_vector[current_index] == knot)
            {
                ++count;
                ++current_index;
            }
            else
                break;
        }
        return count;
    }



    int is_include_konts(T u, T eps = 0.0)
    {
        int knots_size = m_knots_vector.size();
        for (int index = 0; index < knots_size; ++index)
        {
            if (std::abs(u - m_knots_vector[index]) <= eps)
                return index;
        }
        return INDEX_IS_OUTSIDE_OF_KNOTS_VECTOR;
    }

    T get_knot(int index) const
    {
        return m_knots_vector[index];
    }
    ENUM_NURBS set_degree(int degree) { m_degree = degree; return ENUM_NURBS::NURBS_SUCCESS; }

    ENUM_NURBS set_control_points(const Eigen::Matrix<T, rows, Eigen::Dynamic> &points) { m_control_points = points; return ENUM_NURBS::NURBS_SUCCESS; }

    ENUM_NURBS set_knots_vector(const Eigen::VectorX<T> &knots_vector) { m_knots_vector = knots_vector; return ENUM_NURBS::NURBS_SUCCESS; }

    Eigen::Matrix<T, rows, Eigen::Dynamic> get_control_points() const
    {
        return m_control_points;
    }

    Eigen::VectorX<T> get_knots_vector() const
    {
        return m_knots_vector;
    }

    int get_degree() const
    {
        return m_degree;
    }

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
        result.resize(n + 1);
        int index = -1;
        find_span<T>(u, m_degree, m_knots_vector, index);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ders(m_degree + 1, n + 1);
        int new_n = std::min(n, m_degree);
        ders_basis_funs<T>(index, new_n, m_degree, u, m_knots_vector, ders);
        Eigen::Vector<T, dim> zero;
        zero.setConstant(0.0);
        Eigen::VectorX<Eigen::Vector<T, rows>> temp(new_n + 1);

        for (int idx = 0; idx <= new_n; ++idx)
        {
            Eigen::Vector<T, rows> vec = m_control_points.block(0, index - m_degree, rows,m_degree + 1) * ders.col(idx);
            temp[idx] = vec;
        }
        Eigen::VectorX<Eigen::Vector<T, dim>> project_point = rat_curve_derivs_project<T, is_rational, rows>::project_point_to_euclidean_space(temp);
        for (int idx = 0; idx <= new_n; ++idx)
            result[idx] = project_point[idx];
        for (int idx = new_n + 1; idx <= n; ++idx)
        {
            result[idx] = zero;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
    /// @tparam n 求(0, 1, 2 ... n)右导数 < m_degree
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
        Eigen::Vector<Eigen::Vector<T, rows>, Eigen::Dynamic> temp(n + 1);
        // result.resize(n + 1);
        for (int idx = 0; idx <= n; ++idx)
        {
            Eigen::Vector<T, rows> vec = m_control_points.block(0, index - m_degree, rows, m_degree + 1) * ders.col(idx);
            temp[idx] = vec;
            // result[idx] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(vec);
        }
        result = rat_curve_derivs_project<T, is_rational, rows>::project_point_to_euclidean_space(temp);
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
        if (insert_knots.size() == 0)
            return ENUM_NURBS::NURBS_SUCCESS;
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
            // std::cout << new_control_points[index] << std::endl;
            bezier_curves[index] = cur;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 将nurbs曲线分解成bezier曲线
    /// @param bezier_curves 分解的bezier曲线, 内存用户释放
    /// @return ENUM_NURBS错误码
    ENUM_NURBS decompose_to_bezier(std::vector<nurbs_curve<T, dim, is_rational, -1, -1> *> &nurbs_curves)
    {
        // std::cout << m_knots_vector << std::endl;
        // std::cout << m_control_points << std::endl;
        Eigen::VectorX<T> new_knots_vector;
        Eigen::VectorX<Eigen::Matrix<T, rows, Eigen::Dynamic>> new_control_points;
        int interval_count = -1;
        find_interval_segment_count<T>(m_degree, m_knots_vector, interval_count);
        new_control_points.resize(interval_count);
        decompose_curve_to_bezier<T, rows>(m_degree, interval_count, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
        nurbs_curves.resize(interval_count, nullptr);
        T low = m_knots_vector[0];
        T high = m_knots_vector[m_degree + 1];
        int current_index = m_degree + 1;
        int new_knots_count = m_degree * 2 + 2;
        for (int index = 0; index < interval_count; ++index)
        {
            nurbs_curve<T, dim, is_rational, -1, -1> *cur = new nurbs_curve<T, dim, is_rational, -1, -1>();
            cur->set_control_points(new_control_points[index]);
            // std::cout << new_control_points[index] << std::endl;
            cur->set_degree(m_degree);
            Eigen::VectorX<T> knots_vector(new_knots_count);
            int half_count = new_knots_count / 2;
            for (int index = 0; index < half_count; ++index)
            {
                knots_vector[index] = low;
                knots_vector[new_knots_count - 1 - index] = high;
            }
            cur->set_knots_vector(knots_vector);
            nurbs_curves[index] = cur;
            low = high;
            current_index = get_next_different_index(current_index);
            //有空处理一下此处的逻辑
            if (current_index == INDEX_IS_OUTSIDE_OF_KNOTS_VECTOR)
                continue;
            high = m_knots_vector[current_index];
            
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

    /// @brief 将nurbs曲线降低一阶
    /// @param error 误差
    /// @return ENUM_NURBS错误码
    ENUM_NURBS degree_reduce(T error = DEFAULT_ERROR)
    {
        Eigen::VectorX<T> new_knots_vector;
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        ENUM_NURBS flag = degree_reduce_curve<T, rows, is_rational>(m_degree, m_knots_vector, m_control_points, new_knots_vector, new_control_points, error);
        if (flag != ENUM_NURBS::NURBS_SUCCESS)
            return flag;
        m_control_points = new_control_points;
        m_knots_vector = new_knots_vector;
        m_degree -= 1;
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    bool is_closed() const
    {
        int points_count = m_control_points.cols();
        Eigen::Vector<T, dim> start = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
        Eigen::Vector<T, dim> end = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(points_count - 1));
        return (start - end).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR;
    }

    //TODO: 写个多线程, 多线程每次结果可能都会不同
    ENUM_NURBS find_nearst_point_on_curve(const Eigen::Vector<T, rows> &point, T &u, Eigen::Vector<T, rows> &nearst_point, T eps = DEFAULT_ERROR) const
    {
        T min = m_knots_vector[0];
        int knots_size = m_knots_vector.size();
        T max = m_knots_vector(knots_size - 1);

        //将节点分成(max - min + 1) * 100份；以后需要优化
        int step_count = (max - min + 1) * 10;
        T step = (max - min) / step_count;
        T min_length = (point - m_control_points.col(0)).squaredNorm();
        nearst_point = nearst_point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
        u = m_knots_vector[0];
        Eigen::Vector<Eigen::Vector<T, dim>, 3> ders_vec;
        for (int index = 0; index <= step_count; ++index)
        {
            T distance;
            T current_u = index * step + min;
            for (int loop_index = 0; loop_index < MAX_ITERATE_DEEP; ++loop_index)
            {
                curve_derivs_alg2<2>(current_u, ders_vec);
                Eigen::Vector<T, dim> dist_vec = ders_vec[0] - point;
                distance = dist_vec.squaredNorm();
                if (distance < eps * eps)
                {
                    u = current_u;
                    nearst_point = ders_vec[0];
                    return ENUM_NURBS::NURBS_SUCCESS;
                }
                if (min_length > distance)
                {
                    min_length = distance;
                    u = current_u;
                    nearst_point = ders_vec[0];
                }

                T cos_angle = ders_vec[1].dot(dist_vec);
                T tanget_vec_square_len = ders_vec[1].squaredNorm();
                if (cos_angle * cos_angle < tanget_vec_square_len * distance * eps * eps)
                {
                    break;
                }
            
                T next_u = current_u - cos_angle / (ders_vec[2].dot(dist_vec) + tanget_vec_square_len);
                bool is_closed_flag = is_closed();
                if (is_closed_flag)
                {
                    if (next_u < min)
                        next_u = max - (min - next_u);
                    else if (next_u > max)
                        next_u = min + (next_u - max);
                }
                else
                {
                    if (next_u < min)
                        next_u = min;
                    else if (next_u > max)
                        next_u = max;
                }
                if ((next_u - current_u) * (next_u - current_u) * tanget_vec_square_len < eps * eps)
                    break;
                current_u = next_u;
            }
        
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    ENUM_NURBS parallel_projection(const Eigen::Vector<T, dim> &reference_point, const Eigen::Vector<T, dim> &normal, 
        const Eigen::Vector<T, dim> &projection_direction, nurbs_curve<T, dim, is_rational, -1, -1> &project_nurbs)
    {
        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        prallel_projection_curve<T, dim, is_rational>::parallel_projection(reference_point, normal, projection_direction, m_control_points, new_control_points);
        project_nurbs.set_control_points(new_control_points);
        project_nurbs.set_knots_vector(m_knots_vector);
        project_nurbs.set_degree(m_degree);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS perspective_projection(const Eigen::Vector<T, dim> &reference_point, const Eigen::Vector<T, dim> &normal, 
        const Eigen::Vector<T, dim> &eye, nurbs_curve<T, dim, true, -1, -1> &project_nurbs)
    {
        Eigen::Matrix<T, dim + 1, Eigen::Dynamic> new_control_points;
        perspective_projection_curve<T, dim, is_rational>::perspective_projection(reference_point, normal, eye, m_control_points, new_control_points);
        project_nurbs.set_control_points(new_control_points);
        project_nurbs.set_knots_vector(m_knots_vector);
        project_nurbs.set_degree(m_degree);
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    // //bezier_curve
    // // template<bool is_rational_function>
    // ENUM_NURBS bezier_curve_reparameter(const nurbs_curve<T, 1, false, -1, -1> &reparameter_function,
    //     nurbs_curve<T, dim, is_rational, -1, -1> &new_nurbs)
    // {
    //     int function_degree = reparameter_function.get_degree();
    //     T low = reparameter_function.get_knot(0);
    //     T high = reparameter_function.get_knot(function_degree + 1);
     
    //     int new_degree = m_degree * function_degree;
    //     // Eigen::VectorX<T> new_nurbs_vector(2 * (new_degree + 1));
    //     Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
    //     new_control_points.resize(rows, new_degree + 1);
    //     new_control_points.col(0) = m_control_points.col(0);
    //     new_control_points.col(new_degree) = m_control_points.col(m_degree);

    //     int ml = (new_degree + 1) / 2;
    //     int mr = new_degree - new_degree / 2 - 1;
    //     int pqf =  std::tgamma<int>(new_degree + 1);
        
    //     Eigen::MatrixX<int> Bin = binary_coeff(ml + 1);

    //     Eigen::VectorX<Eigen::Vector<T, rows>> ders(ml + 1);
    //     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ders_temp;
    //     int new_n = std::min(ml, m_degree);
    //     ders_basis_funs<T>(m_degree, new_n, m_degree, m_knots_vector[0], m_knots_vector, ders_temp);
    //     Eigen::Vector<T, rows> zero;
    //     zero.setConstant(0.0);
    //     for (int idx = 0; idx <= new_n; ++idx)
    //     {
    //         Eigen::Vector<T, rows> vec = m_control_points.block(0, 0, rows, m_degree + 1) * ders_temp.col(idx);
    //         ders[idx] = vec;
    //     }
    //     for (int idx = new_n + 1; idx <= ml; ++idx)
    //     {
    //         ders[idx] = zero;
    //     }

    //     Eigen::VectorX<Eigen::Vector<T, 1>> function_ders;
    //     reparameter_function.derivative_on_curve(low, ml, function_ders);
  
    //     Eigen::VectorX<Eigen::Vector<T, rows>> ders_s(ml + 1);
    //     for (int index = 0; index <= ml; ++index)
    //     {
    //         ders_s[index].setConstant(0.0);
    //     }
    //     ders_s[0] = ders[0];
        
    //     for (int n = 1; n <= ml; ++n)
    //     {
    //         for (int j = 0; j <= n; ++j)
    //         {
    //             std::vector<std::vector<int>> temp = find_sum_n(n, j, 0, false);
    //             for (auto vec : temp)
    //             {
    //                 int vec_sum = 0;
    //                 for (int index = 0; index < n; ++index)
    //                     vec_sum += (index + 1) * vec[index];
    //                 if (vec_sum == n)
    //                 {
    //                     T coeff = std::tgamma<int>(n + 1) * std::pow(function_ders[1][0], vec[0]);
    //                     int denominator = std::tgamma<int>(vec[0] + 1);
    //                     for (int index = 1; index < n; ++index)
    //                     {
    //                         denominator *= std::tgamma<int>(vec[index] + 1) * std::pow(std::tgamma<int>(index + 2), vec[index]);
    //                         coeff *= std::pow(function_ders[index + 1][0], vec[index]);
    //                     }
    //                     ders_s[n] += ders[j] * (coeff / denominator);
    //                 }

    //             }   
    //         }
    //     }


    //     for (int index = 1; index <= ml; ++index)
    //     {
    //         int num = std::tgamma<int>(new_degree + 1 - index);
    //         T s = std::pow(high - low, index);
    //         new_control_points.col(index) = ((num * s) / pqf) * ders_s[index];
    //         for (int j = 0; j < index; ++j)
    //         {
    //             int coeff = (index + j - 1) % 2 == 0 ? 1 : -1;
    //             new_control_points.col(index) = new_control_points.col(index) + coeff * Bin(index, j) * new_control_points.col(j);
    //         }function_bezier_curves
    //     }

        
        
    //     Bin = binary_coeff(mr + 1);
    //     new_n = std::min(mr, m_degree);
    //     ders_basis_funs<T>(m_degree, new_n, m_degree, m_knots_vector[m_degree + 1], m_knots_vector, ders_temp);
    //     ders.resize(mr + 1);
    //     int knots_index = m_knots_vector.size() - m_degree * 2 - 2;
    //     for (int idx = 0; idx <= new_n; ++idx)
    //     {
    //         Eigen::Vector<T, rows> vec = m_control_points.block(0, knots_index, rows, m_degree + 1) * ders_temp.col(idx);
    //         ders[idx] = vec;
    //     }
    //     for (int idx = new_n + 1; idx <= mr; ++idx)
    //     {
    //         ders[idx] = zero;
    //     }

    //     ders_s.resize(mr + 1);
    //     ders_s[0] = ders[0];
    //     for (int index = 1; index <= mr; ++index)
    //         ders_s[index] = zero;
        

    //     reparameter_function.derivative_on_curve(high, mr, function_ders);
    //     for (int n = 1; n <= mr; ++n)
    //     {
    //         for (int j = 0; j <= n; ++j)
    //         {
    //             std::vector<std::vector<int>> temp = find_sum_n(n, j, 0, false);
    //             for (auto vec : temp)
    //             {
    //                 int vec_sum = 0;
    //                 for (int index = 0; index < mr; ++index)
    //                     vec_sum += (index + 1) * vec[index];
    //                 if (vec_sum == n)
    //                 {
    //                     T coeff = std::tgamma<int>(n + 1) * std::pow(function_ders[1][0], vec[0]);
    //                     int denominator = std::tgamma<int>(vec[0] + 1);
    //                     for (int index = 1; index < n; ++index)
    //                     {
    //                         denominator *= std::tgamma<int>(vec[index] + 1) * std::pow(std::tgamma<int>(index + 2), vec[index]);
    //                         coeff *= std::pow(function_ders[index + 1][0], vec[index]);
    //                     }
    //                     ders_s[n] += ders[j] * (coeff / denominator);
    //                 }

    //             }   
    //         }
    //     }

        
    //     for (int index = 1; index <= mr; ++index)
    //     {
    //         int num = std::tgamma<int>(new_degree + 1 - index);
    //         T s = std::pow(high - low, index);
    //         int coeff = index % 2 == 0 ? 1 : -1;
    //         new_control_points.col(new_degree - index) = ((coeff * num * s) / pqf) * ders_s[index];
    //         for (int j = 0; j < index; ++j)
    //         {
    //             int coeff1 = (index + j - 1) % 2 == 0 ? 1 : -1;
    //             new_control_points.col(new_degree - index) = new_control_points.col(new_degree - index) + coeff1 * Bin(index, j) * new_control_points.col(new_degree - j);
    //         }
    //     }
    //     new_nurbs.set_control_points(new_control_points);
    //     Eigen::VectorX<T> new_knots_vector(2 * new_degree + 2);
    //     for (int index = 0; index <= new_degree; ++index)
    //     {
    //         new_knots_vector[index] = low;
    //         new_knots_vector[index + new_degree + 1] = high;
    //     }
    //     new_nurbs.set_knots_vector(new_knots_vector);
    //     new_nurbs.set_degree(new_degree);
    //     return ENUM_NURBS::NURBS_SUCCESS;
    // }
// private:

    ENUM_NURBS bezier_curve_reparameter(const nurbs_curve<T, 1, false, -1, -1> &reparameter_function,
        nurbs_curve<T, dim, is_rational, -1, -1> &new_nurbs)
    {
        std::array<Eigen::Vector<T, 1>, 2> ends_points;
        int knots_vector_size = m_knots_vector.size();
        reparameter_function.get_ends_point(ends_points);

        if (std::abs(ends_points[0][0] - m_knots_vector[0]) > DEFAULT_ERROR)
            return ENUM_NURBS::NURBS_ERROR;
        if (std::abs(ends_points[1][0] - m_knots_vector[knots_vector_size - 1]) > DEFAULT_ERROR)
           return ENUM_NURBS::NURBS_ERROR;

    
        int function_degree = reparameter_function.get_degree();
        T low = reparameter_function.get_knot(0);
        T high = reparameter_function.get_knot(function_degree + 1);     
        int new_degree = m_degree * function_degree;
        int ml = (new_degree + 1) / 2;
        int mr = new_degree - new_degree / 2 - 1;

        Eigen::VectorX<Eigen::Vector<T, 1>> function_low_ders;
        reparameter_function.derivative_on_curve(low, ml, function_low_ders);

        Eigen::VectorX<Eigen::Vector<T, 1>> function_high_ders;
        reparameter_function.derivative_on_curve(high, mr, function_high_ders);

        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        new_control_points.resize(rows, new_degree + 1);

        new_control_points.col(0) = m_control_points.col(0);
        new_control_points.col(new_degree) = m_control_points.col(m_degree);
        
        reparameter_bezier_curve(m_degree, m_knots_vector, m_control_points, function_degree, 
            function_low_ders, function_high_ders, high - low, new_control_points);
        // new_control_points.col(0) = m_control_points.nurbs_curve<T, dim, is_rational, -1, -1>
        
        new_nurbs.set_control_points(new_control_points);
        Eigen::VectorX<T> new_knots_vector(2 * new_degree + 2);
        for (int index = 0; index <= new_degree; ++index)
        {
            new_knots_vector[index] = low;
            new_knots_vector[index + new_degree + 1] = high;
        }
        new_nurbs.set_knots_vector(new_knots_vector);
        new_nurbs.set_degree(new_degree);
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    ENUM_NURBS curve_reparameter_with_polynomial(const nurbs_curve<T, 1, false, -1, -1> &reparameter_function,
        nurbs_curve<T, dim, is_rational, -1, -1> &new_nurbs)
    {
        std::array<Eigen::Vector<T, 1>, 2> ends_points;
        int knots_vector_size = m_knots_vector.size();
        reparameter_function.get_ends_point(ends_points);

        if (std::abs(ends_points[0][0] - m_knots_vector[0]) > DEFAULT_ERROR)
            return ENUM_NURBS::NURBS_ERROR;
        if (std::abs(ends_points[1][0] - m_knots_vector[knots_vector_size - 1]) > DEFAULT_ERROR)
           return ENUM_NURBS::NURBS_ERROR;

        int this_nurbs_interval_count;
        find_interval_segment_count(m_degree, m_knots_vector, this_nurbs_interval_count);
        Eigen::VectorX<T> function_knots_vector = reparameter_function.get_knots_vector();
        int reparameter_function_interval_count;
        int reparameter_function_degree = reparameter_function.get_degree();
        find_interval_segment_count(reparameter_function_degree, function_knots_vector, reparameter_function_interval_count);
        if (this_nurbs_interval_count == 1 && reparameter_function_interval_count == 1)
            return bezier_curve_reparameter(reparameter_function, new_nurbs);

        
        T current_knot = function_knots_vector[0];
        std::vector<T> insert_knots;
        std::vector<int> insert_knots_degree;
        std::vector<int> equal_knots;
        std::vector<int> equal_knots_degree;
        

        int new_degree = reparameter_function_degree * m_degree;
        int function_knots_vector_size = function_knots_vector.size();
        
        int reparameter_function_control_points_count = function_knots_vector_size - reparameter_function_degree - 1;
        for (int index = reparameter_function_degree + 1; index < reparameter_function_control_points_count; ++index)
        {
            if (current_knot == function_knots_vector[index])
                continue;
            current_knot = function_knots_vector[index];
            Eigen::Vector<T, 1> image_knot;
            reparameter_function.point_on_curve(current_knot, image_knot);
            int knots_index = is_include_konts(image_knot[0], DEFAULT_ERROR);
            if (knots_index == INDEX_IS_OUTSIDE_OF_KNOTS_VECTOR)
            {
                insert_knots.push_back(image_knot[0]);
                int repeat = reparameter_function.get_knots_multiplicity(index);
                insert_knots_degree.push_back(reparameter_function_degree - repeat - 1);
            }
            else
            {
                equal_knots.push_back(index);

                int repeat1 = reparameter_function.get_knots_multiplicity(index);
                int repeat2 = get_knots_multiplicity(knots_index);
                int repeat = std::min(reparameter_function_degree - repeat1 - 1, m_degree - repeat2 - 1);
                equal_knots_degree.push_back(repeat);
            }
        }

        nurbs_curve<T, dim, is_rational, -1, -1> *insert_knots_nurbs = new nurbs_curve<T, dim, is_rational, -1, -1>(*this);
        std::vector<nurbs_curve<T, dim, is_rational, -1, -1> *> bezier_curves;
        int insert_knots_size = insert_knots.size();
        Eigen::VectorX<T> temp_insert_knots(insert_knots_size);
        for (int index = 0; index < insert_knots_size; ++index)
            temp_insert_knots[index] = insert_knots[index];
        
        insert_knots_nurbs->refine_knots_vector(temp_insert_knots);
        insert_knots_nurbs->decompose_to_bezier(bezier_curves);
            
        delete insert_knots_nurbs;


        std::vector<T> reparameter_function_insert_knots;
        std::vector<int> reparameter_function_insert_knots_degree;

        current_knot = m_knots_vector[0];
        auto next_equal_knots = equal_knots.begin();
        auto eqaul_knots_end = equal_knots.end();
        int control_points_count = knots_vector_size - m_degree - 1;
        for (int index = m_degree + 1; index < control_points_count; ++index)
        {
            if (next_equal_knots != eqaul_knots_end)
            {
                if (index == *next_equal_knots)
                {
                    ++next_equal_knots;
                    current_knot = m_knots_vector[index];
                    continue;
                }
            }

            if (current_knot == m_knots_vector[index])
                continue;
            current_knot = m_knots_vector[index];

            Eigen::Vector<T, 1> point{current_knot}, nearst_point;
            T u;
            reparameter_function.find_nearst_point_on_curve(point, u, nearst_point, KNOTS_VECTOR_ERROR);
            T distacne = (point - nearst_point).squaredNorm();
            if (distacne > KNOTS_VECTOR_ERROR * KNOTS_VECTOR_ERROR)
            {
                int curve_count = bezier_curves.size();
                for (int curve_index = 0; curve_index < curve_count; ++curve_index)
                    delete bezier_curves[index];
                return ENUM_NURBS::NURBS_ERROR;
            }
            reparameter_function_insert_knots.push_back(u);
            int repeat = get_knots_multiplicity(index);
            reparameter_function_insert_knots_degree.push_back(m_degree - repeat - 1);
        }

        nurbs_curve<T, 1, false, -1, -1> *insert_knots_function_nurbs = new nurbs_curve<T, 1, false, -1, -1>(reparameter_function);
        std::vector<nurbs_curve<T, 1, false, -1, -1> *> function_bezier_curves;
        int reparameter_function_insert_knots_size = reparameter_function_insert_knots.size();
        Eigen::VectorX<T> temp_insert_knots_vector(reparameter_function_insert_knots_size);
        for (int index = 0; index < reparameter_function_insert_knots_size; ++index)
            temp_insert_knots_vector[index] = reparameter_function_insert_knots[index];
        insert_knots_function_nurbs->refine_knots_vector(temp_insert_knots_vector);
        insert_knots_function_nurbs->decompose_to_bezier(function_bezier_curves);
        delete insert_knots_function_nurbs;

        int curves_count = function_bezier_curves.size();
        std::vector<nurbs_curve<T, dim, is_rational, -1, -1>*> new_bezier_curves(curves_count, nullptr);
        for (int index = 0; index < curves_count; ++index)
        {
            new_bezier_curves[index] = new nurbs_curve<T, dim, is_rational, -1, -1>();
            bezier_curves[index]->bezier_curve_reparameter(*function_bezier_curves[index], *(new_bezier_curves[index]));
            
            //debug
            // int a1 = bezier_curves[index]->get_knots_count();
            // T low = bezier_curves[index]->get_knot(0);
            // T high = bezier_curves[index]->get_knot(a1 - 1);
            // T step = (high - low) / 100;
            // std::vector<Eigen::Vector<T, 2>> pointss;
            // std::string dir("view" + std::to_string(index) + ".obj");
            // std::ofstream outfile(dir);
            // for (int i = 0; i < 100; ++i)
            // {
            //     T u = low + i * step;

            //     Eigen::Vector2d point;
            //     bezier_curves[index]->point_on_curve(u, point);
            //     pointss.push_back(point);
            // }
            // for (auto point : pointss)
            // {
            //     outfile << "v " << point[0] << " " <<
            //     point[1] << " "<< 0.0 << std::endl;
            // }
            // outfile.close();

            delete bezier_curves[index];
            delete function_bezier_curves[index];
        }

        Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
        Eigen::VectorX<T> new_knots_vector;

        merge_two_curve(new_degree, new_degree, new_bezier_curves[0]->get_knots_vector(), new_bezier_curves[1]->get_knots_vector(),
           new_bezier_curves[0]->get_control_points(), new_bezier_curves[1]->get_control_points(),
           new_knots_vector, new_control_points);
        
        for (int index = 2; index < curves_count; ++index)
        {
            Eigen::Matrix<T, rows, Eigen::Dynamic> nnew_control_points;
            Eigen::VectorX<T> nnew_knots_vector;
            merge_two_curve(new_degree, new_degree, new_knots_vector, new_bezier_curves[index]->get_knots_vector(),
                new_control_points, new_bezier_curves[index]->get_control_points(),
                nnew_knots_vector, nnew_control_points);
            new_knots_vector = nnew_knots_vector;
            new_control_points = nnew_control_points;
        }
        new_nurbs.set_control_points(new_control_points);
        new_nurbs.set_degree(new_degree);
        new_nurbs.set_knots_vector(new_knots_vector);

        int remove_count = insert_knots.size();
        for (int index = 0; index < remove_count; ++index)
        {
            int time;
            new_nurbs.remove_knots(insert_knots[index], insert_knots_degree[index], time);
            if (time != insert_knots_degree[index])
                return ENUM_NURBS::NURBS_ERROR;
        }

        remove_count = equal_knots.size();
        for (int index = 0; index < remove_count; ++index)
        {
            int time;
            new_nurbs.remove_knots(function_knots_vector[equal_knots[index]], equal_knots_degree[index], time);
            if (time != equal_knots_degree[index])
                return ENUM_NURBS::NURBS_ERROR;
        }

        remove_count = reparameter_function_insert_knots.size();
        for (int index = 0; index < remove_count; ++index)
        {
            int time;
            new_nurbs.remove_knots(reparameter_function_insert_knots[index], reparameter_function_insert_knots_degree[index], time);
            if (time != reparameter_function_insert_knots_degree[index])
                return ENUM_NURBS::NURBS_ERROR;
        }

        return ENUM_NURBS::NURBS_SUCCESS;
    }



};