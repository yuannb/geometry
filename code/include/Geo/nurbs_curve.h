#pragma once
#include "nurbs_tool.h"
#include "bezier_curve.h"
#include <array>
#include "curve.h"
#include <concepts>
// #include "ThreadPool.h"
//TODO:将knots_vector的类型将Eigen::vector换成std::vector
namespace tnurbs
{
    using namespace tnurbs;
    /// @brief nurbs curve class
    /// @tparam T : double float int ...
    /// @tparam dim : nurbs所在欧氏空间的维数
    /// @tparam is_rational : 是否是有理的
    /// @tparam points_count : 控制点个数
    /// @tparam degree : nurbs的阶数
    template<typename T, int dim, bool is_rational, int points_count, int degree>
    class nurbs_curve : public curve<nurbs_curve<T, dim, is_rational, points_count, degree>>
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

        
        /// @brief 给定误差, 在误差范围内消去所有可以消去的节点; 此函数仅在非有理nurbs曲线是合法的(要求参数从小到大排列, 且不重复, 稍加改造应该可以使得消去误差在整个参数域都成立)
        /// @param params 消去的曲线在参数params上误差小于给定的误差; 即(new_curve(params[i]) - curve(params[i])).norm() < E;
        /// @param error 返回的各个点的误差(如果params为空, 则返回各个节点区间的误差)
        /// @param new_nurbs 节点消去后新的nurbs
        /// @param E 最大误差
        /// @return 错误码
        ENUM_NURBS remove_knots_bound_curve(const std::vector<T> &params, std::vector<T> &error, nurbs_curve<T, dim, false, -1, -1> &new_nurbs, T E = DEFAULT_ERROR) const;
        int get_control_points_count() const;
        ENUM_NURBS get_control_point(int index, Eigen::Vector<T, dim> &point) const;
        Eigen::VectorX<T> get_knots_vector() const;
        ENUM_NURBS degree_elevate(int t);
        ENUM_NURBS curve_reparameter_with_linear_function(T alpha, T beta,
            nurbs_curve<T, dim, is_rational, -1, -1> &new_nurbs) const;
        Eigen::Matrix<T, rows, Eigen::Dynamic> get_control_points() const;
        ENUM_NURBS get_weight(int index, T &w) const;
        ENUM_NURBS set_nonhome_control_points(const Eigen::Matrix<T, dim, points_count> &control_points);
        
        
        
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
            Eigen::Matrix<T, degree + 1, Eigen::Dynamic> ders(degree + 1, n + 1);
            ders_basis_funs<T, degree, points_count>(index, n, u, m_knots_vector, ders);

            Eigen::VectorX<Eigen::Vector<T, rows>> temp(n + 1);
            for (int idx = 0; idx <= n; ++idx)
            {
                Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, degree + 1) * ders.col(idx);
                temp[idx] = vec;
            }
            result = rat_curve_derivs_project<T, is_rational, rows>::project_point_to_euclidean_space(temp);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
        /// @tparam n 求(0, 1, 2 ... n)右导数
        /// @param u 曲线参数
        /// @param result out_put_pram, result[i]为nurbs曲线的第n次导数
        /// @return ENUM_NURBS错误码
        template<int n, ENUM_LIMITDIRECTION flag = ENUM_LIMITDIRECTION::RIGHT>
        ENUM_NURBS derivative_on_curve(T u, Eigen::Vector<Eigen::Vector<T, dim>, n + 1> &result) const
        {
            int index = -1;
            find_span<T, points_count, degree, flag>(u, m_knots_vector, index);
            Eigen::Matrix<T, degree + 1, n + 1> ders(degree + 1, n + 1);
            ders_basis_funs<T, degree, points_count, n>(index, u, m_knots_vector, ders);

            Eigen::Vector<Eigen::Vector<T, rows>, n + 1> temp;
            for (int idx = 0; idx <= n; ++idx)
            {
                Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, degree + 1) * ders.col(idx);
                temp[idx] = vec;
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
            decompose_curve_to_bezier<T, rows, points_count, degree>(interval_count, m_knots_vector, m_control_points, new_knots_vector, new_control_points);
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
            nearst_point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
            T min_length = (point - nearst_point).squaredNorm();
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
    class nurbs_curve<T, dim, is_rational, -1, degree> : public curve<nurbs_curve<T, dim, is_rational, -1, degree>>
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
            Eigen::Matrix<T, degree + 1, Eigen::Dynamic> ders(degree + 1, n + 1);
            ders_basis_funs<T, degree>(index, n, u, m_knots_vector, ders);

            Eigen::VectorX<Eigen::Vector<T, rows>> temp(n + 1);
            for (int idx = 0; idx <= n; ++idx)
            {
                Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, degree + 1) * ders.col(idx);
                temp[idx] = vec;
            }
            result = rat_curve_derivs_project<T, is_rational, rows>::project_point_to_euclidean_space(temp);
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
            find_span<T, degree>(u, m_knots_vector, index);
            Eigen::Matrix<T, degree + 1, n + 1> ders(degree + 1, n + 1);
            ders_basis_funs<T, degree, n>(index, u, m_knots_vector, ders);
            Eigen::Vector<Eigen::Vector<T, rows>, n + 1> temp;
            for (int idx = 0; idx <= n; ++idx)
            {
                Eigen::Vector<T, rows> vec = m_control_points.block(0, index - degree, rows, degree + 1) * ders.col(idx);
                temp[idx] = vec;
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
        ENUM_NURBS find_nearst_point_on_curve(const Eigen::Vector<T, dim> &point, T &u, Eigen::Vector<T, dim> &nearst_point, T eps = DEFAULT_ERROR) const
        {
            T min = m_knots_vector[0];
            int knots_size = m_knots_vector.size();
            T max = m_knots_vector(knots_size - 1);

            //将节点分成(max - min + 1) * 100份；以后需要优化
            int step_count = (max - min + 1) * 10;
            T step = (max - min) / step_count;
            nearst_point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
            T min_length = (point - m_control_points.col(0)).squaredNorm();
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
    class nurbs_curve<T, dim, is_rational, -1, -1> : public curve<nurbs_curve<T, dim, is_rational, -1, -1>>
    {

    private:
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
                        m_degree = m_knots_vector.size() - m_control_points.cols() - 1; 
                    }
        nurbs_curve(const nurbs_curve<T, dim, is_rational, -1, -1> &nurbs_curve_to_copy)
        {
            m_knots_vector = nurbs_curve_to_copy.get_knots_vector();
            m_degree = nurbs_curve_to_copy.get_degree();
            m_control_points = nurbs_curve_to_copy.get_control_points();
        }

        // //生成值全为1的一阶的nurbs曲线
        // static nurbs_curve<T, dim, false, -1, -1> *make_identity_nurbs(T low, T high) const
        // {
        //     Eigen::VectorX<T> knots(2);
        //     knots << low, high;
        //     Eigen::Matrix<T, dim, Eigen::Dynamic> control_points(dim, 1);
        //     control_points.col(0) = Eigen::Vector<T, dim>(1.0, 1.0, 1.0);
        //     nurbs_curve<T, dim, false, -1, -1> *id_nurbs = new nurbs_curve<T, dim, false, -1, -1>(knots, control_points);
        //     return id_nurbs;
        // }

        ENUM_NURBS get_ends_knots(std::array<T, 2> &ends_knots) const
        {
            int konts_size = m_knots_vector.size();
            ends_knots[0] = m_knots_vector[0];
            ends_knots[1] = m_knots_vector[konts_size - 1];
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        // /// @brief 将节点矢量去重, 并且记录各个节点的重复度
        // /// @param different_knots 去重后的节点矢量
        // /// @param multiples 各个节点的重复度
        // /// @return 错误码
        // ENUM_NURBS get_different_knots(std::vector<T> &different_knots, std::vector<int> &multiples) const
        // {
        //     int knots_count = m_knots_vector.size();
        //     different_knots.reserve(knots_count - 2 * m_degree);
        //     multiples.reserve(knots_count - 2 * m_degree);
        //     int current_knots_multiple = 0;
        //     T current_knots = m_knots_vector[0];
        //     for (int index = 0; index < knots_count; ++index)
        //     {
        //         if (current_knots != m_knots_vector[index])
        //         {
        //             different_knots.push_back(current_knots);
        //             multiples.push_back(current_knots_multiple);
        //             current_knots = m_knots_vector[index];
        //             current_knots_multiple = 1;
        //         }
        //         else
        //             current_knots_multiple += 1;
        //     }
        //     different_knots.push_back(current_knots);
        //     multiples.push_back(current_knots_multiple);
        //     return ENUM_NURBS::NURBS_SUCCESS;
        // }


        ENUM_NURBS get_ends_point(std::array<Eigen::Vector<T, dim>, 2> &points) const
        {
            points[0] = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
            int points_count = m_control_points.cols();
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

        ENUM_NURBS get_weight(int index, T &w) const
        {
            if (index < 0 || index >= m_control_points.cols())
                return ENUM_NURBS::NURBS_ERROR;
            if constexpr (is_rational == true)
            {
                w = m_control_points(dim, index);
                return  ENUM_NURBS::NURBS_SUCCESS;
            }
            //else
            w = 1.0;
            return ENUM_NURBS::NURBS_SUCCESS;

        }

        ENUM_NURBS get_control_point(int index, Eigen::Vector<T, dim> &point) const
        {
            if (index < 0 || index >= m_control_points.cols())
                return ENUM_NURBS::NURBS_ERROR;
            if constexpr (is_rational == true)
            {
                point = m_control_points.template block<dim, 1>(0, index) / m_control_points(dim, index);
            }
            else
            {
                point = m_control_points.col(index);
            }

            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS get_homo_control_point(int index, Eigen::Vector<T, dim + 1> &point) const
        {
            if (index < 0 || index >= m_control_points.cols())
                return ENUM_NURBS::NURBS_ERROR;
            if constexpr (is_rational == true)
            {
                point = m_control_points.col(index);
            }
            else
            {
                point.block(0, 0, dim, 1) = m_control_points.col(index);
                point[dim] = 1.0;
            }

            return ENUM_NURBS::NURBS_SUCCESS;
        }


        int is_include_konts(T u, T eps = 0.0) const
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


        ENUM_NURBS set_control_points(const Eigen::Matrix<T, dim, Eigen::Dynamic> &points, const Eigen::VectorX<T> &weights)
        {
            static_assert(is_rational == true, "cannot set weight of non-rational b-nurbs");
            int points_count = points.cols();
            int weights_count = weights.rows();
            if (points_count != weights_count)
                return ENUM_NURBS::NURBS_ERROR;
            m_control_points.resize(rows, points_count);
            for (int index = 0; index < points_count; ++index)
            {
                m_control_points.block(0, index, dim, 1) = points.col(index) * weights[index];
                m_control_points(dim, index) = weights[index];
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS set_control_points(const Eigen::Matrix<T, rows, Eigen::Dynamic> &points) { m_control_points = points; return ENUM_NURBS::NURBS_SUCCESS; }


        ENUM_NURBS set_nonhome_control_points(const Eigen::Matrix<T, dim, Eigen::Dynamic> &control_points)
        {
            int cols = control_points.cols();
            m_control_points.resize(rows, cols);
            m_control_points.block(dim, cols, 0, 0) = control_points;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS set_knots_vector(const Eigen::VectorX<T> &knots_vector) 
        {
            m_knots_vector = knots_vector;
            return ENUM_NURBS::NURBS_SUCCESS; 
        }

        ENUM_NURBS to_rational_nurbs(nurbs_curve<T, dim, true, -1, -1> &new_nurbs)
        {
            Eigen::Matrix<T, dim + 1, Eigen::Dynamic> new_control_points;
            to_ratioanl_contrl_points<T, is_rational, rows>::convert(m_control_points, new_control_points);
            new_nurbs.set_control_points(new_control_points);
            new_nurbs.set_knots_vector(m_knots_vector);
            new_nurbs.set_degree(m_degree);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        Eigen::Matrix<T, rows, Eigen::Dynamic> get_control_points() const
        {
            return m_control_points;
        }

        Eigen::Matrix<T, dim, Eigen::Dynamic> get_nonhomo_control_points() const
        {
            if constexpr (is_rational == false)
                return m_control_points;
            else
            {
                int points_count = m_control_points.cols();
                Eigen::Matrix<T, dim, Eigen::Dynamic> cp(dim, points_count);
                for (int index = 0; index < points_count; ++index)
                {
                    cp.col(index) = m_control_points.block(0, index, dim, 1) / m_control_points(index, dim);
                }
                return std::move(cp);
            }
        }

        int get_control_points_count() const
        {
            return m_control_points.cols();
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

        /// @brief 计算nurbs曲线上在参数u处的权重
        /// @param u 曲线参数
        /// @param w out_put_param nurbs曲线参数u对应的权重
        /// @return ENUM_NURBS错误码
        ENUM_NURBS weight_on_curve(T u, T &w) const
        {
            if constexpr (is_rational == false)
            {
                w = 1.0;
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            else
            {
                int index = -1;
                find_span<T>(u,m_degree, m_knots_vector, index);
                Eigen::Vector<T, Eigen::Dynamic> coeff(m_degree + 1);
                basis_functions<T>(index, u, m_degree, m_knots_vector, coeff);
                Eigen::Vector<T, rows> vec = m_control_points.block(0, index - m_degree, rows, m_degree + 1) * coeff;
                w = vec[dim];
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            return ENUM_NURBS::NURBS_ERROR;
        }


        /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
        /// @param u 曲线参数
        /// @param n 求(0, 1, 2 ... n)右导数
        /// @param result out_put_pram, result[i]为nurbs曲线的第n次导数
        /// @return ENUM_NURBS错误码
        template<ENUM_LIMITDIRECTION flag = ENUM_LIMITDIRECTION::RIGHT>
        ENUM_NURBS derivative_on_curve(T u, int n, Eigen::Vector<Eigen::Vector<T, dim>, Eigen::Dynamic> &result) const
        {
            int index = -1;
            find_span<T, flag>(u, m_degree, m_knots_vector, index);
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ders(m_degree + 1, n + 1);
            ders_basis_funs<T>(index, n, m_degree, u, m_knots_vector, ders);
            Eigen::VectorX<Eigen::Vector<T, rows>> temp(n + 1);

            for (int idx = 0; idx <= n; ++idx)
            {
                Eigen::Vector<T, rows> vec = m_control_points.block(0, index - m_degree, rows,m_degree + 1) * ders.col(idx);
                temp[idx] = vec;
            }
            result = rat_curve_derivs_project<T, is_rational, rows>::project_point_to_euclidean_space(temp);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 计算nurbs曲线的右导数(结果向量不单位化)
        /// @tparam n 求(0, 1, 2 ... n)右导数 < m_degree
        /// @param u 曲线参数
        /// @param result out_put_pram, result[i]为nurbs曲线的第n次导数
        /// @return ENUM_NURBS错误码
        template<int n, ENUM_LIMITDIRECTION flag = ENUM_LIMITDIRECTION::RIGHT>
        ENUM_NURBS derivative_on_curve(T u, Eigen::Vector<Eigen::Vector<T, dim>, n + 1> &result) const
        {
            int index = -1;
            find_span<T, flag>(u, m_degree, m_knots_vector, index);
            Eigen::Matrix<T, Eigen::Dynamic, n + 1> ders(m_degree + 1, n + 1);
            ders_basis_funs<T, n>(index, m_degree, u, m_knots_vector, ders);
            Eigen::Vector<Eigen::Vector<T, rows>, n + 1> temp;
            for (int idx = 0; idx <= n; ++idx)
            {
                Eigen::Vector<T, rows> vec = m_control_points.block(0, index - m_degree, rows, m_degree + 1) * ders.col(idx);
                temp[idx] = vec;
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
        /// @param insert_knots 插入节点数组(要求升序)
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
        ENUM_NURBS decompose_to_bezier(std::vector<bezier_curve<T, dim, is_rational, -1> *> &bezier_curves) const
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

        /// @brief 将nurbs曲线分解成bezier曲线
        /// @param bezier_curves 分解的bezier曲线, 内存用户释放
        /// @return ENUM_NURBS错误码
        ENUM_NURBS decompose_to_bezier(std::vector<nurbs_curve<T, dim, is_rational, -1, -1> *> &nurbs_curves) const
        {
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
        ENUM_NURBS find_nearst_point_on_curve(const Eigen::Vector<T, dim> &point, T &u, Eigen::Vector<T, dim> &nearst_point, T eps = DEFAULT_ERROR) const
        {
            T min = m_knots_vector[0];
            int knots_size = m_knots_vector.size();
            T max = m_knots_vector(knots_size - 1);

            //将节点分成(max - min + 1) * 100份；以后需要优化
            int step_count = (max - min + 1) * 10;
            T step = (max - min) / step_count;
            nearst_point = project_point<T, is_rational, rows>::project_point_to_euclidean_space(m_control_points.col(0));
            T min_length = (point - nearst_point).squaredNorm();
            u = m_knots_vector[0];
            Eigen::Vector<Eigen::Vector<T, dim>, 3> ders_vec;
            for (int index = 0; index <= step_count; ++index)
            {
                T distance;
                T current_u = index * step + min;
                for (int loop_index = 0; loop_index < MAX_ITERATE_DEEP; ++loop_index)
                {
                    derivative_on_curve<2>(current_u, ders_vec);
                    // curve_derivs_alg2<2>(current_u, ders_vec);
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
            const Eigen::Vector<T, dim> &projection_direction, nurbs_curve<T, dim, is_rational, -1, -1> &project_nurbs) const
        {
            Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
            prallel_projection_curve<T, dim, is_rational>::parallel_projection(reference_point, normal, projection_direction, m_control_points, new_control_points);
            project_nurbs.set_control_points(new_control_points);
            project_nurbs.set_knots_vector(m_knots_vector);
            project_nurbs.set_degree(m_degree);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS perspective_projection(const Eigen::Vector<T, dim> &reference_point, const Eigen::Vector<T, dim> &normal, 
            const Eigen::Vector<T, dim> &eye, nurbs_curve<T, dim, true, -1, -1> &project_nurbs) const
        {
            Eigen::Matrix<T, dim + 1, Eigen::Dynamic> new_control_points;
            perspective_projection_curve<T, dim, is_rational>::perspective_projection(reference_point, normal, eye, m_control_points, new_control_points);
            project_nurbs.set_control_points(new_control_points);
            project_nurbs.set_knots_vector(m_knots_vector);
            project_nurbs.set_degree(m_degree);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        ENUM_NURBS curve_reparameter(const nurbs_curve<T, 1, false, -1, -1> &reparameter_function,
            nurbs_curve<T, dim, is_rational, -1, -1> &new_nurbs) const
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
            std::vector<T> refine_knots;
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
                    insert_knots.push_back(current_knot);
                    refine_knots.push_back(image_knot[0]);
                    int repeat = reparameter_function.get_knots_multiplicity(index);
                    insert_knots_degree.push_back(reparameter_function_degree - repeat);
                }
                else
                {
                    equal_knots.push_back(index);

                    int repeat1 = reparameter_function.get_knots_multiplicity(index);
                    int repeat2 = get_knots_multiplicity(knots_index);
                    int repeat = std::min(reparameter_function_degree - repeat1, m_degree - repeat2);
                    equal_knots_degree.push_back(repeat);
                }
            }

            nurbs_curve<T, dim, is_rational, -1, -1> *insert_knots_nurbs = new nurbs_curve<T, dim, is_rational, -1, -1>(*this);
            std::vector<nurbs_curve<T, dim, is_rational, -1, -1> *> bezier_curves;
            int insert_knots_size = refine_knots.size();
            Eigen::VectorX<T> temp_insert_knots(insert_knots_size);
            for (int index = 0; index < insert_knots_size; ++index)
                temp_insert_knots[index] = refine_knots[index];
            
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
                reparameter_function_insert_knots_degree.push_back(m_degree - repeat);
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

        ENUM_NURBS curve_reparameter(const nurbs_curve<T, 1, true, -1, -1> &reparameter_function,
            nurbs_curve<T, dim, true, -1, -1> &new_nurbs) const
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
            std::vector<T> refine_knots;
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
                    insert_knots.push_back(current_knot);
                    refine_knots.push_back(image_knot[0]);
                    int repeat = reparameter_function.get_knots_multiplicity(index);
                    insert_knots_degree.push_back(reparameter_function_degree - repeat);
                }
                else
                {
                    equal_knots.push_back(index);

                    int repeat1 = reparameter_function.get_knots_multiplicity(index);
                    int repeat2 = get_knots_multiplicity(knots_index);
                    int repeat = std::min(reparameter_function_degree - repeat1, m_degree - repeat2);
                    equal_knots_degree.push_back(repeat);
                }
            }

            nurbs_curve<T, dim, is_rational, -1, -1> *insert_knots_nurbs = new nurbs_curve<T, dim, is_rational, -1, -1>(*this);
            std::vector<nurbs_curve<T, dim, is_rational, -1, -1> *> bezier_curves;
            int insert_knots_size = refine_knots.size();
            Eigen::VectorX<T> temp_insert_knots(insert_knots_size);
            for (int index = 0; index < insert_knots_size; ++index)
                temp_insert_knots[index] = refine_knots[index];
            
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
                        delete bezier_curves[curve_index];
                    return ENUM_NURBS::NURBS_ERROR;
                }
                reparameter_function_insert_knots.push_back(u);
                int repeat = get_knots_multiplicity(index);
                reparameter_function_insert_knots_degree.push_back(m_degree - repeat);
            }

            nurbs_curve<T, 1, true, -1, -1> *insert_knots_function_nurbs = new nurbs_curve<T, 1, true, -1, -1>(reparameter_function);
            std::vector<nurbs_curve<T, 1, true, -1, -1> *> function_bezier_curves;
            int reparameter_function_insert_knots_size = reparameter_function_insert_knots.size();
            Eigen::VectorX<T> temp_insert_knots_vector(reparameter_function_insert_knots_size);
            for (int index = 0; index < reparameter_function_insert_knots_size; ++index)
                temp_insert_knots_vector[index] = reparameter_function_insert_knots[index];
            insert_knots_function_nurbs->refine_knots_vector(temp_insert_knots_vector);
            insert_knots_function_nurbs->decompose_to_bezier(function_bezier_curves);
            delete insert_knots_function_nurbs;

            int curves_count = function_bezier_curves.size();
            std::vector<nurbs_curve<T, dim, true, -1, -1>*> new_bezier_curves(curves_count, nullptr);
            for (int index = 0; index < curves_count; ++index)
            {
                new_bezier_curves[index] = new nurbs_curve<T, dim, true, -1, -1>();
                bezier_curves[index]->bezier_curve_reparameter(*function_bezier_curves[index], *(new_bezier_curves[index]));
                delete bezier_curves[index];
                delete function_bezier_curves[index];
            }

            Eigen::Matrix<T, dim + 1, Eigen::Dynamic> new_control_points;
            Eigen::VectorX<T> new_knots_vector;

            merge_two_curve(new_degree, new_degree, new_bezier_curves[0]->get_knots_vector(), new_bezier_curves[1]->get_knots_vector(),
            new_bezier_curves[0]->get_control_points(), new_bezier_curves[1]->get_control_points(),
            new_knots_vector, new_control_points);
            
            for (int index = 2; index < curves_count; ++index)
            {
                Eigen::Matrix<T, dim + 1, Eigen::Dynamic> nnew_control_points;
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

        // u = f(s) = alpha * s + beta
        ENUM_NURBS curve_reparameter_with_linear_function(T alpha, T beta,
            nurbs_curve<T, dim, is_rational, -1, -1> &new_nurbs) const
        {
            if (alpha == 0)
                return ENUM_NURBS::NURBS_ERROR;
            new_nurbs.set_control_points(m_control_points);
            int knots_vector_size = m_knots_vector.size();
            Eigen::VectorX<T> new_knots_vector(knots_vector_size);
            for (int index = 0; index < knots_vector_size; ++index)
            {
                new_knots_vector[index] = (m_knots_vector[index] - beta) / alpha;
            }
            new_nurbs.set_knots_vector(new_knots_vector);
            new_nurbs.set_degree(m_degree);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        // s = g(u) = (alpha * u + beta) / (gamma * u + delta)
        // alpha = a11, beta = a12, gamma = a21, dalta = a22
        ENUM_NURBS curve_reparameter_with_linear_function(Eigen::Matrix<T, 2, 2> reparameter_function,
            nurbs_curve<T, dim, true, -1, -1> &new_nurbs) const
        {
            if (reparameter_function.determinant() < DEFAULT_ERROR)
                return ENUM_NURBS::NURBS_ERROR;


            Eigen::Matrix<T, dim + 1, Eigen::Dynamic> new_control_points;
            to_ratioanl_contrl_points<T, is_rational, rows>::convert(m_control_points, new_control_points);
            int knots_vector_size = m_knots_vector.size();
            Eigen::VectorX<T> new_knots_vector(knots_vector_size);

            new_knots_vector[0] = reparameter_function(0, 0) * m_knots_vector[0] + reparameter_function(0, 1);
            new_knots_vector[0] /= (reparameter_function(1, 0) * m_knots_vector[0] + reparameter_function(1, 1));
            int current_index = 0;
            for (int index = 1; index <  knots_vector_size; ++index)
            {
                if (m_knots_vector[index] == m_knots_vector[current_index])
                    new_knots_vector[index] = new_knots_vector[current_index];
                else
                {
                    new_knots_vector[index] = reparameter_function(0, 0) * m_knots_vector[index] + reparameter_function(0, 1);
                    new_knots_vector[index] /= (reparameter_function(1, 0) * m_knots_vector[index] + reparameter_function(1, 1));
                    current_index = index;
                }    
            }
            int control_points_size = knots_vector_size - m_degree - 1;
            for (int index = 0; index < control_points_size; ++index)
            {
                T lamda_s = 1.0;
                for (int j = 1; j <= m_degree; ++j)
                    lamda_s *= (reparameter_function(1, 0) * m_knots_vector[index + j] + reparameter_function(1, 1));
                new_control_points(dim, index) /= lamda_s;
                new_control_points.block(0, index, dim, 1) /= lamda_s;
            }
            if (new_control_points(dim, 0) < 0)
                new_control_points.block(dim, 0, 1, control_points_size) = -1.0 * new_control_points.block(dim, 0, 1, control_points_size);
            new_nurbs.set_control_points(new_control_points);
            new_nurbs.set_knots_vector(new_knots_vector);
            new_nurbs.set_degree(m_degree);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        ENUM_NURBS curve_reverse()
        {
            int knots_vector_size = m_knots_vector.size();
            Eigen::VectorX<T> new_knots_vector(knots_vector_size);
            int control_points_count = knots_vector_size - m_degree - 1;
            Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points;
            new_control_points.resize(rows, control_points_count);
            T interval_begin = m_knots_vector[0];
            T interval_end = m_knots_vector[knots_vector_size - 1];
            for (int index = 0; index < knots_vector_size; ++index)
            {
                new_knots_vector[index] = interval_end + interval_begin - m_knots_vector[knots_vector_size - index - 1];
            }
            for (int index = 0; index < control_points_count; ++index)
            {
                new_control_points.col(index) = m_control_points.col(control_points_count - index - 1);
            }
            m_knots_vector = new_knots_vector;
            m_control_points = new_control_points;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        /// @brief 给定误差, 在误差范围内消去所有可以消去的节点; 此函数仅在非有理nurbs曲线是合法的(稍微修改下一误差即可在有理曲线下也是合法的, 要求参数从小到大排列, 且不重复, 稍加改造应该可以使得消去误差在整个参数域都成立)
        /// @param params 消去的曲线在参数params上误差小于给定的误差; 即(new_curve(params[i]) - curve(params[i])).norm() < E;
        /// @param error 返回的各个点的误差
        /// @param new_nurbs 节点消去后新的nurbs
        /// @param E 最大误差
        /// @return 错误码
        ENUM_NURBS remove_knots_bound_curve(const std::vector<T> &params, std::vector<T> &error, nurbs_curve<T, dim, false, -1, -1> &new_nurbs, T E = DEFAULT_ERROR) const
        {
            static_assert(is_rational == false, "rational b spline curve can not use this function");
            if constexpr (is_rational == true)
                return ENUM_NURBS::NURBS_ERROR;
            int params_count = params.size();
            for (int index = 1; index < params_count; ++index)
            {
                if (params[index] <= params[index - 1])
                    return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
            }
            
            // T inf = -1.0;
            std::vector<T> different_knots;
            std::vector<int> multiple;
            get_different_knots(m_knots_vector, m_degree, different_knots, multiple);
            if (params_count > 0)
                if (params[0] < different_knots[0] || params[params_count - 1] > different_knots.back())
                    return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
            int current_index = m_degree;
            int different_knots_count = different_knots.size();
            std::vector<T> Brs(different_knots_count - 2, 0.0);
            //计算内节点Br的值
            for (int index = 1; index < different_knots_count - 1; ++index)
            {
                int r = current_index + multiple[index];
                T Br;
                get_removal_bnd_curve<T, dim>(m_control_points, m_knots_vector, m_degree, different_knots[index], r, multiple[index], Br);
                Brs[index - 1] = Br;
                current_index = r;
            }
            //对每个基函数, 计算相关参数params的下标范围
            int points_count = m_control_points.cols();
            std::list<std::array<int, 2>> basis_params_relation;
            int first_index = 0;
            int last_index = m_degree + 1;
            int params_first_index = 0;
            int params_last_index = 0;
            for (int index = 0; index < points_count; ++index)
            {
                T u_min = m_knots_vector[first_index];
                T u_max = m_knots_vector[last_index];
                while (params[params_first_index] < u_min && params_first_index < params_count)
                {
                    ++params_first_index;
                }
                while (params[params_last_index] < u_max && params_last_index < params_count - 1)
                {
                    ++params_last_index;
                }
                
                if (params[params_last_index] >= u_max)
                {
                    --params_last_index;
                }

                if (params_first_index <= params_last_index)
                {
                    std::array<int , 2> index_segment{params_first_index, params_last_index};
                    basis_params_relation.push_back(std::move(index_segment));
                }
                else
                {
                    std::array<int , 2> index_segment{-1, -1};
                    basis_params_relation.push_back(std::move(index_segment));
                    // params_last_index = params_first_index;
                }

                ++first_index;
                ++last_index;
            }

            error.resize(params_count, 0.0);
 

            Eigen::VectorX<T> new_knots_vector = m_knots_vector;
            Eigen::Matrix<T, dim, Eigen::Dynamic> new_control_points = m_control_points;
            while (true)
            {
                auto min_it = Brs.begin();
                for (auto it = Brs.begin(); it < Brs.end(); ++it)
                {
                    if ((*it < *min_it && *it >= 0.0) || *min_it == -1.0)
                        min_it = it;
                }
                if (*min_it == -1.0)
                    break;
                int r = std::accumulate(multiple.begin(), multiple.begin() + (min_it - Brs.begin() + 2), -1);
                int &mul = multiple[min_it - Brs.begin() + 1];
                std::vector<T> temp_error = error;
                if ((m_degree + mul) % 2 == 0)
                {
                    int k = (m_degree + mul) / 2;
                    int i = r - k;
                    auto params_it = basis_params_relation.begin();
                    for (int index = 0; index < i; ++index)
                        ++params_it;
                    std::array<int , 2> &indexs = *params_it;
                    if (indexs[0] != -1 && indexs[1] != -1)
                    {
                        for (int index = indexs[0]; index <= indexs[1]; ++index)
                        {
                            T basis;
                            one_basis_function<T>(i, params[index], m_degree, new_knots_vector, basis);
                            temp_error[index] += basis * (*min_it);
                        }
                    }
                }
                else
                {
                    int k = (m_degree + mul + 1) / 2;
                    int i = r - k + 1;
                    auto params_it = basis_params_relation.begin();
                    for (int index = 0; index < i; ++index)
                        ++params_it;
                    std::array<int , 2> &indexs = *params_it;
                    if (indexs[0] != -1 && indexs[1] != -1)
                    {
                        T alpha = (new_knots_vector[r] - new_knots_vector[r - k + 1]) / (new_knots_vector[r - k + m_degree + 2] - new_knots_vector[r - k + 1]);
                        for (int index = indexs[0]; index <= indexs[1]; ++index)
                        {
                            T basis;
                            one_basis_function<T>(i, params[index], m_degree, new_knots_vector, basis);
                            temp_error[index] += basis * (*min_it) * (1.0 - alpha);
                        }
                    }
                }
            
                bool flag = true;
                for (auto &e : temp_error)
                {
                    if (e > E)
                        flag = false;
                }
                if (flag == true)
                {
                    error = temp_error;
                    T remove_knot = new_knots_vector[r];
                    int new_points_count = new_control_points.cols();
                    int beg = std::max(r - m_degree, m_degree + 1);
                    int end = std::min(r + m_degree - mul + 1, new_points_count - 1);

                    remove_curve_knots_no_check_error<T, dim, false>(r, 1, m_degree, mul, new_knots_vector, new_control_points);
                    if (new_knots_vector.size() == 2 * (m_degree + 1))
                        break;
                    //更新基函数对应的参数的下标
                    auto params_it = basis_params_relation.begin();
                    for (int index = 0; index < r - m_degree - 1; ++index)
                        ++params_it;
                    for (int index = r - m_degree - 1; index <= r - mul; ++index)
                    {
                        if ((*params_it)[0] == -1)
                        {
                            auto current_it = params_it++;
                            *current_it = *(params_it);
                        }
                        else
                        {
                            auto current_it = params_it++;
                            (*current_it)[1] = (*params_it)[1];
                        }
                        // ++params_it;
                    }
                    basis_params_relation.erase(params_it);


                    int different_knots_size = different_knots.size() - 1;
                    int knot_index_start = m_degree + 1;
                    int kont_index_end = m_degree + 1 + multiple[1] - 1;

                    for (int index = 1; index < different_knots_size; ++index)
                    {
                        if ((knot_index_start >= beg && knot_index_start <= end) || (kont_index_end >= beg && kont_index_end <= end))
                        {
                            if (different_knots[index] == remove_knot && 1 == mul)
                            {
                                end = std::min(r + m_degree - mul + 1, new_points_count - 2);
                                --kont_index_end;
                                --mul;
                                knot_index_start = kont_index_end + 1;
                                kont_index_end = knot_index_start + multiple[index + 1] - 1;
                                continue;
                            }
                            if (different_knots[index] == remove_knot)
                            {
                                end = std::min(r + m_degree - mul + 1, new_points_count - 2);
                                kont_index_end -= 1;
                                --multiple[index];
                            }

                            T Br;
                            get_removal_bnd_curve<T, dim>(new_control_points, new_knots_vector, m_degree, different_knots[index], kont_index_end, multiple[index], Br);
                            Brs[index - 1] = Br;
                        }
                        knot_index_start = kont_index_end + 1;
                        kont_index_end = knot_index_start + multiple[index + 1] - 1;
                    }

                    if (mul == 0)
                    {
                        auto knot_it = std::find(different_knots.begin(), different_knots.end(), remove_knot);
                        multiple.erase(multiple.begin() + (knot_it - different_knots.begin()));
                        different_knots.erase(knot_it);
                        Brs.erase(min_it);
                    }

                }
                else
                {
                    *min_it = -1.0;
                }
            }
            new_nurbs.set_control_points(new_control_points);
            new_nurbs.set_knots_vector(new_knots_vector);
            new_nurbs.set_degree(m_degree);
            
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 将nurbs曲线的维数提升一维
        /// @param index 在第index个坐标之前插入0
        /// @return 错误码
        ENUM_NURBS dimension_elevate(int index, nurbs_curve<T, dim + 1, is_rational, -1, -1> &new_nurbs) const
        {
            int points_count = m_control_points.cols();
            Eigen::Matrix<T, rows + 1, Eigen::Dynamic> new_control_points(rows + 1, points_count);

            // if (index > 0)
            new_control_points.block(0, 0, index, points_count) = m_control_points.block(0, 0, index, points_count);
            new_control_points.block(index, 0, 1, points_count).setConstant(0.0);
            new_control_points.block(index + 1, 0, rows - index, points_count) = m_control_points.block(index, 0, rows - index, points_count);


            new_nurbs.set_knots_vector(m_knots_vector);
            new_nurbs.set_control_points(new_control_points);
            new_nurbs.set_degree(m_degree);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS move(const Eigen::Vector<T, dim> &translate_vector, nurbs_curve<T, dim, is_rational, -1, -1> &new_nurbs) const
        {
            Eigen::Matrix<T, rows, Eigen::Dynamic> new_control_points = m_control_points;
            int points_count = m_control_points.cols();
            for (int index = 0; index < points_count; ++index)
            {
                T w = 1.0;
                if constexpr (is_rational == true)
                    w = new_control_points(dim, index);
                new_control_points.template block<dim, 1>(0, index) += (translate_vector * w);
            }
            new_nurbs.set_control_points(new_control_points);
            new_nurbs.set_degree(m_degree);
            new_nurbs.set_knots_vector(m_knots_vector);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


    private:

        ENUM_NURBS bezier_curve_reparameter(const nurbs_curve<T, 1, false, -1, -1> &reparameter_function,
            nurbs_curve<T, dim, is_rational, -1, -1> &new_nurbs) const
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

        ENUM_NURBS bezier_curve_reparameter(const nurbs_curve<T, 1, true, -1, -1> &reparameter_function,
            nurbs_curve<T, dim, true, -1, -1> &new_nurbs) const 
        {
            Eigen::Matrix<T, dim + 1, Eigen::Dynamic> rational_control_points;
            to_ratioanl_contrl_points<T, is_rational, rows>::convert(m_control_points, rational_control_points);

            std::array<Eigen::Vector<T, 1>, 2> ends_points;
            int knots_vector_size = m_knots_vector.size();
            reparameter_function.get_ends_point(ends_points);

            if (std::abs(ends_points[0][0] - m_knots_vector[0]) > DEFAULT_ERROR)
                return ENUM_NURBS::NURBS_ERROR;
            if (std::abs(ends_points[1][0] - m_knots_vector[knots_vector_size - 1]) > DEFAULT_ERROR)
            return ENUM_NURBS::NURBS_ERROR;

            Eigen::VectorX<T> function_knots_vector = reparameter_function.get_knots_vector();
            Eigen::Matrix<T, 2, Eigen::Dynamic> function_control_points = reparameter_function.get_control_points();
            int function_degree = reparameter_function.get_degree();

            Eigen::Matrix<T, dim + 1, Eigen::Dynamic> new_control_points;
            int new_degree = m_degree * function_degree;
            new_control_points.resize(dim + 1, new_degree + 1);
            new_control_points.setConstant(0.0);
            new_control_points.col(0) = rational_control_points.col(0);
            new_control_points.col(new_degree) = rational_control_points.col(m_degree);
            T low = function_knots_vector[0];
            T high = function_knots_vector[function_degree + 1];
            reparameter_bezier_curve<T, dim + 1>(m_degree, m_knots_vector, rational_control_points, function_degree, function_knots_vector,
                function_control_points, function_knots_vector[function_degree + 1] - function_knots_vector[0], new_control_points);

            
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

    };

    template<typename T, int dim, bool is_rational, int points_count, int degree>
    struct geo_traits<nurbs_curve<T, dim, is_rational, points_count, degree> >
    {
        static constexpr int point_size = is_rational ? dim + 1 : dim;
        using type = nurbs_curve<T, dim, is_rational, points_count, degree>;
        using point_number_type = T;
        using point_type = typename  Eigen::Vector<T, dim> ;
    };



}

