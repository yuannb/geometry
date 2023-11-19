#pragma once
#include "nurbs_tool.h"
#include "bezier_surface.h"
#include "nurbs_curve.h"
#include <concepts>
#include "surface.h"
namespace tnurbs
{
    using namespace tnurbs;
    //TODO: U向控制点或者V向控制点个数不固定, 另一个方向固定

    /*
                                [P_00, P_01 ,..., P_0n        [u_0  
                                P_10, P_11 ,..., P_1n         u_1
    [v_0, v_1 ,..., v_m]   *     .     .     .    .      *      .
                                .     .     .    .             .  
                                .     .     .    .             .
                                P_m0, P_m1 ,..., P_mn]        u_n]

    */


    /// @brief nurbs surface class
    /// @tparam T : double float int...
    /// @tparam dim : nurbs所在欧氏空间的维数
    /// @tparam rows : U向控制点的个数(不严谨)
    /// @tparam cols : V向控制点的个数(不严谨)
    /// @tparam u_degree : nurbs曲面的u向阶数
    /// @tparam v_degree : nurbs曲面的v向阶数
    /// @tparam is_rational : 是否时有理nurbs
    template<typename T, int dim, int rows, int cols, int u_degree, int v_degree, bool is_rational>
    class nurbs_surface
    {
        static constexpr int point_size = is_rational ? dim + 1 : dim;
        //TODO:
        int get_u_degree();
        int get_v_degree();
        ENUM_NURBS reverse_uv();
        ENUM_NURBS scale_surface(const Eigen::Vector<T, dim> &center, const Eigen::Matrix<T, dim, dim> &mat, const Eigen::Vector<T, dim> &scale_factory);
        Eigen::VectorX<T> get_u_knots() const;
        Eigen::VectorX<T> get_v_knots() const;
        ENUM_NURBS set_u_degree(int degree);
        ENUM_NURBS set_v_degree(int degree);
        ENUM_NURBS get_u_different_knots(std::vector<T> &vec) const;
        ENUM_NURBS get_v_different_knots(std::vector<T> &vec) const;
        ENUM_NURBS set_uv_degree(int u_degree_, int v_degree_);
        ENUM_NURBS set_control_points(Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> &points);
        ENUM_NURBS point_on_surface(T u, T v, Eigen::Vector<T, dim> &point) const;
        ENUM_NURBS set_uv_knots(const Eigen::VectorX<T> &u_knots, const Eigen::VectorX<T> &v_knots);
        ENUM_NURBS degree_elevate(int t, ENUM_DIRECTION direction);
        ENUM_NURBS refine_knots_vector(const Eigen::VectorX<T> &insert_knots, ENUM_DIRECTION direction);
        ENUM_NURBS get_uv_knots_end(std::array<T, 2> &u_knots_end, std::array<T, 2> &v_knots_end) const;
         // u = f(s) = alpha * s + beta
        ENUM_NURBS surface_reparameter_with_linear_function(T alpha, T beta, ENUM_DIRECTION direction, 
            nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &new_nurbs_surface);
        bool is_u_closed() const;
        bool is_v_closed() const;
        ENUM_NURBS get_uv_box(Box<T, 2> &uv_box) const;
        template<int n>
        ENUM_NURBS derivative_on_surface(T u, T v, Eigen::Matrix<Eigen::Vector<T, dim>, n + 1, n + 1> &result) const;
    };


    /// @brief nurbs surface class
    /// @tparam T : double float int...
    /// @tparam dim : nurbs所在欧氏空间的维数
    /// @tparam u_degree : nurbs曲面的u向阶数
    /// @tparam v_degree : nurbs曲面的v向阶数
    /// @tparam is_rational : 是否时有理nurbs
    template<typename T, int dim, int u_degree, int v_degree, bool is_rational>
    class nurbs_surface<T, dim, -1, -1, u_degree, v_degree, is_rational>
    {
    private:
        static constexpr int point_size = is_rational ? dim + 1 : dim;
        Eigen::VectorX<T> m_u_knots_vector;
        Eigen::VectorX<T> m_v_knots_vector;

        //m_control_points[i](, j)为P_ij
        Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> m_control_points;

    public:
        nurbs_surface(const Eigen::VectorX<T> &u_knots_vector, const Eigen::VectorX<T> &v_knots_vector,
        const Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> &control_points) :
        m_u_knots_vector(u_knots_vector), m_v_knots_vector(v_knots_vector), m_control_points(control_points)
        {
        }

        /// @brief 计算曲面上的点
        /// @param u u向参数
        /// @param v v向参数
        /// @param point 曲面上的点
        /// @return ENUM_NURBS错误码
        ENUM_NURBS point_on_surface(T u, T v, Eigen::Vector<T, dim> &point) const
        {
            int uspan = -1, vspan = -1;
            find_span<T, u_degree>(u, m_u_knots_vector, uspan);
            find_span<T, v_degree>(v, m_v_knots_vector, vspan);
            Eigen::Vector<T, u_degree + 1> nu;
            Eigen::Vector<T, v_degree + 1> nv;
            basis_functions<T, u_degree>(uspan, u, m_u_knots_vector, nu);
            basis_functions<T, v_degree>(vspan, v, m_v_knots_vector, nv);

            Eigen::Matrix<T, point_size, v_degree + 1> temp;
            for (int index = 0; index <= v_degree; ++index)
            {
                temp.col(index) = m_control_points[vspan + index - v_degree].block(0, uspan - u_degree, point_size, u_degree + 1) * nu;
            }
            point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(temp * nv);

            return ENUM_NURBS::NURBS_SUCCESS;
        }

        /// @brief 计算nurbs surface在参数(u, v)处直到n阶偏导数
        /// @tparam n 最高阶数
        /// @param u 参数u
        /// @param v 参数v
        /// @param result out_put_param result(i, j)表示S(u, v)处的u向i次偏导数, 
        /// remark : v向j次偏导数;当(i + j) > m_v_degree + m_u_degree || i > m_u_degree ||j > m_v_degree 时, 原则上来说是非法的, 此处赋值为零向量
        /// @return ENUM_NURBS错误码
        template<int n>
        ENUM_NURBS derivative_on_surface(T u, T v, Eigen::Matrix<Eigen::Vector<T, dim>, n + 1, n + 1> &result) const
        {
            Eigen::Matrix<Eigen::Vector<T, point_size>, n + 1, n + 1> ration_result;
            constexpr int du = n > u_degree ? u_degree : n;
            constexpr int dv = n > v_degree ? v_degree : n;
            constexpr int u_order = u_degree + 1;
            constexpr int v_order = v_degree + 1;
            Eigen::Vector<T, point_size> zero_vector;
            zero_vector.setConstant(0.0);
            for (int k = u_order; k <= n; ++k)
            {
                int count = n - k;
                for (int l = 0; l <= count; ++l)
                    ration_result(k, l) = zero_vector;
            }

            for (int l = v_order; l <= n; ++l)
            {
                int count = n - l;
                for (int k = 0; k <= count; ++k)
                    ration_result(k, l) = zero_vector;
            }

            int uspan = -1;
            int vspan = -1;
            find_span<T, u_degree>(u, m_u_knots_vector, uspan);
            find_span<T, v_degree>(v, m_v_knots_vector, vspan);
            Eigen::Matrix<T, u_degree + 1, n + 1> nu;
            Eigen::Matrix<T, v_degree + 1, n + 1> nv;
            ders_basis_funs<T, u_degree, du>(uspan, u, m_u_knots_vector, nu);
            ders_basis_funs<T, v_degree, dv>(vspan, v, m_v_knots_vector, nv);
            Eigen::Matrix<T, point_size,  v_degree + 1> temps;
            for (int k = 0; k <= du; ++k)
            {
                temps.setConstant(0.0);
                for (int s = 0; s <= v_degree; ++s)
                {
                    temps.col(s) = m_control_points[vspan - v_degree + s].template block<point_size, u_degree + 1>(0, uspan - u_degree) * nu.col(k);
                }
                int dd = std::min(n - k, dv);
                for (int l = 0; l <= dd; ++l)
                {
                    ration_result(k, l) = temps * nv.col(l);
                }
            }
            result = project_derivs_point<T, is_rational, point_size, n>::project_point_to_euclidean_space(ration_result);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        /// @brief 计算nurbs surface在参数(u, v)处直到n阶偏导数
        /// @param n 最高阶数
        /// @param u 参数u
        /// @param v 参数v
        /// @param result out_put_param result(i, j)表示S(u, v)处的u向i次偏导数, 
        /// remark : v向j次偏导数;当(i + j) > m_v_degree + m_u_degree || i > m_u_degree ||j > m_v_degree 时, 原则上来说是非法的, 此处赋值为零向量
        /// @return ENUM_NURBS错误码
        ENUM_NURBS derivative_on_surface(int n, T u, T v, Eigen::MatrixX<Eigen::Vector<T, dim>> &result) const
        {
            Eigen::MatrixX<Eigen::Vector<T, point_size>> ration_result(n + 1, n + 1);
            int du = std::min(n, u_degree);
            int dv = std::min(n, v_degree);
            Eigen::Vector<T, point_size> zero_vector;
            zero_vector.setConstant(0.0);
            for (int k = u_degree + 1; k <= n; ++k)
            {
                int count = n - k;
                for (int l = 0; l <= count; ++l)
                    ration_result(k, l) = zero_vector;
            }

            for (int l = v_degree + 1; l <= n; ++l)
            {
                int count = n - l;
                for (int k = 0; k <= count; ++k)
                    ration_result(k, l) = zero_vector;
            }

            int uspan = -1;
            int vspan = -1;
            find_span<T, u_degree>(u, m_u_knots_vector, uspan);
            find_span<T, v_degree>(v, m_v_knots_vector, vspan);
            Eigen::Matrix<T, u_degree + 1, Eigen::Dynamic> nu;
            Eigen::Matrix<T, v_degree + 1, Eigen::Dynamic> nv;
            ders_basis_funs<T, u_degree>(uspan, du, u, m_u_knots_vector, nu);
            ders_basis_funs<T, v_degree>(vspan, dv, v, m_v_knots_vector, nv);
            Eigen::Matrix<T, point_size,  v_degree + 1> temps;
            for (int k = 0; k <= du; ++k)
            {
                temps.setConstant(0.0);
                for (int s = 0; s <= v_degree; ++s)
                {
                    temps.col(s) = m_control_points[vspan - v_degree + s].template block<point_size, u_degree + 1>(0, uspan - u_degree) * nu.col(k);
                }
                int dd = std::min(n - k, dv);
                for (int l = 0; l <= dd; ++l)
                {
                    ration_result(k, l) = temps * nv.col(l);
                }
            }
            result = project_derivs_point<T, is_rational, point_size, -1>::project_point_to_euclidean_space(ration_result);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 计算nurbs surface在参数(u, v)处直到n阶偏导数
        /// @param n1 u向求导最高阶数
        /// @param n2 v向求导最高阶数
        /// @param u 参数u
        /// @param v 参数v
        /// @param result out_put_param result(i, j)表示S(u, v)处的u向i次偏导数, 
        /// remark : v向j次偏导数;当(i + j) > m_v_degree + m_u_degree || i > m_u_degree ||j > m_v_degree 时, 原则上来说是非法的, 此处赋值为零向量
        /// @return ENUM_NURBS错误码
        ENUM_NURBS derivative_on_surface(int n1, int n2, T u, T v, Eigen::MatrixX<Eigen::Vector<T, dim>> &result) const
        {
            Eigen::MatrixX<Eigen::Vector<T, point_size>> ration_result(n1 + 1, n2 + 1);

            int uspan = -1;
            int vspan = -1;
            find_span<T>(u, u_degree, m_u_knots_vector, uspan);
            find_span<T>(v, v_degree, m_v_knots_vector, vspan);
            Eigen::MatrixX<T> nu, nv;
            ders_basis_funs<T>(uspan, n1, u_degree, u, m_u_knots_vector, nu);
            ders_basis_funs<T>(vspan, n2, v_degree, v, m_v_knots_vector, nv);
            Eigen::Matrix<T, point_size,  Eigen::Dynamic> temps;
            temps.resize(point_size, v_degree + 1);
            for (int k = 0; k <= n1; ++k)
            {
                temps.setConstant(0.0);
                for (int s = 0; s <= v_degree; ++s)
                {
                    temps.col(s) = m_control_points[vspan - v_degree + s].block(0, uspan - u_degree, point_size, u_degree + 1) * nu.col(k);
                }
                for (int l = 0; l <= n2; ++l)
                {
                    ration_result(k, l) = temps * nv.col(l);
                }
            }
            result = project_derivs_point<T, is_rational, point_size, -1>::project_point_to_euclidean_space(ration_result);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        /// @brief 计算nurbs surface在参数(u, v)处直到n阶偏导数
        /// @tparam n1 u向求导最高阶数
        /// @tparam n2 v向求导最高阶数
        /// @param u 参数u
        /// @param v 参数v
        /// @param result out_put_param result(i, j)表示S(u, v)处的u向i次偏导数, 
        /// remark : v向j次偏导数;当(i + j) > m_v_degree + m_u_degree || i > m_u_degree ||j > m_v_degree 时, 原则上来说是非法的, 此处赋值为零向量
        /// @return ENUM_NURBS错误码
        template<int n1, int n2>
        ENUM_NURBS derivative_on_surface(T u, T v, Eigen::Matrix<Eigen::Vector<T, dim>, n1 + 1, n2 + 1> &result) const
        {
            Eigen::Matrix<Eigen::Vector<T, point_size>, n1 + 1, n2 + 1> ration_result;

            int uspan = -1;
            int vspan = -1;
            find_span<T, u_degree>(u, m_u_knots_vector, uspan);
            find_span<T, v_degree>(v, m_v_knots_vector, vspan);
            Eigen::Matrix<T, Eigen::Dynamic, n1 + 1> nu;
            Eigen::Matrix<T, Eigen::Dynamic, n2 + 1> nv;
            ders_basis_funs<T, u_degree, n1>(uspan, u, m_u_knots_vector, nu);
            ders_basis_funs<T, v_degree, n2>(vspan, v, m_v_knots_vector, nv);
            Eigen::Matrix<T, point_size,  Eigen::Dynamic> temps;
            temps.resize(point_size, v_degree + 1);
            for (int k = 0; k <= n1; ++k)
            {
                temps.setConstant(0.0);
                for (int s = 0; s <= v_degree; ++s)
                {
                    temps.col(s) = m_control_points[vspan - v_degree + s].block(0, uspan - u_degree, point_size, u_degree + 1) * nu.col(k);
                }
                for (int l = 0; l <= n2; ++l)
                {
                    ration_result(k, l) = temps * nv.col(l);
                }
            }
            result = project_derivs_point<T, is_rational, point_size, -1>::project_point_to_euclidean_space(ration_result);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        // template<int r1, int r2, int s1, int s2, int d>
        // ENUM_NURBS surface_derivs_alg2(T u, T v, Eigen::Matrix<Eigen::Vector<T, dim>, d + 1, d + 1> &SKL)
        // {
        //     //先将导数为0位置设置为0
        //     constexpr int du = std::min(d, u_degree);
        //     constexpr int dv = std::min(d, v_degree);
        //     Eigen::Matrix<Eigen::Vector<T, point_size>, d + 1, d + 1> &SKLT;
        //     Eigen::Vector<T, point_size> zero;
        //     zero.setConstant(0.0);
        //     for (int k = u_degree + 1; k <= d; ++k)
        //         for (int l = 0; l <= d - k; ++l)
        //             SKLT(k, l) = zero;
        //     for (int l = v_degree + 1; l <= d; ++l)
        //         for(int k = 0; k <= d - l; ++k)
        //             SKLT(k, l) = zero;
        //     int uspan = -1, vspan = -1;
        //     find_span<T, u_degree>(u, m_u_knots_vector, uspan);
        //     find_span<T, v_degree>(v, m_v_knots_vector, vspan);
        //     Eigen::Matrix<T, u_degree + 1, u_degree + 1> nu;
        //     Eigen::Matrix<T, v_degree + 1, v_degree + 1> nv;
        //     all_basis_functions<T, u_degree>(uspan, u, m_u_knots_vector, nu);
        //     all_basis_functions<T, v_degree>(vspan, v, m_v_knots_vector, nv);
        //     surface_derv_cpts_matrix<T, point_size, u_degree, v_degree, d> PKL;
        //     surface_deriv_cpts<T, u_degree, v_degree, r1, r2, s1, s2, point_size, d>(m_u_knots_vector.size(), m_v_knots_vector.size(),
        //                 m_control_points, m_u_knots_vector, m_v_knots_vector, PKL);
        //     for (int k = 0; k <= du; ++k)
        //     {
        //         int dd = std::min(d - k, dv);
        //         for (int l = 0; l <= dd; ++l)
        //         {
        //             Eigen::Matrix<T, point_size, Eigen::Dynamic>temp;
        //             temp.resize(point_size, v_degree - l + 1);
        //             int col_size = v_degree - l + 1;
        //             for (int j = 0; j < col_size; ++j)
        //             {
        //                 temp.col(j) = PKL(k, l)[j] * nu.col(k);
        //             }
        //             SKLT(k, l) = temp * nv.block(0, l, v_degree - l + 1, 1)col();
        //         }
        //     }
        //     SKL = project_derivs_point<T, is_rational, point_size, d>::project_point_to_euclidean_space(SKLT);
        //     return ENUM_NURBS::NURBS_SUCCESS;
        // }

        //此函数仅内部使用(因为u_degree和v_degree都是常量, 因此实际上不能交换曲面的u向和v向(u_degree = v_degree时可以交换))
        //这个函数转置控制点矩阵,交换u向和v向节点矢量, 但是不改变u_degree的值和v_degree的值
        ENUM_NURBS reverse_uv()
        {
            int rows = m_control_points.rows();
            int cols = m_control_points[0].cols();
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points;
            new_control_points.resize(cols);
            for (int col_index = 0; col_index < cols; ++col_index)
            {
                new_control_points[col_index].resize(point_size, rows);
                for (int row_index = 0; row_index < rows; ++row_index)
                {
                    new_control_points[col_index].col(row_index) = m_control_points[row_index].col(col_index);
                }
            }
            m_u_knots_vector.swap(m_v_knots_vector);
            m_control_points = new_control_points;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        /// @brief 将nurbs曲面的一个方向插入节点
        /// @param uv 插入的节点
        /// @param r 插入节点的次数
        /// @param direction 插入节点的方向
        /// @return ENUM_NURBS错误码
        ENUM_NURBS surface_knots_insert(T uv, int r, ENUM_DIRECTION direction)
        {
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            int degree = direction == ENUM_DIRECTION::V_DIRECTION ? v_degree : u_degree;
            int rows = m_control_points.rows();
            int cols = m_control_points[0].cols();
            int knots_size = m_u_knots_vector.size();
            // if (std::abs(uv - m_u_knots_vector[knots_size - 1]) < KNOTS_VECTOR_EPS)
            if (uv == m_u_knots_vector[knots_size - 1])
                return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
            // if (std::abs(uv - m_u_knots_vector[0]) < KNOTS_VECTOR_EPS)
            if (uv == m_u_knots_vector[0])
                return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
            Eigen::Matrix<T, point_size, Eigen::Dynamic> new_control_points;
            int span = -1;
            if(find_span<T>(uv, degree, m_u_knots_vector, span) != ENUM_NURBS::NURBS_SUCCESS)
                return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
            int current_index = span;
            int repeat = 0;
            while (current_index >= 0)
            {
                // if (std::abs(m_u_knots_vector[current_index--] - uv) < KNOTS_VECTOR_EPS)
                if (m_u_knots_vector[current_index--] == uv)
                {
                    ++repeat;
                    continue;
                }
                break;
            }

            for (int row_index = 0; row_index < rows; ++row_index)
            {
                auto control_points = m_control_points[row_index];
                curve_knots_insert<T, point_size>(r, cols, degree, m_u_knots_vector, control_points, 
                    uv, span, repeat, new_control_points);
                m_control_points[row_index] = new_control_points;
            }

            Eigen::VectorX<T> new_knots_vector(knots_size + r);
            new_knots_vector.block(0, 0, span + 1, 1) = m_u_knots_vector.block(0, 0, span + 1, 1);
            for (int i = 1; i <= r; ++i)
                new_knots_vector[span + i] = uv;
            new_knots_vector.block(span + 1 + r, 0, knots_size - span - 1, 1) = m_u_knots_vector.block(span + 1, 0, knots_size - span - 1, 1);
            m_u_knots_vector = new_knots_vector;
            
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            return ENUM_NURBS::NURBS_SUCCESS;
        
        }


        ENUM_NURBS refine_knots_vector(const Eigen::VectorX<T> &insert_knots, ENUM_DIRECTION direction)
        {
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            int degree = direction == ENUM_DIRECTION::V_DIRECTION ? v_degree : u_degree;
            int rows = m_control_points.rows();
            // int cols = m_control_points[0].cols();
            // int knots_size = m_u_knots_vector.size();
            Eigen::Matrix<T, point_size, Eigen::Dynamic> new_control_points;
            Eigen::VectorX<T> new_knots_vector;
            int knots_vector_size = m_u_knots_vector.size();
            for (int row_index = 0; row_index < rows; ++row_index)
            {
                auto control_points = m_control_points[row_index];
                refine_knots_vector_curve<T, point_size>(knots_vector_size, degree, m_u_knots_vector, control_points,
                    insert_knots, new_knots_vector, new_control_points);
                m_control_points[row_index] = new_control_points;
            }
            m_u_knots_vector = new_knots_vector;
            
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            return ENUM_NURBS::NURBS_SUCCESS;
        
        }


        /// @brief 将nurbs曲线分解成bezier曲面
        /// @param bezier_curves 分解的bezier曲线, 内存用户释放
        /// @return ENUM_NURBS错误码
        ENUM_NURBS decompose_to_bezier(Eigen::MatrixX<bezier_surface<T, dim, u_degree + 1, v_degree + 1, is_rational> *> &bezier_surfaces)
        {
            //先u向加细成bezier
            Eigen::VectorX<T> new_knots_vector;
            Eigen::VectorX<Eigen::Matrix<T, point_size, u_degree + 1>> new_control_points;
            int u_interval_count = -1, v_interval_count = -1;
            find_interval_segment_count<T>(u_degree, m_u_knots_vector, u_interval_count);
            find_interval_segment_count<T>(v_degree, m_v_knots_vector, v_interval_count);
            new_control_points.resize(u_interval_count);
            int rows_count = m_control_points.rows() - 1;

            //(B_ij)(k,l)表示第(i, j)段beizer曲面的第(k, l)个控制点
            Eigen::MatrixX<Eigen::Matrix<Eigen::Vector<T, point_size>, v_degree + 1, u_degree + 1>> bezier_surface_control_points(u_interval_count, v_interval_count);

            //T[i].col(j)表示第i列曲线的第j个控制点
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> temp_control_points;
            int cols_count = u_interval_count * u_degree;
            temp_control_points.resize(cols_count + 1);
            int rsize = m_v_knots_vector.size() - v_degree - 1;
            for (int i = 0; i <= cols_count; ++i)
            {
                temp_control_points[i].resize(point_size, rsize);
            }

            for (int index = 0; index <= rows_count; ++index)
            {
                decompose_curve_to_bezier<T, point_size, u_degree>(u_interval_count, m_u_knots_vector, m_control_points[index], new_knots_vector, new_control_points);
                for (int i = 0; i < u_interval_count; ++i)
                {
                    for (int k = 0; k < u_degree; ++k)
                    {       
                        temp_control_points[i * u_degree + k].col(index) = new_control_points[i].col(k);
                    }
                }
                temp_control_points[cols_count].col(index) = new_control_points[u_interval_count - 1].col(u_degree);           
            }

            //在加细v向
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> v_new_control_points;
            v_new_control_points.resize(v_interval_count);
            for (int index = 0; index < cols_count; ++index)
            {
                int br = index % v_degree;
                int r = index / v_degree;
                decompose_curve_to_bezier<T, point_size>(v_degree, v_interval_count, m_v_knots_vector, temp_control_points[index], new_knots_vector, v_new_control_points);
                
                for (int i = 0; i < v_interval_count; ++i)
                {
                    for (int k = 0; k <= v_degree; ++k)
                    {
                        bezier_surface_control_points(r, i)(k, br) = v_new_control_points[i].col(k);
                    }
                    if (r != 0 && br == 0)
                    {
                        for (int k = 0; k <= v_degree; ++k)
                        {
                            bezier_surface_control_points(r - 1, i)(k, u_degree) = v_new_control_points[i].col(k);
                        }    
                    }
                }        
            }
            decompose_curve_to_bezier<T, point_size>(v_degree, v_interval_count, m_v_knots_vector, temp_control_points[cols_count], new_knots_vector, v_new_control_points);
            for (int i = 0; i < v_interval_count; ++i)
            {
                for (int k = 0; k <= v_degree; ++k)
                {
                    bezier_surface_control_points(u_interval_count - 1, i)(k, u_degree) = v_new_control_points[i].col(k);
                }
            }

            bezier_surfaces.resize(u_interval_count, v_interval_count);
            for (int i = 0; i < u_interval_count; ++i)
            {
                for (int j = 0; j < v_interval_count; ++j)
                {
                    bezier_surface<T, dim, u_degree + 1, v_degree + 1, is_rational> *bs = new bezier_surface<T, dim, u_degree + 1, v_degree + 1, is_rational>(bezier_surface_control_points(i, j));
                    bezier_surfaces(i, j) = bs;
                }
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
        ENUM_NURBS remove_knots(T u, int count, ENUM_DIRECTION direction, int &time, T error = DEFAULT_ERROR)
        {
            int degree = u_degree;
            if (direction == ENUM_DIRECTION::V_DIRECTION)
            {
                reverse_uv();
                degree = u_degree;
            }
    
            int knots_size = m_u_knots_vector.size();
            if (u == m_u_knots_vector[0] || u == m_u_knots_vector[knots_size - 1])
            {
                time = 0;
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            int span = -1;
            find_span<T>(u, degree, m_u_knots_vector, span);
            if (m_u_knots_vector[span] != u)
                return ENUM_NURBS::NURBS_ERROR;
            
            int repeat = 1;
            int current = span - 1;
            while (current >= 0 && m_u_knots_vector[current] == u)
            {
                repeat += 1;
                current -= 1;
            }
            count = std::min(repeat, count);
            Eigen::VectorX<T> new_knots_vector;
            int rows = m_control_points.rows();
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points = m_control_points;
            for (time = 0; time < count; ++time)
            {
                Eigen::VectorX<T> knots_vector;
                for(int row = 0; row < rows; ++row)
                {
                    knots_vector = m_u_knots_vector;
                    int ture_or_false = -1;
                    Eigen::Matrix<T, point_size, Eigen::Dynamic> points = new_control_points[row];
                    remove_curve_knots<T, point_size, is_rational>(span, 1, degree, repeat - time, knots_vector, points, ture_or_false, error);
                    if (ture_or_false != 1)
                    {
                        if (direction == ENUM_DIRECTION::V_DIRECTION)
                            reverse_uv();
                        return ENUM_NURBS::NURBS_SUCCESS;
                    }
                    new_control_points[row] = points;
                }
                m_u_knots_vector = knots_vector;
                m_control_points = new_control_points;
            }
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        
        /// @brief 升阶
        /// @param t 提升的阶数
        /// @param direction 提升的方向
        /// @param surface 提升后的曲面
        /// @return ENUM_NURBS错误码
        ENUM_NURBS degree_elevate(int t, ENUM_DIRECTION direction, nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &surface)
        {
            int degree = u_degree;
            if (direction == ENUM_DIRECTION::V_DIRECTION)
            {
                degree = v_degree;
                reverse_uv();
            }

            Eigen::VectorX<T> new_knots_vector;
            int count = m_control_points.rows();
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points(count);
            ENUM_NURBS flag;
            for (int index = 0; index < count; ++index)
            {
                flag = degree_elevate_curve<T, point_size>(t, degree, m_u_knots_vector, m_control_points[index], new_knots_vector, new_control_points[index]);
                if (flag != ENUM_NURBS::NURBS_SUCCESS)
                {
                    if (direction == ENUM_DIRECTION::V_DIRECTION)
                        reverse_uv();
                    return flag;
                }
            }

            if (direction == ENUM_DIRECTION::V_DIRECTION)
            {
                surface.set_uv_degree(v_degree + t, u_degree);
                surface.set_control_points(new_control_points);
                surface.set_uv_knots(new_knots_vector, m_v_knots_vector);
                surface.reverse_uv();
                reverse_uv();
                return ENUM_NURBS::NURBS_SUCCESS;
            }

            //else
            surface.set_uv_degree(u_degree + t, v_degree);
            surface.set_control_points(new_control_points);
            surface.set_uv_knots(new_knots_vector, m_v_knots_vector);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 将nurbs surface的某个方向降低一阶
        /// @param direction 降阶的方向
        /// @param surface 降阶后的nurbs surface
        /// @param error 误差
        /// @return ENUM_NURBS错误码
        ENUM_NURBS degree_reduce(ENUM_DIRECTION direction , nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &surface, T error = DEFAULT_ERROR)
        {

            int degree = u_degree;
            if (direction == ENUM_DIRECTION::V_DIRECTION)
            {
                degree = v_degree;
                reverse_uv();
            }

            Eigen::VectorX<T> new_knots_vector;
            int count = m_control_points.rows();
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points(count);
            ENUM_NURBS flag;
            for (int index = 0; index < count; ++index)
            {
                flag = degree_reduce_curve<T, point_size, is_rational>(degree, m_u_knots_vector, m_control_points[index], new_knots_vector, new_control_points[index], error);
                if (flag != ENUM_NURBS::NURBS_SUCCESS)
                {
                    if (direction == ENUM_DIRECTION::V_DIRECTION)
                        reverse_uv();
                    return flag;
                }
            }

            if (direction == ENUM_DIRECTION::V_DIRECTION)
            {
                surface.set_uv_degree(v_degree -1, u_degree);
                surface.set_control_points(new_control_points);
                surface.set_uv_knots(new_knots_vector, m_v_knots_vector);
                surface.reverse_uv();
                reverse_uv();
                return ENUM_NURBS::NURBS_SUCCESS;
            }

            //else
            surface.set_uv_degree(u_degree - 1, v_degree);
            surface.set_control_points(new_control_points);
            surface.set_uv_knots(new_knots_vector, m_v_knots_vector);
            return ENUM_NURBS::NURBS_SUCCESS;

        }


    };

    /// @brief nurbs surface class
    /// @tparam T : double float int...
    /// @tparam dim : nurbs所在欧氏空间的维数
    /// @tparam is_rational : 是否时有理nurbs
    template<typename T, int dim, bool is_rational>
    class nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> : public surface<nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>>
    {
    public:
        static constexpr int point_size = is_rational ? dim + 1 : dim;
        static constexpr int dimension = dim;
        using Type = T;
    private:
        
        Eigen::VectorX<T> m_u_knots_vector;
        Eigen::VectorX<T> m_v_knots_vector;
        int m_v_degree;  // = m_v_knots_vector.size() - m_control_points.size() - 1
        int m_u_degree;  // = m_u_knots_vector.size() - m_control_points[0].size() - 1
        //m_control_points[i](, j)为P_ij
        Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> m_control_points;

    public:

        nurbs_surface() = default;
        nurbs_surface(const Eigen::VectorX<T> &u_knots_vector, const Eigen::VectorX<T> &v_knots_vector,
        const Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> &control_points) :
        m_u_knots_vector(u_knots_vector), m_v_knots_vector(v_knots_vector), m_control_points(control_points)
        {
            m_v_degree = m_v_knots_vector.size() - m_control_points.rows() - 1;
            m_u_degree = m_u_knots_vector.size() - m_control_points[0].cols() - 1;
        }

        nurbs_surface(const Eigen::VectorX<T> &u_knots_vector, const Eigen::VectorX<T> &v_knots_vector,
        const Eigen::MatrixX<Eigen::Vector<T, point_size>> &control_points) :
        m_u_knots_vector(u_knots_vector), m_v_knots_vector(v_knots_vector)
        {
            int rows = control_points.rows();
            int cols = control_points.cols();
            m_control_points.resize(rows);
            for (int row_index = 0; row_index < rows; ++row_index)
            {
                m_control_points[row_index].resize(point_size, cols);
                for (int col_index = 0; col_index < cols; ++col_index)
                {
                    m_control_points[row_index].col(col_index) = control_points(row_index, col_index);
                }
            }
            m_v_degree = m_v_knots_vector.size() - m_control_points.rows() - 1;
            m_u_degree = m_u_knots_vector.size() - m_control_points[0].cols() - 1;
        }

        ENUM_NURBS get_uv_knots_end(std::array<T, 2> &u_knots_end, std::array<T, 2> &v_knots_end) const
        {
            u_knots_end[0] = m_u_knots_vector[0];
            u_knots_end[1] = m_u_knots_vector[m_u_knots_vector.size() - 1];
            v_knots_end[0] = m_v_knots_vector[0];
            v_knots_end[1] = m_v_knots_vector[m_v_knots_vector.size() - 1];
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS get_uv_box(Box<T, 2> &uv_box) const
        {
            uv_box.Min[0] = m_u_knots_vector[0];
            uv_box.Min[1] = m_v_knots_vector[0];
            uv_box.Max[0] = m_u_knots_vector[m_u_knots_vector.size() - 1];
            uv_box.Max[1] = m_v_knots_vector[m_v_knots_vector.size() - 1];
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        nurbs_surface(const nurbs_surface &surf)
        {
            m_u_degree = surf.get_u_degree();
            m_v_degree = surf.get_v_degree();
            m_control_points = surf.get_control_points();
            m_u_knots_vector = surf.get_u_knots();
            m_v_knots_vector = surf.get_v_knots();
        }

        int get_u_degree() const 
        {
            return m_u_degree;
        }

        int get_v_degree() const
        {
            return m_v_degree;
        }

        ENUM_NURBS reverse_uv()
        {
            int rows = m_control_points.rows();
            int cols = m_control_points[0].cols();
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points;
            new_control_points.resize(cols);
            for (int col_index = 0; col_index < cols; ++col_index)
            {
                new_control_points[col_index].resize(point_size, rows);
                for (int row_index = 0; row_index < rows; ++row_index)
                {
                    new_control_points[col_index].col(row_index) = m_control_points[row_index].col(col_index);
                }
            }
            m_u_knots_vector.swap(m_v_knots_vector);
            std::swap(m_u_degree, m_v_degree);
            m_control_points = new_control_points;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        Eigen::VectorX<T> get_u_knots() const { return m_u_knots_vector; }

        Eigen::VectorX<T> get_v_knots() const { return m_v_knots_vector; }
        
        ENUM_NURBS set_u_degree(int degree)
        {
            m_u_degree = degree;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS set_v_degree(int degree)
        {
            m_v_degree = degree;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        

        ENUM_NURBS get_u_different_knots(std::vector<T> &vec) const
        {
            return  get_all_defference_knots(m_u_degree, m_u_knots_vector, vec);
        }

        ENUM_NURBS get_v_different_knots(std::vector<T> &vec) const
        {
            return  get_all_defference_knots(m_v_degree, m_v_knots_vector, vec);
        }

        ENUM_NURBS set_uv_degree(int u_degree, int v_degree)
        {
            m_u_degree = u_degree;
            m_v_degree = v_degree;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS set_uv_knots(const Eigen::VectorX<T> &u_knots, const Eigen::VectorX<T> &v_knots)
        {
            m_u_knots_vector = u_knots;
            m_v_knots_vector = v_knots;
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS set_control_points(const Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> &points)
        {
            m_control_points = points;
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> get_nonhomo_control_points() const
        {
            int v_control_points_count = m_control_points.rows();
            int u_control_point_count = m_control_points[0].cols();
            Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> reuslt(v_control_points_count);
            if constexpr(is_rational == false)
            {
                reuslt = m_control_points;
                return reuslt;
            }
            else
            {

                for (int v_index = 0; v_index < v_control_points_count; ++v_index)
                {
                    reuslt[v_index].resize(dim, u_control_point_count);
                    for (int u_index = 0; u_index < u_control_point_count; ++u_index)
                    {
                        reuslt[v_index].col(u_index) = m_control_points[v_index].template block<dim, 1>(0, u_index) / m_control_points[v_index](dim, u_index);
                    }
                }
            }
            return reuslt; 
            
        }


        ENUM_NURBS get_control_points_row(int row, Eigen::Matrix<T, point_size, Eigen::Dynamic> &points) const
        {
            if (row > m_control_points.rows() || row < 0)
                return ENUM_NURBS::NURBS_ERROR;
            points = m_control_points[row];
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS get_control_points_col(int col, Eigen::Matrix<T, point_size, Eigen::Dynamic> &points) const
        {
            int cols = m_control_points[0].cols();
            if (col > cols || col < 0)
                return ENUM_NURBS::NURBS_ERROR;
            int rows = m_control_points.rows();
            points.resize(point_size, rows);
            for (int index = 0; index < rows; ++index)
            {
                points.col(index) = m_control_points[index].col(col);
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> get_control_points() const
        {
            return m_control_points;
        }

        /// @brief 计算曲面上的点
        /// @param u u向参数
        /// @param v v向参数
        /// @param point 曲面上的点
        /// @return ENUM_NURBS错误码
        ENUM_NURBS point_on_surface(T u, T v, Eigen::Vector<T, dim> &point) const
        {
            int uspan = -1, vspan = -1;
            find_span<T>(u, m_u_degree, m_u_knots_vector, uspan);
            find_span<T>(v, m_v_degree,  m_v_knots_vector, vspan);
            Eigen::VectorX<T> nu(m_u_degree + 1);
            Eigen::VectorX<T> nv(m_v_degree + 1);
            basis_functions<T>(uspan, u, m_u_degree, m_u_knots_vector, nu);
            basis_functions<T>(vspan, v, m_v_degree,  m_v_knots_vector, nv);

            Eigen::Matrix<T, point_size, Eigen::Dynamic> temp;
            temp.resize(point_size, m_v_degree + 1);
            for (int index = 0; index <= m_v_degree; ++index)
            {
                temp.col(index) = m_control_points[vspan + index - m_v_degree].block(0, uspan - m_u_degree, point_size, m_u_degree + 1) * nu;
            }
            point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(temp * nv);

            return ENUM_NURBS::NURBS_SUCCESS;
        }

        /// @brief 计算nurbs surface在参数(u, v)处直到n阶偏导数
        /// @tparam n 最高阶数
        /// @param u 参数u
        /// @param v 参数v
        /// @param result out_put_param result(i, j)表示S(u, v)处的u向i次偏导数, 
        /// remark : v向j次偏导数;当(i + j) > m_v_degree + m_u_degree || i > m_u_degree ||j > m_v_degree 时, 原则上来说是非法的, 此处赋值为零向量
        /// @return ENUM_NURBS错误码
        template<int n>
        ENUM_NURBS derivative_on_surface(T u, T v, Eigen::Matrix<Eigen::Vector<T, dim>, n + 1, n + 1> &result) const
        {
            Eigen::Matrix<Eigen::Vector<T, point_size>, n + 1, n + 1> ration_result;
            int du = std::min(n, m_u_degree);
            int dv = std::min(n, m_v_degree);
            Eigen::Vector<T, point_size> zero_vector;
            zero_vector.setConstant(0.0);
            for (int k = m_u_degree + 1; k <= n; ++k)
            {
                int count = n - k;
                for (int l = 0; l <= count; ++l)
                    ration_result(k, l) = zero_vector;
            }

            for (int l = m_v_degree + 1; l <= n; ++l)
            {
                int count = n - l;
                for (int k = 0; k <= count; ++k)
                    ration_result(k, l) = zero_vector;
            }

            int uspan = -1;
            int vspan = -1;
            find_span<T>(u, m_u_degree, m_u_knots_vector, uspan);
            find_span<T>(v, m_v_degree, m_v_knots_vector, vspan);
            Eigen::MatrixX<T> nu, nv;
            ders_basis_funs<T>(uspan, du, m_u_degree, u, m_u_knots_vector, nu);
            ders_basis_funs<T>(vspan, dv, m_v_degree, v, m_v_knots_vector, nv);
            Eigen::Matrix<T, point_size,  Eigen::Dynamic> temps;
            temps.resize(point_size, m_v_degree + 1);
            for (int k = 0; k <= du; ++k)
            {
                temps.setConstant(0.0);
                for (int s = 0; s <= m_v_degree; ++s)
                {
                    temps.col(s) = m_control_points[vspan - m_v_degree + s].block(0, uspan - m_u_degree, point_size, m_u_degree + 1) * nu.col(k);
                }
                int dd = std::min(n - k, dv);
                for (int l = 0; l <= dd; ++l)
                {
                    ration_result(k, l) = temps * nv.col(l);
                }
            }
            result = project_derivs_point<T, is_rational, point_size, n>::project_point_to_euclidean_space(ration_result);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 计算nurbs surface在参数(u, v)处直到n阶偏导数
        /// @param n 最高阶数
        /// @param u 参数u
        /// @param v 参数v
        /// @param result out_put_param result(i, j)表示S(u, v)处的u向i次偏导数, 
        /// remark : v向j次偏导数;当(i + j) > m_v_degree + m_u_degree || i > m_u_degree ||j > m_v_degree 时, 原则上来说是非法的, 此处赋值为零向量
        /// @return ENUM_NURBS错误码
        ENUM_NURBS derivative_on_surface(int n, T u, T v, Eigen::MatrixX<Eigen::Vector<T, dim>> &result) const
        {
            Eigen::MatrixX<Eigen::Vector<T, point_size>> ration_result(n + 1, n + 1);
            int du = std::min(n, m_u_degree);
            int dv = std::min(n, m_v_degree);
            Eigen::Vector<T, point_size> zero_vector;
            zero_vector.setConstant(0.0);
            for (int k = m_u_degree + 1; k <= n; ++k)
            {
                int count = n - k;
                for (int l = 0; l <= count; ++l)
                    ration_result(k, l) = zero_vector;
            }

            for (int l = m_v_degree + 1; l <= n; ++l)
            {
                int count = n - l;
                for (int k = 0; k <= count; ++k)
                    ration_result(k, l) = zero_vector;
            }

            int uspan = -1;
            int vspan = -1;
            find_span<T>(u, m_u_degree, m_u_knots_vector, uspan);
            find_span<T>(v, m_v_degree, m_v_knots_vector, vspan);
            Eigen::MatrixX<T> nu, nv;
            ders_basis_funs<T>(uspan, du, m_u_degree, u, m_u_knots_vector, nu);
            ders_basis_funs<T>(vspan, dv, m_v_degree, v, m_v_knots_vector, nv);
            Eigen::Matrix<T, point_size,  Eigen::Dynamic> temps;
            temps.resize(point_size, m_v_degree + 1);
            for (int k = 0; k <= du; ++k)
            {
                temps.setConstant(0.0);
                for (int s = 0; s <= m_v_degree; ++s)
                {
                    temps.col(s) = m_control_points[vspan - m_v_degree + s].block(0, uspan - m_u_degree, point_size, m_u_degree + 1) * nu.col(k);
                }
                int dd = std::min(n - k, dv);
                for (int l = 0; l <= dd; ++l)
                {
                    ration_result(k, l) = temps * nv.col(l);
                }
            }
            result = project_derivs_point<T, is_rational, point_size, -1>::project_point_to_euclidean_space(ration_result);
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 计算nurbs surface在参数(u, v)处直到n阶偏导数
        /// @param n1 u向求导最高阶数
        /// @param n2 v向求导最高阶数
        /// @param u 参数u
        /// @param v 参数v
        /// @param result out_put_param result(i, j)表示S(u, v)处的u向i次偏导数, 
        /// remark : v向j次偏导数;当(i + j) > m_v_degree + m_u_degree || i > m_u_degree ||j > m_v_degree 时, 原则上来说是非法的, 此处赋值为零向量
        /// @return ENUM_NURBS错误码
        ENUM_NURBS derivative_on_surface(int n1, int n2, T u, T v, Eigen::MatrixX<Eigen::Vector<T, dim>> &result) const
        {
            Eigen::MatrixX<Eigen::Vector<T, point_size>> ration_result(n1 + 1, n2 + 1);

            int uspan = -1;
            int vspan = -1;
            find_span<T>(u, m_u_degree, m_u_knots_vector, uspan);
            find_span<T>(v, m_v_degree, m_v_knots_vector, vspan);
            Eigen::MatrixX<T> nu, nv;
            ders_basis_funs<T>(uspan, n1, m_u_degree, u, m_u_knots_vector, nu);
            ders_basis_funs<T>(vspan, n2, m_v_degree, v, m_v_knots_vector, nv);
            Eigen::Matrix<T, point_size,  Eigen::Dynamic> temps;
            temps.resize(point_size, m_v_degree + 1);
            for (int k = 0; k <= n1; ++k)
            {
                temps.setConstant(0.0);
                for (int s = 0; s <= m_v_degree; ++s)
                {
                    temps.col(s) = m_control_points[vspan - m_v_degree + s].block(0, uspan - m_u_degree, point_size, m_u_degree + 1) * nu.col(k);
                }
                for (int l = 0; l <= n2; ++l)
                {
                    ration_result(k, l) = temps * nv.col(l);
                }
            }
            result = project_derivs_point<T, is_rational, point_size, -1>::project_point_to_euclidean_space(ration_result);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        /// @brief 计算nurbs surface在参数(u, v)处直到n阶偏导数
        /// @tparam n1 u向求导最高阶数
        /// @tparam n2 v向求导最高阶数
        /// @param u 参数u
        /// @param v 参数v
        /// @param result out_put_param result(i, j)表示S(u, v)处的u向i次偏导数, 
        /// remark : v向j次偏导数;当(i + j) > m_v_degree + m_u_degree || i > m_u_degree ||j > m_v_degree 时, 原则上来说是非法的, 此处赋值为零向量
        /// @return ENUM_NURBS错误码
        template<int n1, int n2>
        ENUM_NURBS derivative_on_surface(T u, T v, Eigen::Matrix<Eigen::Vector<T, dim>, n1 + 1, n2 + 1> &result) const
        {
            Eigen::Matrix<Eigen::Vector<T, point_size>, n1 + 1, n2 + 1> ration_result;

            int uspan = -1;
            int vspan = -1;
            find_span<T>(u, m_u_degree, m_u_knots_vector, uspan);
            find_span<T>(v, m_v_degree, m_v_knots_vector, vspan);
            Eigen::Matrix<T, Eigen::Dynamic, n1 + 1> nu;
            Eigen::Matrix<T, Eigen::Dynamic, n2 + 1> nv;
            ders_basis_funs<T, n1>(uspan, m_u_degree, u, m_u_knots_vector, nu);
            ders_basis_funs<T, n2>(vspan, m_v_degree, v, m_v_knots_vector, nv);
            Eigen::Matrix<T, point_size,  Eigen::Dynamic> temps;
            temps.resize(point_size, m_v_degree + 1);
            for (int k = 0; k <= n1; ++k)
            {
                temps.setConstant(0.0);
                for (int s = 0; s <= m_v_degree; ++s)
                {
                    temps.col(s) = m_control_points[vspan - m_v_degree + s].block(0, uspan - m_u_degree, point_size, m_u_degree + 1) * nu.col(k);
                }
                for (int l = 0; l <= n2; ++l)
                {
                    ration_result(k, l) = temps * nv.col(l);
                }
            }
            result = project_derivs_point<T, is_rational, point_size, -1>::project_point_to_euclidean_space(ration_result);
            return ENUM_NURBS::NURBS_SUCCESS;
        }



        ENUM_NURBS tangent_vector_u(double u, double v,  Eigen::Vector<T, dim> &tangent_vector) const
        {
            Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> result;
            derivative_on_surface<1>(u, v, result);
            tangent_vector = result(1, 0);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS tangent_vector_v(double u, double v,  Eigen::Vector<T, dim> &tangent_vector) const
        {
            Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> result;
            derivative_on_surface<1>(u, v, result);
            tangent_vector = result(0, 1);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        /// @brief 将nurbs曲面的一个方向插入节点
        /// @param uv 插入的节点
        /// @param r 插入节点的次数
        /// @param direction 插入节点的方向
        /// @return ENUM_NURBS错误码
        ENUM_NURBS surface_knots_insert(T uv, int r, ENUM_DIRECTION direction)
        {
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            
            int rows = m_control_points.rows();
            int cols = m_control_points[0].cols();
            int knots_size = m_u_knots_vector.size();
            // if (std::abs(uv - m_u_knots_vector[knots_size - 1]) < KNOTS_VECTOR_EPS)
            if (uv == m_u_knots_vector[knots_size - 1])
                return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
            // if (std::abs(uv - m_u_knots_vector[0]) < KNOTS_VECTOR_EPS)
            if (uv == m_u_knots_vector[0])
                return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
            Eigen::Matrix<T, point_size, Eigen::Dynamic> new_control_points;
            int span = -1;
            if(find_span<T>(uv, m_u_degree, m_u_knots_vector, span) != ENUM_NURBS::NURBS_SUCCESS)
                return ENUM_NURBS::NUBRS_CANNOT_INSERT_KNOTS;
            int current_index = span;
            int repeat = 0;
            while (current_index >= 0)
            {
                // if (std::abs(m_u_knots_vector[current_index--] - uv) < KNOTS_VECTOR_EPS)
                if (m_u_knots_vector[current_index--] == uv)
                {
                    ++repeat;
                    continue;
                }
                break;
            }

            for (int row_index = 0; row_index < rows; ++row_index)
            {
                auto control_points = m_control_points[row_index];
                curve_knots_insert<T, point_size>(r, cols, m_u_degree,m_u_knots_vector, control_points, 
                    uv, span, repeat, new_control_points);
                m_control_points[row_index] = new_control_points;
            }

            Eigen::VectorX<T> new_knots_vector(knots_size + r);
            new_knots_vector.block(0, 0, span + 1, 1) = m_u_knots_vector.block(0, 0, span + 1, 1);
            for (int i = 1; i <= r; ++i)
                new_knots_vector[span + i] = uv;
            new_knots_vector.block(span + 1 + r, 0, knots_size - span - 1, 1) = m_u_knots_vector.block(span + 1, 0, knots_size - span - 1, 1);
            m_u_knots_vector = new_knots_vector;
            
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            return ENUM_NURBS::NURBS_SUCCESS;
        
        }

        ENUM_NURBS refine_knots_vector(const Eigen::VectorX<T> &insert_knots, ENUM_DIRECTION direction)
        {
            if (insert_knots.size() == 0)
                return ENUM_NURBS::NURBS_SUCCESS;
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            int degree = direction == ENUM_DIRECTION::V_DIRECTION ? m_v_degree : m_u_degree;
            int rows = m_control_points.rows();
            // int cols = m_control_points[0].cols();
            // int knots_size = m_u_knots_vector.size();
            Eigen::Matrix<T, point_size, Eigen::Dynamic> new_control_points;
            Eigen::VectorX<T> new_knots_vector;
            int knots_vector_size = m_u_knots_vector.size();
            for (int row_index = 0; row_index < rows; ++row_index)
            {
                auto control_points = m_control_points[row_index];
                refine_knots_vector_curve<T, point_size>(knots_vector_size, degree, m_u_knots_vector, control_points,
                    insert_knots, new_knots_vector, new_control_points);
                m_control_points[row_index] = new_control_points;
            }
            m_u_knots_vector = new_knots_vector;
            
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            return ENUM_NURBS::NURBS_SUCCESS;
        
        }

        ENUM_NURBS get_box(Box<T, dim> &box) const
        {
            int v_points_count = m_control_points.rows();
            int u_points_count = m_control_points[0].cols();
            if constexpr (is_rational == false)
            {
                Eigen::Vector<T, dim> min = m_control_points[0].col(0);
                Eigen::Vector<T, dim> max = min;
                for (int v_index = 0; v_index < v_points_count; ++v_index)
                {
                    for (int u_index = 0; u_index < u_points_count; ++u_index)
                    {
                        Eigen::Vector<T, dim> current_point = m_control_points[v_index].col(u_index);
                        for (int index = 0; index < dim; ++index)
                        {
                            if (min[index] > current_point[index])
                                min[index] = current_point[index];
                            if (max[index] < current_point[index])
                                max[index] = current_point[index];
                        }
                    }
                }
                box.Min = min;
                box.Max = max;
            }
            else
            {
                Eigen::Vector<T, dim> min = m_control_points[0].template block<dim, 1>(0, 0) / m_control_points[0](dim, 0);
                Eigen::Vector<T, dim> max = min;
                for (int v_index = 0; v_index < v_points_count; ++v_index)
                {
                    for (int u_index = 0; u_index < u_points_count; ++u_index)
                    {
                        Eigen::Vector<T, dim> current_point = m_control_points[v_index].template block<dim, 1>(0, u_index) / m_control_points[v_index](dim, u_index);
                        for (int index = 0; index < dim; ++index)
                        {
                            if (min[index] > current_point[index])
                                min[index] = current_point[index];
                            if (max[index] < current_point[index])
                                max[index] = current_point[index];
                        }
                    }
                }
                box.Min = min;
                box.Max = max;
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        ENUM_NURBS sub_divide(Box<T, 2> &uv_box, nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &sub_nurbs) const
        {
            sub_nurbs.set_control_points(m_control_points);
            sub_nurbs.set_uv_degree(m_u_degree, m_v_degree);
            sub_nurbs.set_uv_knots(m_u_knots_vector, m_v_knots_vector);
            return sub_nurbs.sub_divide(uv_box);
            // nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> temp_nurbs(m_u_knots_vector, m_v_knots_vector, m_control_points);
            // int u_begin_index;
            // int u_begin_mul = konts_multiple<T>(uv_box.Min[0],  m_u_knots_vector, m_u_degree, u_begin_index);
            // int u_begin_insert_num = std::max(0, m_u_degree - u_begin_mul);
            // int u_end_index;
            // int u_end_mul = konts_multiple<T>(uv_box.Max[0], m_u_knots_vector, m_u_degree, u_end_index);
            // int u_end_insert_num = std::max(m_u_degree - u_end_mul, 0);
            // Eigen::VectorX<T> insert_knots(u_begin_insert_num + u_end_insert_num);
            // insert_knots.block(0, 0, u_begin_insert_num, 1).setConstant(uv_box.Min[0]);
            // insert_knots.block(u_begin_insert_num, 0, u_end_insert_num, 1).setConstant(uv_box.Max[0]);
            // temp_nurbs.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);

            // int v_begin_index;
            // int v_begin_mul = konts_multiple<T>(uv_box.Min[1], m_v_knots_vector, m_v_degree, v_begin_index);
            // int v_begin_insert_num = std::max(m_v_degree - v_begin_mul, 0);
            // int v_end_index;
            // int v_end_mul = konts_multiple<T>(uv_box.Max[1], m_v_knots_vector, m_v_degree, v_end_index);
            // int v_end_insert_num = std::max(m_v_degree - v_end_mul, 0);
            // insert_knots.resize(v_begin_insert_num + v_end_insert_num);
            // insert_knots.block(0, 0, v_begin_insert_num, 1).setConstant(uv_box.Min[1]);
            // insert_knots.block(v_begin_insert_num, 0, v_end_insert_num, 1).setConstant(uv_box.Max[1]);
            // temp_nurbs.refine_knots_vector(insert_knots, ENUM_DIRECTION::V_DIRECTION);


            // int u_new_knots_count = 2 * m_u_degree + 2 + u_end_index - u_begin_index - u_begin_mul;
            // int v_new_knots_count = 2 * m_v_degree + 2 + v_end_index - v_begin_index - v_begin_mul;
            // Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points(v_new_knots_count - m_v_degree - 1);
            // Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> points = temp_nurbs.get_control_points();
            // Eigen::VectorX<T> u_new_knots(u_new_knots_count);
            // Eigen::VectorX<T> temp_u_knots = temp_nurbs.get_u_knots();
            // u_new_knots[0] = uv_box.Min[0];
            // u_new_knots[u_new_knots_count - 1] = m_u_knots_vector[u_end_index];
            // u_new_knots.block(1, 0, u_new_knots_count - 2, 1) = temp_u_knots.block(std::max(u_begin_index, 1), 0, u_new_knots_count - 2, 1);

            // Eigen::VectorX<T> v_new_knots(v_new_knots_count);
            // Eigen::VectorX<T> temp_v_knots = temp_nurbs.get_v_knots();
            // v_new_knots[0] = m_v_knots_vector[v_begin_index];
            // v_new_knots[v_new_knots_count - 1] = m_v_knots_vector[v_end_index];
            // v_new_knots.block(1, 0, v_new_knots_count - 2, 1) = temp_v_knots.block(std::max(v_begin_index, 1), 0, v_new_knots_count - 2, 1);
            // for (int v_index = 0; v_index < v_new_knots_count - m_v_degree - 1; ++v_index)
            // {
            //     new_control_points[v_index] = points[v_index + v_begin_index].block(0, u_begin_index, point_size, u_new_knots_count - m_u_degree - 1);
            // }
            // sub_nurbs.set_control_points(new_control_points);
            // sub_nurbs.set_uv_degree(m_u_degree, m_v_degree);
            // sub_nurbs.set_uv_knots(u_new_knots, v_new_knots);
            // return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS sub_divide(Box<T, 2> &uv_box)
        {
            nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> temp_nurbs(m_u_knots_vector, m_v_knots_vector, m_control_points);
            int u_begin_index;
            int u_begin_mul = konts_multiple<T>(uv_box.Min[0],  m_u_knots_vector, m_u_degree, u_begin_index);
            int u_begin_insert_num = std::max(0, m_u_degree - u_begin_mul);
            int u_end_index;
            int u_end_mul = konts_multiple<T>(uv_box.Max[0], m_u_knots_vector, m_u_degree, u_end_index);
            int u_end_insert_num = std::max(m_u_degree - u_end_mul, 0);
            Eigen::VectorX<T> insert_knots(u_begin_insert_num + u_end_insert_num);
            insert_knots.block(0, 0, u_begin_insert_num, 1).setConstant(uv_box.Min[0]);
            insert_knots.block(u_begin_insert_num, 0, u_end_insert_num, 1).setConstant(uv_box.Max[0]);
            temp_nurbs.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);

            int v_begin_index;
            int v_begin_mul = konts_multiple<T>(uv_box.Min[1], m_v_knots_vector, m_v_degree, v_begin_index);
            int v_begin_insert_num = std::max(m_v_degree - v_begin_mul, 0);
            int v_end_index;
            int v_end_mul = konts_multiple<T>(uv_box.Max[1], m_v_knots_vector, m_v_degree, v_end_index);
            int v_end_insert_num = std::max(m_v_degree - v_end_mul, 0);
            insert_knots.resize(v_begin_insert_num + v_end_insert_num);
            insert_knots.block(0, 0, v_begin_insert_num, 1).setConstant(uv_box.Min[1]);
            insert_knots.block(v_begin_insert_num, 0, v_end_insert_num, 1).setConstant(uv_box.Max[1]);
            temp_nurbs.refine_knots_vector(insert_knots, ENUM_DIRECTION::V_DIRECTION);


            int u_new_knots_count = 2 * m_u_degree + 2 + u_end_index - u_begin_index - u_begin_mul;
            int v_new_knots_count = 2 * m_v_degree + 2 + v_end_index - v_begin_index - v_begin_mul;
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points(v_new_knots_count - m_v_degree - 1);
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> points = temp_nurbs.get_control_points();
            Eigen::VectorX<T> u_new_knots(u_new_knots_count);
            Eigen::VectorX<T> temp_u_knots = temp_nurbs.get_u_knots();
            u_new_knots[0] = uv_box.Min[0];
            u_new_knots[u_new_knots_count - 1] = uv_box.Max[0];
            u_new_knots.block(1, 0, u_new_knots_count - 2, 1) = temp_u_knots.block(std::max(u_begin_index, 1), 0, u_new_knots_count - 2, 1);

            Eigen::VectorX<T> v_new_knots(v_new_knots_count);
            Eigen::VectorX<T> temp_v_knots = temp_nurbs.get_v_knots();
            v_new_knots[0] = uv_box.Min[1];// m_v_knots_vector[v_begin_index];
            v_new_knots[v_new_knots_count - 1] = uv_box.Max[1];//  m_v_knots_vector[v_end_index];
            v_new_knots.block(1, 0, v_new_knots_count - 2, 1) = temp_v_knots.block(std::max(v_begin_index, 1), 0, v_new_knots_count - 2, 1);
            
            int v_start = std::max(v_begin_index - 1, 0);
            int u_start = std::max(u_begin_index - 1, 0);
            for (int v_index = 0; v_index < v_new_knots_count - m_v_degree - 1; ++v_index)
            {
                // std::cout << points[v_index + v_start] << std::endl;
                new_control_points[v_index] = points[v_index + v_start].block(0, u_start, point_size, u_new_knots_count - m_u_degree - 1);
            }
            m_control_points = new_control_points;
            m_u_knots_vector = u_new_knots;
            m_v_knots_vector = v_new_knots;
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 将nurbs曲线分解成bezier曲面
        /// @param bezier_curves 分解的bezier曲线, 内存用户释放
        /// @return ENUM_NURBS错误码
        ENUM_NURBS decompose_to_bezier(Eigen::MatrixX<bezier_surface<T, dim,-1, -1, is_rational> *> &bezier_surfaces) const
        {
            //先u向加细成bezier
            Eigen::VectorX<T> new_knots_vector;
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points;
            int u_interval_count = -1, v_interval_count = -1;
            find_interval_segment_count<T>(m_u_degree, m_u_knots_vector, u_interval_count);
            find_interval_segment_count<T>(m_v_degree, m_v_knots_vector, v_interval_count);
            new_control_points.resize(u_interval_count);
            int rows_count = m_control_points.rows() - 1;

            //(B_ij)(k,l)表示第(i, j)段beizer曲面的第(k, l)个控制点
            Eigen::MatrixX<Eigen::MatrixX<Eigen::Vector<T, point_size>>> bezier_surface_control_points(u_interval_count, v_interval_count);
            for (int i = 0; i < u_interval_count; ++i)
            {
                for (int j = 0; j < v_interval_count; ++j)
                {
                    bezier_surface_control_points(i,j).resize(m_v_degree + 1, m_u_degree + 1);
                }
            }
            
            //T[i].col(j)表示第i列曲线的第j个控制点
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> temp_control_points;
            int cols_count = u_interval_count * m_u_degree;
            temp_control_points.resize(cols_count + 1);
            int rsize = m_v_knots_vector.size() - m_v_degree - 1;
            for (int i = 0; i <= cols_count; ++i)
            {
                temp_control_points[i].resize(point_size, rsize);
            }

            for (int index = 0; index <= rows_count; ++index)
            {
                decompose_curve_to_bezier<T, point_size>(m_u_degree, u_interval_count, m_u_knots_vector, m_control_points[index], new_knots_vector, new_control_points);
                for (int i = 0; i < u_interval_count; ++i)
                {
                    for (int k = 0; k < m_u_degree; ++k)
                    {       
                        temp_control_points[i * m_u_degree + k].col(index) = new_control_points[i].col(k);
                    }
                }
                temp_control_points[cols_count].col(index) = new_control_points[u_interval_count - 1].col(m_u_degree);           
            }

            //在加细v向
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> v_new_control_points;
            v_new_control_points.resize(v_interval_count);
            for (int index = 0; index < cols_count; ++index)
            {
                int br = index % m_v_degree;
                int r = index / m_v_degree;
                decompose_curve_to_bezier<T, point_size>(m_v_degree, v_interval_count, m_v_knots_vector, temp_control_points[index], new_knots_vector, v_new_control_points);
                
                for (int i = 0; i < v_interval_count; ++i)
                {
                    for (int k = 0; k <= m_v_degree; ++k)
                    {
                        bezier_surface_control_points(r, i)(k, br) = v_new_control_points[i].col(k);
                    }
                    if (r != 0 && br == 0)
                    {
                        for (int k = 0; k <= m_v_degree; ++k)
                        {
                            bezier_surface_control_points(r - 1, i)(k, m_u_degree) = v_new_control_points[i].col(k);
                        }    
                    }
                }        
            }
            decompose_curve_to_bezier<T, point_size>(m_v_degree, v_interval_count, m_v_knots_vector, temp_control_points[cols_count], new_knots_vector, v_new_control_points);
            for (int i = 0; i < v_interval_count; ++i)
            {
                for (int k = 0; k <= m_v_degree; ++k)
                {
                    bezier_surface_control_points(u_interval_count - 1, i)(k, m_u_degree) = v_new_control_points[i].col(k);
                }
            }

            bezier_surfaces.resize(u_interval_count, v_interval_count);
            for (int i = 0; i < u_interval_count; ++i)
            {
                for (int j = 0; j < v_interval_count; ++j)
                {
                    bezier_surface<T, dim, -1, -1, is_rational> *bs = new bezier_surface<T, dim, -1, -1, is_rational>(bezier_surface_control_points(i, j));
                    bezier_surfaces(i, j) = bs;
                }
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 将nurbs曲线分解成bezier曲面
        /// @param bezier_curves 分解的bezier曲线, 内存用户释放
        /// @return ENUM_NURBS错误码
        ENUM_NURBS decompose_to_bezier(Eigen::MatrixX<nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> *> &bezier_surfaces) const
        {
            //先u向加细成bezier
            Eigen::VectorX<T> new_knots_vector;
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points;
            int u_interval_count = -1, v_interval_count = -1;
            find_interval_segment_count<T>(m_u_degree, m_u_knots_vector, u_interval_count);
            find_interval_segment_count<T>(m_v_degree, m_v_knots_vector, v_interval_count);
            new_control_points.resize(u_interval_count);
            int rows_count = m_control_points.rows() - 1;

            //(B_ij)(k,l)表示第(i, j)段beizer曲面的第(k, l)个控制点
            Eigen::MatrixX<Eigen::MatrixX<Eigen::Vector<T, point_size>>> bezier_surface_control_points(u_interval_count, v_interval_count);
            for (int i = 0; i < u_interval_count; ++i)
            {
                for (int j = 0; j < v_interval_count; ++j)
                {
                    bezier_surface_control_points(i,j).resize(m_v_degree + 1, m_u_degree + 1);
                }
            }
            
            //T[i].col(j)表示第i列曲线的第j个控制点
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> temp_control_points;
            int cols_count = u_interval_count * m_u_degree;
            temp_control_points.resize(cols_count + 1);
            int rsize = m_v_knots_vector.size() - m_v_degree - 1;
            for (int i = 0; i <= cols_count; ++i)
            {
                temp_control_points[i].resize(point_size, rsize);
            }

            for (int index = 0; index <= rows_count; ++index)
            {
                decompose_curve_to_bezier<T, point_size>(m_u_degree, u_interval_count, m_u_knots_vector, m_control_points[index], new_knots_vector, new_control_points);
                for (int i = 0; i < u_interval_count; ++i)
                {
                    for (int k = 0; k < m_u_degree; ++k)
                    {       
                        temp_control_points[i * m_u_degree + k].col(index) = new_control_points[i].col(k);
                    }
                }
                temp_control_points[cols_count].col(index) = new_control_points[u_interval_count - 1].col(m_u_degree);           
            }

            //在加细v向
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> v_new_control_points;
            v_new_control_points.resize(v_interval_count);
            for (int index = 0; index < cols_count; ++index)
            {
                int br = index % m_v_degree;
                int r = index / m_v_degree;
                decompose_curve_to_bezier<T, point_size>(m_v_degree, v_interval_count, m_v_knots_vector, temp_control_points[index], new_knots_vector, v_new_control_points);
                
                for (int i = 0; i < v_interval_count; ++i)
                {
                    for (int k = 0; k <= m_v_degree; ++k)
                    {
                        bezier_surface_control_points(r, i)(k, br) = v_new_control_points[i].col(k);
                    }
                    if (r != 0 && br == 0)
                    {
                        for (int k = 0; k <= m_v_degree; ++k)
                        {
                            bezier_surface_control_points(r - 1, i)(k, m_u_degree) = v_new_control_points[i].col(k);
                        }    
                    }
                }        
            }
            decompose_curve_to_bezier<T, point_size>(m_v_degree, v_interval_count, m_v_knots_vector, temp_control_points[cols_count], new_knots_vector, v_new_control_points);
            for (int i = 0; i < v_interval_count; ++i)
            {
                for (int k = 0; k <= m_v_degree; ++k)
                {
                    bezier_surface_control_points(u_interval_count - 1, i)(k, m_u_degree) = v_new_control_points[i].col(k);
                }
            }

            std::vector<Eigen::VectorX<T>> new_u_knots_vectors(u_interval_count);
            int current_u_knots_index = 0;
            int u_interval_index = 0;
            int u_knots_count = m_u_knots_vector.size();
            for (int index = 0; index < u_knots_count; ++index)
            {
                if (m_u_knots_vector[current_u_knots_index] != m_u_knots_vector[index])
                {
                    new_u_knots_vectors[u_interval_index].resize(2 * m_u_degree + 2);
                    new_u_knots_vectors[u_interval_index].block(0, 0, m_u_degree + 1, 1).setConstant(m_u_knots_vector[current_u_knots_index]);
                    new_u_knots_vectors[u_interval_index].block(m_u_degree + 1, 0, m_u_degree + 1, 1).setConstant(m_u_knots_vector[index]);
                    ++u_interval_index;
                    current_u_knots_index = index;
                }
            }

            std::vector<Eigen::VectorX<T>> new_v_knots_vectors(v_interval_count);
            int current_v_knots_index = 0;
            int v_interval_index = 0;
            int v_knots_count = m_v_knots_vector.size();
            for (int index = 0; index < v_knots_count; ++index)
            {
                if (m_v_knots_vector[current_v_knots_index] != m_v_knots_vector[index])
                {
                    new_v_knots_vectors[v_interval_index].resize(2 * m_v_degree + 2);
                    new_v_knots_vectors[v_interval_index].block(0, 0, m_v_degree + 1, 1).setConstant(m_v_knots_vector[current_v_knots_index]);
                    new_v_knots_vectors[v_interval_index].block(m_v_degree + 1, 0, m_v_degree + 1, 1).setConstant(m_v_knots_vector[index]);
                    ++v_interval_index;
                    current_v_knots_index = index;
                }
            }

            bezier_surfaces.resize(u_interval_count, v_interval_count);
            for (int i = 0; i < u_interval_count; ++i)
            {
                for (int j = 0; j < v_interval_count; ++j)
                {
                    nurbs_surface<T, dim,-1, -1, -1, -1, is_rational> *bs = new nurbs_surface<T, dim,-1, -1, -1, -1, is_rational>(new_u_knots_vectors[i], new_v_knots_vectors[j], bezier_surface_control_points(i, j));
                    bezier_surfaces(i, j) = bs;
                }
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
        ENUM_NURBS remove_knots(T u, int count, ENUM_DIRECTION direction, int &time, T error = DEFAULT_ERROR)
        {
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            
            int knots_size = m_u_knots_vector.size();
            if (u == m_u_knots_vector[0] || u == m_u_knots_vector[knots_size - 1])
            {
                time = 0;
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            int span = -1;
            find_span<T>(u, m_u_degree, m_u_knots_vector, span);
            if (m_u_knots_vector[span] != u)
                return ENUM_NURBS::NURBS_ERROR;
            
            int repeat = 1;
            int current = span - 1;
            while (current >= 0 && m_u_knots_vector[current] == u)
            {
                repeat += 1;
                current -= 1;
            }
            count = std::min(repeat, count);
            Eigen::VectorX<T> new_knots_vector;
            int rows = m_control_points.rows();
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points = m_control_points;
            for (time = 0; time < count; ++time)
            {
                Eigen::VectorX<T> knots_vector;
                for(int row = 0; row < rows; ++row)
                {
                    knots_vector = m_u_knots_vector;
                    int ture_or_false = -1;
                    Eigen::Matrix<T, point_size, Eigen::Dynamic> points = new_control_points[row];
                    remove_curve_knots<T, point_size, is_rational>(span, 1, m_u_degree, repeat - time, knots_vector, points, ture_or_false, error);
                    if (ture_or_false != 1)
                    {
                        if (direction == ENUM_DIRECTION::V_DIRECTION)
                            reverse_uv();
                        return ENUM_NURBS::NURBS_SUCCESS;
                    }
                    new_control_points[row] = points;
                }
                m_u_knots_vector = knots_vector;
                m_control_points = new_control_points;
            }
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 升阶
        /// @param t 提升的阶数
        /// @param direction 提升的方向
        /// @return ENUM_NURBS错误码
        ENUM_NURBS degree_elevate(int t, ENUM_DIRECTION direction)
        {
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            Eigen::VectorX<T> new_knots_vector;
            int count = m_control_points.rows();
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points(count);
            ENUM_NURBS flag;
            for (int index = 0; index < count; ++index)
            {
                flag = degree_elevate_curve<T, point_size>(t, m_u_degree, m_u_knots_vector, m_control_points[index], new_knots_vector, new_control_points[index]);
                if (flag != ENUM_NURBS::NURBS_SUCCESS)
                {
                    if (direction == ENUM_DIRECTION::V_DIRECTION)
                        reverse_uv();
                    return flag;
                }
            }

            m_u_degree += t;
            m_u_knots_vector = new_knots_vector;
            m_control_points = new_control_points;

            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 将nurbs surface的某个方向降低一阶
        /// @param direction 降阶的方向
        /// @param error 误差
        /// @return ENUM_NURBS错误码
        ENUM_NURBS degree_reduce(ENUM_DIRECTION direction , T error = DEFAULT_ERROR)
        {
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();

            int count = m_control_points.rows();
            ENUM_NURBS flag;
            Eigen::VectorX<T> new_knots_vector;
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points(count);
            for (int index = 0; index < count; ++index)
            {
                flag = degree_reduce_curve<T, point_size, is_rational>(m_u_degree, m_u_knots_vector, m_control_points[index], new_knots_vector, new_control_points[index], error);
                if (flag != ENUM_NURBS::NURBS_SUCCESS)
                {
                    if (direction == ENUM_DIRECTION::V_DIRECTION)
                        reverse_uv();
                    return flag;
                }
            }
            m_control_points = new_control_points;
            m_u_knots_vector = new_knots_vector;
            m_u_degree -= 1;
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();

            return ENUM_NURBS::NURBS_SUCCESS;
        }


        bool is_u_closed() const
        {
            int v_count = m_control_points.rows();
            int u_count = m_control_points[0].cols();
            for (int index = 0; index < v_count; ++index)
            {
                Eigen::Vector<T, dim> start_point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(m_control_points[index].col(0));
                Eigen::Vector<T, dim> end_point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(m_control_points[index].col(u_count - 1));
                if ((start_point - end_point).squaredNorm() > DEFAULT_ERROR * DEFAULT_ERROR)
                    return false;
            }
            return true;
        }

        bool is_v_closed() const
        {
            int v_count = m_control_points.rows();
            int u_count = m_control_points[0].cols();
            for (int index = 0; index < u_count; ++index)
            {
                Eigen::Vector<T, dim> start_point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(m_control_points[0].col(index));
                Eigen::Vector<T, dim> end_point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(m_control_points[v_count - 1].col(index));

                if ((start_point - end_point).squaredNorm() > DEFAULT_ERROR * DEFAULT_ERROR)
                    return false;
            }
            return true;
        }

        //TODO: 写个多线程, 多线程每次结果可能都会不同
        ENUM_NURBS find_nearst_point_on_surface(const Eigen::Vector<T, dim> &point, T &u, T &v, Eigen::Vector<T, dim> &nearst_point) const
        {
            // TODO : 包围盒加速
            T u_min = m_u_knots_vector[0];
            int u_knots_size = m_u_knots_vector.size();
            T u_max = m_u_knots_vector[u_knots_size - 1];

            T v_min = m_v_knots_vector[0];
            int v_knots_size = m_v_knots_vector.size();
            T v_max = m_v_knots_vector[v_knots_size - 1];

            //将节点分成(max - min + 1) * 100份；以后需要优化
            // TODO: 曲率大的地方分的多一些，曲率小的地方分少一些
            int u_step_count = (u_max - u_min + 1) * 10;
            int v_step_count = (v_max - v_min + 1) * 10;


            T u_step = (u_max - u_min) / u_step_count;
            T v_step = (v_max - v_min) / v_step_count;
            nearst_point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(m_control_points[0].col(0));
            T min_length = (point - nearst_point).squaredNorm();
            u = m_u_knots_vector[0];
            v = m_v_knots_vector[0];
            Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> ders_vec;
            std::vector<T> u_param(u_step_count + 1);
            std::vector<T> v_param(v_step_count + 1);
            u_param[0] = m_u_knots_vector[0];
            for (int u_index = 1; u_index < u_step_count; ++u_index)
                u_param[u_index] = u_param[u_index - 1] + u_step;
            u_param[u_step_count] = u_max;

            v_param[0] = m_v_knots_vector[0];
            for (int v_index = 1; v_index < v_step_count; ++v_index)
                v_param[v_index] = v_param[v_index - 1] + v_step;
            v_param[v_step_count] = v_max;
            
            bool is_u_closed_flag = is_u_closed();
            bool is_v_closed_flag = is_v_closed();

            for (int u_index = 0; u_index <= u_step_count; ++u_index)
            {
                for (int v_index = 0; v_index <= v_step_count; ++v_index)
                {
                    T current_u = u_param[u_index];
                    T current_v = v_param[v_index];
                    T distance;
                    
                    for (int loop_index = 0; loop_index < MAX_SURFACE_ITERATE_DEEP; ++loop_index)
                    {
                        derivative_on_surface<2>(current_u, current_v, ders_vec);
                        Eigen::Vector<T, dim> dist_vec = ders_vec(0, 0) - point;
                        distance = dist_vec.squaredNorm();
                        if (distance < DEFAULT_ERROR * DEFAULT_ERROR)
                        {
                            u = current_u;
                            v = current_v;
                            nearst_point = ders_vec(0, 0);
                            return ENUM_NURBS::NURBS_SUCCESS;
                        }
                        if (distance < min_length)
                        {
                            min_length = distance;
                            u = current_u;
                            v = current_v;
                            nearst_point = ders_vec(0, 0);
                        }

                        T cos_angle_1 = ders_vec(1, 0).dot(dist_vec);
                        T cos_angle_2 = ders_vec(0, 1).dot(dist_vec);
                        T tanget_vec_square_len_1 = ders_vec(1, 0).squaredNorm();
                        T tanget_vec_square_len_2 = ders_vec(0, 1).squaredNorm();
                        bool flag1 = cos_angle_1 * cos_angle_1 < tanget_vec_square_len_1 * distance * ANGLE_ERROR * ANGLE_ERROR;
                        bool flag2 = cos_angle_2 * cos_angle_2 < tanget_vec_square_len_2 * distance * ANGLE_ERROR * ANGLE_ERROR;
                        if ((flag1 == true) && (flag2 == true))
                        {
                            break;
                        }
                        Eigen::Matrix<T, 2, 2> J;
                        J(0, 0) = tanget_vec_square_len_1 + dist_vec.dot(ders_vec(2, 0));
                        J(1, 0) = ders_vec(1, 0).dot(ders_vec(0, 1)) + dist_vec.dot(ders_vec(1, 1));
                        J(0, 1) = J(1, 0);
                        J(1, 1) = tanget_vec_square_len_2 + dist_vec.dot(ders_vec(0, 2));

                        Eigen::Vector<T, 2> K;
                        K[0] = -1 * dist_vec.dot(ders_vec(1, 0));
                        K[1] = -1 * dist_vec.dot(ders_vec(0, 1));

                        double det = J.determinant();
                        if (std::abs(det) < DEFAULT_ERROR)
                        {
                            //退化矩阵, break
                            break;
                        }
                        Eigen::JacobiSVD<Eigen::MatrixX<T>,  Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(J);
                        Eigen::Vector<T, 2> delta_param = matSvd.solve(K);
                        
                        // TODO ：根据曲线弧长和区间长度放缩delta_param
                        T next_u = current_u + delta_param[0];
                        T next_v = current_v + delta_param[1];
                        
                        if (is_u_closed_flag)
                        {
                            if (next_u < u_min)
                                next_u = u_max - (u_min - next_u);
                            else if (next_u > u_max)
                                next_u = u_min + (next_u - u_max);
                        }
                        else
                        {
                            if (next_u < u_min)
                                next_u = u_min;
                            else if (next_u > u_max)
                                next_u = u_max;
                        }

                        if (is_v_closed_flag)
                        {
                            if (next_v < v_min)
                                next_v = v_max - (v_min - next_v);
                            else if (next_v > v_max)
                                next_v = v_min + (next_v - v_max);
                        }
                        else
                        {
                            if (next_v < v_min)
                                next_v = v_min;
                            else if (next_v > v_max)
                                next_v = v_max;
                        }
                        if (((next_u - current_u) * ders_vec(1, 0) + (next_v - next_u) * ders_vec(0, 1)).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR)
                        {
                            break;
                        }
                        current_u = next_u;
                        current_v = next_v;
                    }
                }
                
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        //反求参数空间的切向量,没有测试
        ENUM_NURBS surface_tangent_vector_inversion(const Eigen::Vector<T, dim> &tangent_on_surface, 
            const Eigen::Vector<T, dim> &point ,Eigen::Vector<T, 2> &tangent_on_param_space) const
        {
            T u, v;
            Eigen::Vector<T, dim> nearst_point;
            find_nearst_point_on_surface(point, u, v, nearst_point);
            if ((nearst_point - point).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR)
            {
                return NURBS_POINT_IS_NOT_ON_CURVE;
            }
            Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> ders_vec;
            derivative_on_surface<1>(u, v, ders_vec);

            Eigen::Matrix<T, 2, 2> M;
            Eigen::Vector<T, 2> vec;
            M(0, 0) = ders_vec(1, 0).squaredNorm();
            M(0, 1) = ders_vec(1, 0).dot(ders_vec(0, 1));
            M(1, 0) = M(0, 1);
            M(1, 1) = ders_vec(0, 1).squaredNorm();
            vec[0] = ders_vec(1, 0).dot(point);
            vec[1] = ders_vec(0, 1).dot(point);
            Eigen::JacobiSVD<Eigen::MatrixX<T>,  Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(M);
            tangent_on_param_space = matSvd.solve(vec);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS surface_reparameter(const nurbs_curve<T, 1, true, -1, -1> &reparameter_function, ENUM_DIRECTION direction,
            nurbs_surface<T, dim, -1, -1, -1, -1, true> &new_nurbs_surface) const
        {
            nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> *copy_surface = new nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>(*this);
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                copy_surface->reverse_uv();
            Eigen::VectorX<Eigen::Matrix<T, point_size, -1>> copy_nurbs_surface_control_points = copy_surface->get_control_points();
            int v_control_points_count = copy_nurbs_surface_control_points.rows();
            Eigen::VectorX<Eigen::Matrix<T, dim + 1, Eigen::Dynamic>> new_control_points(v_control_points_count);
            int new_u_degree = copy_surface->get_u_degree() * reparameter_function.get_degree();
            Eigen::VectorX<T> new_u_knots_vector;
            for (int v_index = 0; v_index < v_control_points_count; ++v_index)
            {
                Eigen::Matrix<T, dim + 1, Eigen::Dynamic> rational_control_points;
                to_ratioanl_contrl_points<T, is_rational, point_size>::convert(copy_nurbs_surface_control_points[v_index], rational_control_points);
                nurbs_curve<T, dim, true, -1, -1> *rational_nurbs = new nurbs_curve<T, dim, true, -1, -1>(copy_surface->get_u_knots(), rational_control_points);
                nurbs_curve<T, dim, true, -1, -1> *reparamter_nurbs = new nurbs_curve<T, dim, true, -1, -1>();
                rational_nurbs->curve_reparameter(reparameter_function, *reparamter_nurbs);
                new_control_points[v_index] = reparamter_nurbs->get_control_points();
                if (v_index == 0)
                    new_u_knots_vector = reparamter_nurbs->get_knots_vector();
                delete rational_nurbs;
                delete reparamter_nurbs;
            }
            new_nurbs_surface.set_u_degree(new_u_degree);
            new_nurbs_surface.set_control_points(new_control_points);
            new_nurbs_surface.set_v_degree(copy_surface->get_v_degree());
            new_nurbs_surface.set_uv_knots(new_u_knots_vector, copy_surface->get_v_knots());

            delete copy_surface;
            if (direction == ENUM_DIRECTION::V_DIRECTION)
            {
                new_nurbs_surface.reverse_uv();
            }

            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS surface_reparameter(const nurbs_curve<T, 1, false, -1, -1> &reparameter_function, ENUM_DIRECTION direction,
            nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &new_nurbs_surface) const
        {
            nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> *copy_surface = new nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>(*this);
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                copy_surface.reverse_uv();

            Eigen::VectorX<Eigen::Matrix<T, point_size, -1>> copy_nurbs_surface_control_points = copy_surface->get_control_points();
            int v_control_points_count = copy_nurbs_surface_control_points.rows();
            Eigen::VectorX<Eigen::Matrix<T, dim + 1, Eigen::Dynamic>> new_control_points(v_control_points_count);
            int new_u_degree = copy_surface->get_u_degree() * reparameter_function.get_degree();
            Eigen::VectorX<T> new_u_knots_vector;
            Eigen::VectorX<T> copy_surface_u_knots = copy_surface->get_u_knots();
            for (int v_index = 0; v_index < v_control_points_count; ++v_index)
            {
                nurbs_curve<T, dim, is_rational, -1, -1> *rational_nurbs = new nurbs_curve<T, dim, true, -1, -1>(copy_surface_u_knots, copy_nurbs_surface_control_points[v_index]);
                nurbs_curve<T, dim, is_rational, -1, -1> *reparamter_nurbs = new nurbs_curve<T, dim, true, -1, -1>();
                rational_nurbs->curve_reparameter(reparameter_function, *reparamter_nurbs);
                new_control_points[v_index] = reparamter_nurbs->get_control_points();
                if (v_index == 0)
                    new_u_knots_vector = reparamter_nurbs->get_knots_vector();
                delete rational_nurbs;
                delete reparamter_nurbs;
            }
            new_nurbs_surface.set_u_degree(new_u_degree);
            new_nurbs_surface.set_control_points(new_control_points);
            new_nurbs_surface.set_v_degree(copy_surface->get_v_degree());
            new_nurbs_surface.set_uv_knots(new_u_knots_vector, copy_surface->get_v_knots());
            delete copy_surface;

            if (direction == ENUM_DIRECTION::V_DIRECTION)
            {
                new_nurbs_surface.reverse_uv();
            }

            return ENUM_NURBS::NURBS_SUCCESS;
        }


        // u = f(s) = alpha * s + beta
        ENUM_NURBS surface_reparameter_with_linear_function(T alpha, T beta, ENUM_DIRECTION direction, 
            nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &new_nurbs_surface) const
        {
            new_nurbs_surface.set_control_points(m_control_points);
            Eigen::VectorX<T> old_knots_vector = direction == ENUM_DIRECTION::U_DIRECTION ? m_u_knots_vector : m_v_knots_vector;
            int old_knots_vector_size = old_knots_vector.size();
            Eigen::VectorX<T> new_knots_vector(old_knots_vector_size);
            new_knots_vector[0] = (old_knots_vector[0] - beta) / alpha;
            int current_index = 0;
            for (int index = 1; index < old_knots_vector_size; ++index)
            {
                if (old_knots_vector[index] != old_knots_vector[current_index])
                {
                    new_knots_vector[index] = (old_knots_vector[index] - beta) / alpha;
                    current_index = index;
                }
                else
                {
                    new_knots_vector[index] = new_knots_vector[current_index];
                }
            }
            new_nurbs_surface.set_uv_degree(m_u_degree, m_v_degree);
            if (direction == ENUM_DIRECTION::U_DIRECTION)
            {
                new_nurbs_surface.set_uv_knots(new_knots_vector, m_v_knots_vector);
            }
            else
            {
                new_nurbs_surface.set_uv_knots(m_u_knots_vector, m_v_knots_vector);
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        // s = g(u) = (alpha * u + beta) / (gamma * u + delta)
        // alpha = a11, beta = a12, gamma = a21, dalta = a22
        ENUM_NURBS surface_reparameter_with_linear_function(Eigen::Matrix<T, 2, 2> reparameter_function, ENUM_DIRECTION direction,
            nurbs_surface<T, dim, -1, -1, -1, -1, true> &new_nurbs_surface) const
        {
            if (reparameter_function.determinant() < DEFAULT_ERROR)
                return ENUM_NURBS::NURBS_ERROR;

            nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> *copy_surface = new nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>(*this);

            if (direction == ENUM_DIRECTION::V_DIRECTION)
                copy_surface->reverse_uv();

            int v_control_points_count = copy_surface->m_control_points.rows();
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> copy_surface_control_points = copy_surface->get_control_points();
            Eigen::VectorX<Eigen::Matrix<T, dim + 1, Eigen::Dynamic>> new_control_points(v_control_points_count);
            Eigen::VectorX<T> new_u_knots_vector;
            for (int v_index = 0; v_index < v_control_points_count; ++v_index)
            {
                Eigen::Matrix<T, dim + 1, Eigen::Dynamic> rational_control_points;
                to_ratioanl_contrl_points<T, is_rational, point_size>::convert(copy_surface_control_points[v_index], rational_control_points);
                nurbs_curve<T, dim, true, -1, -1> *rational_nurbs = new nurbs_curve<T, dim, true, -1, -1>(m_u_knots_vector, rational_control_points);
                nurbs_curve<T, dim, true, -1, -1> *reparamter_nurbs = new nurbs_curve<T, dim, true, -1, -1>();
                rational_nurbs->curve_reparameter_with_linear_function(reparameter_function, *reparamter_nurbs);
                new_control_points[v_index] = reparamter_nurbs->get_control_points();
                if (v_index == 0)
                    new_u_knots_vector = reparamter_nurbs->get_knots_vector();
                delete rational_nurbs;
                delete reparamter_nurbs;
            }
            new_nurbs_surface.set_u_degree(copy_surface->get_u_degree());
            new_nurbs_surface.set_control_points(new_control_points);
            new_nurbs_surface.set_v_degree(copy_surface->get_v_degree());
            new_nurbs_surface.set_uv_knots(new_u_knots_vector, copy_surface->get_v_knots());
            delete copy_surface;
            if (direction == ENUM_DIRECTION::V_DIRECTION)
            {
                new_nurbs_surface.reverse_uv();
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS surface_reverse(ENUM_DIRECTION direction)
        {
            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();
            
            int v_control_point_count = m_control_points.rows();
            int u_control_point_count = m_control_points[0].cols();
            int u_knots_vector_count = m_u_knots_vector.size();
            Eigen::VectorX<T> new_knots_vector(u_knots_vector_count);
            T interval_begin = m_u_knots_vector[0];
            T interval_end = m_u_knots_vector[u_knots_vector_count - 1];
            Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> new_control_points(v_control_point_count);
            for (int index = 0; index < u_knots_vector_count; ++index)
            {
                new_knots_vector[index] = interval_end + interval_begin - m_u_knots_vector[u_knots_vector_count - index - 1];
            }
            for (int v_index = 0; v_index < v_control_point_count; ++v_index)
            {
                for (int index = 0; index < u_control_point_count; ++index)
                {
                    new_control_points[index].col(index) = m_control_points[index].col(u_control_point_count - index - 1);
                }
            }
            m_control_points = new_control_points;
            m_u_knots_vector = new_knots_vector;

            if (direction == ENUM_DIRECTION::V_DIRECTION)
                reverse_uv();

            return ENUM_NURBS::NURBS_SUCCESS;
        }


        ENUM_NURBS scale_surface(const Eigen::Vector<T, dim> &center, const Eigen::Matrix<T, dim, dim> &mat, const Eigen::Vector<T, dim> &scale_factory)
        {
            int rows = m_control_points.rows();
            int cols = m_control_points[0].cols();
            for (int row_index = 0; row_index < rows; ++row_index)
            {
                for (int cols_index = 0; cols_index < cols; ++cols_index)
                {
                    Eigen::Vector<T, point_size> vec = m_control_points[row_index].col(cols_index);
                    sacle_point<T, dim, is_rational>(vec, mat, scale_factory, center);
                    m_control_points[row_index].col(cols_index) = vec;
                }
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        /// @brief 提取等参线(还没有测试)
        /// @tparam direction  direction = ENUM_DIRECTION::U_DIRECTION表示等u线; direction = ENUM_DIRECTION::V_DIRECTION表示等v线
        /// @param param S(direction = param)
        /// @return 错误码
        template<int direction> 
        ENUM_NURBS get_isoparameter_curve(T param, nurbs_curve<T, dim, is_rational, -1, -1> &nurbs) const
        {
            if constexpr (direction == ENUM_DIRECTION::U_DIRECTION)
            {
                int span = -1;
                find_span<T>(param, m_u_knots_vector, span);
                int cols = m_control_points.rows();
                Eigen::Matrix<T, point_size, Eigen::Dynamic> new_control_points(point_size, cols);
                Eigen::VectorX<T> nu(m_u_degree + 1);
                basis_functions<T>(span, param, m_u_knots_vector, nu);
                for (int index = 0; index < cols; ++index)
                { 
                    new_control_points.col(index) = m_control_points[index].block(0, span - m_u_degree, point_size, m_u_degree + 1) * nu;
                }
                nurbs.set_control_points(new_control_points);
                nurbs.set_knots_vector(m_v_knots_vector);
                nurbs.set_degree(m_v_degree);
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            else if constexpr (direction == ENUM_DIRECTION::V_DIRECTION)
            {
                //reverse_uv
                int rows = m_control_points.rows();
                int cols = m_control_points[0].cols();
                Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> reverse_control_points;
                reverse_control_points.resize(cols);
                for (int col_index = 0; col_index < cols; ++col_index)
                {
                    reverse_control_points[col_index].resize(point_size, rows);
                    for (int row_index = 0; row_index < rows; ++row_index)
                    {
                        reverse_control_points[col_index].col(row_index) = m_control_points[row_index].col(col_index);
                    }
                }
                
                int span = -1;
                find_span<T>(param, m_v_knots_vector, span);
                Eigen::Matrix<T, point_size, Eigen::Dynamic> new_control_points(point_size, cols);
                Eigen::VectorX<T> nv(m_v_degree + 1);
                basis_functions<T>(span, param, m_v_knots_vector, nv);
                for (int index = 0; index < cols; ++index)
                {
                    new_control_points.col(index) = new_control_points[index].block(0, span - m_v_degree, point_size, m_v_degree + 1) * nv;
                }
                nurbs.set_control_points(new_control_points);
                nurbs.set_knots_vector(m_u_knots_vector);
                nurbs.set_degree(m_u_degree);
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            //else
            return ENUM_NURBS::NURBS_ERROR;
        }

        /// @brief 给定误差, 在误差范围内消去所有可以消去的节点; 此函数仅在非有理nurbs曲面是合法的(稍微修改下一误差即可在有理曲面下也是合法的)
        /// @param params 消去的曲线在参数params上误差小于给定的误差; 即(new_nurbs(params[i]) - nurbs(params[i])).norm() < E;
        /// @param error 返回的各个点的误差
        /// @param new_nurbs 节点消去后新的nurbs
        /// @param E 最大误差
        /// @return 错误码
        ENUM_NURBS remove_knots_bound_surface(const std::vector<std::array<T, 2>> &params, std::vector<T> &errors, 
                nurbs_surface<T, dim, -1, -1, -1, -1, false> &new_nurbs, T E = DEFAULT_ERROR) const
        {
            static_assert(is_rational ==false, "rational b spline surface can not use this function");
            if constexpr(is_rational == true)
            {
                return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
            }

            std::vector<T> u_different_knots;
            std::vector<int> u_multiple;
            std::vector<T> v_different_knots;
            std::vector<int> v_multiple;
            get_different_knots(m_u_knots_vector, m_u_degree, u_different_knots, u_multiple);
            get_different_knots(m_v_knots_vector, m_v_degree, v_different_knots, v_multiple);
            //获得每一行和每一列的Br
            int u_different_knots_count = u_different_knots.size();
            int v_different_knots_count = v_different_knots.size();
            std::vector<std::vector<T>> u_Brs;
            std::vector<std::vector<T>> v_Brs;
            int u_points_count = m_control_points[0].cols();
            int v_points_count = m_control_points.rows();            
            u_Brs.reserve(v_points_count);
            v_Brs.reserve(u_points_count);

            for (int v_index = 0; v_index < v_points_count; ++v_index)
            {
                int current_index = m_u_degree;
                std::vector<T> Brs(u_different_knots_count - 2, 0.0);
                //计算内节点Br的值
                for (int u_index = 1; u_index < u_different_knots_count - 1; ++u_index)
                {
                    int r = current_index + u_multiple[u_index];
                    T Br;
                    get_removal_bnd_curve<T, dim>(m_control_points[v_index], m_u_knots_vector, m_u_degree, 
                                            u_different_knots[u_index], r, u_multiple[u_index], Br);
                    Brs[u_index - 1] = Br;
                    current_index = r;
                }
                u_Brs.push_back(std::move(Brs));
            }

            for (int u_index = 0; u_index < u_points_count; ++u_index)
            {
                int current_index = m_v_degree;
                std::vector<T> Brs(v_different_knots_count - 2, 0.0);
                Eigen::Matrix<T, dim, Eigen::Dynamic> v_control_ponts(dim, v_points_count);
                for (int index = 0; index < v_points_count; ++index)
                {
                    v_control_ponts.col(index) = m_control_points[index].col(u_index);
                }

                //计算内节点Br的值
                for (int v_index = 1; v_index < v_different_knots_count - 1; ++v_index)
                {
                    int r = current_index + v_multiple[v_index];
                    T Br;
                    get_removal_bnd_curve<T, dim>(v_control_ponts, m_v_knots_vector, m_v_degree, 
                                            v_different_knots[v_index], r, v_multiple[v_index], Br);
                    Brs[v_index - 1] = Br;
                    current_index = r;
                }
                v_Brs.push_back(std::move(Brs));
            }
        
            //对每个基函数, 计算相关参数params的下标
            std::list<std::vector<int>> u_knots_params_relation;
            std::list<std::vector<int>> v_knots_params_relation;

            int params_count = static_cast<int> (params.size());
            for (int u_index = 0; u_index < u_points_count; ++u_index)
            {
                T u_min = m_u_knots_vector[u_index];
                T u_max = m_u_knots_vector[u_index + m_u_degree + 1];
                std::vector<int> indexs;
                for (int index = 0; index < params_count; ++index)
                {
                    if (params[index][0] >= u_min && params[index][0] < u_max)
                        indexs.push_back(index);
                }
                u_knots_params_relation.push_back(std::move(indexs));
            }
            for (int v_index = 0; v_index < v_points_count; ++v_index)
            {
                T v_min = m_v_knots_vector[v_index];
                T v_max = m_v_knots_vector[v_index + m_v_degree + 1];
                std::vector<int> indexs;
                for (int index = 0; index < params_count; ++index)
                {
                    if (params[index][1] >= v_min && params[index][1] < v_max)
                        indexs.push_back(index);
                }
                v_knots_params_relation.push_back(std::move(indexs));
            }
            errors.resize(params_count, 0.0);
            Eigen::VectorX<T> u_new_knots = m_u_knots_vector;
            Eigen::VectorX<T> v_new_knots = m_v_knots_vector;
            Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> new_control_points = m_control_points;
            while (true)
            {
                //找到Brs和最小的u
                T u_min_sum = -1.0;
                int u_min_index = -1;
                for (int u_index = 0; u_index < u_different_knots_count - 2; ++u_index)
                {
                    bool flag = false;
                    T tem_min_sum = 0.0;
                    for (int v_index = 0; v_index < v_points_count; ++v_index)
                    {
                        flag = true;
                        if (u_Brs[v_index][u_index] == -1.0)
                        {
                            tem_min_sum = -1.0;
                            break;
                        }
                        tem_min_sum += u_Brs[v_index][u_index];
                    }
                    if (tem_min_sum != -1.0)
                    {
                        if (u_min_sum > tem_min_sum || (u_min_sum == -1.0 && flag == true))
                        {
                            u_min_sum = tem_min_sum;
                            u_min_index = u_index;
                        }
                    }
                }
                //找到Brs和最小的v
                T v_min_sum = -1.0;
                int v_min_index = -1;
                for (int v_index = 0; v_index < v_different_knots_count - 2; ++v_index)
                {
                    bool flag = false;
                    T tem_min_sum = 0.0;
                    for (int u_index = 0; u_index < u_points_count; ++u_index)
                    {
                        flag = true;
                        if (v_Brs[u_index][v_index] == -1.0)
                        {
                            tem_min_sum = -1.0;
                            break;
                        }
                        tem_min_sum += v_Brs[u_index][v_index];
                    }
                    if (tem_min_sum != -1.0)
                    {
                        if (v_min_sum > tem_min_sum || (v_min_sum == -1.0 && flag == true))
                        {
                            v_min_sum = tem_min_sum;
                            v_min_index = v_index;
                        }
                    }
                }

                //消去最小的u或者v
                if (u_min_sum == -1.0 && v_min_sum == -1.0)
                    break;
                ENUM_DIRECTION direction;
                if (u_min_sum == -1.0)
                    direction = ENUM_DIRECTION::V_DIRECTION;
                else if (v_min_index == -1.0)
                    direction = ENUM_DIRECTION::U_DIRECTION;
                else
                {
                    T min_sum = std::min(u_min_sum, v_min_sum);
                    direction = min_sum == u_min_sum ? ENUM_DIRECTION::U_DIRECTION : ENUM_DIRECTION::V_DIRECTION;
                }

                int new_v_points_count = new_control_points.rows();
                int new_u_points_count = new_control_points[0].cols();

                std::vector<T> temp_error = errors;
                if (direction == ENUM_DIRECTION::U_DIRECTION)
                {
                    if (u_Brs.size() == 0)
                        break;
                    Eigen::Matrix<T, 1, Eigen::Dynamic> temp_control_points(1, new_v_points_count);
                    for (int i = 0; i < new_v_points_count; ++i)
                        temp_control_points(0, i) = u_Brs[i][u_min_index];
                    nurbs_curve<T, 1, false, -1, -1> temp_nurbs(v_new_knots, temp_control_points);
                    int r = std::accumulate(u_multiple.begin(), u_multiple.begin() + u_min_index + 2, -1);
                    
                    int &mul = u_multiple[u_min_index + 1];
                    bool flag = true;
                    if ((m_u_degree + mul) % 2 == 0)
                    {
                        int k = (m_u_degree + mul) / 2;
                        int i = r - k;
                        auto current_params_it = u_knots_params_relation.begin();
                        for (int index = 1; index <= i; ++index)
                            ++current_params_it;
                        const std::vector<int> &current_params = *current_params_it;//u_knots_params_relation[u_min_index];
                        for (const int index : current_params)
                        {                      
                            int temp_index = -1;
                            find_span<T>(params[index][1], m_v_degree, v_new_knots, temp_index);
                            for (int j = temp_index - m_v_degree; j <= temp_index; ++j)
                            {
                                if (temp_control_points(0, j) == -1.0)
                                {
                                    for (int l = 1; l < v_different_knots_count - 1; ++l)
                                        u_Brs[l - 1][u_min_index] = -1.0;
                                    flag = false;
                                    break;
                                }
                            }
                            if (flag == false)
                                break;
                            Eigen::Vector<T, 1> right;
                            temp_nurbs.point_on_curve(params[index][1], right);
                            T basis;
                            one_basis_function<T>(i, params[index][0], m_u_degree, u_new_knots, basis);
                            temp_error[index] += basis * right[0];
                        }
                    }
                    else
                    {
                        int k = (mul + m_u_degree + 1) / 2;
                        int i = r - k + 1;
                        auto current_params_it = u_knots_params_relation.begin();
                        for (int index = 1; index <= i; ++index)
                            ++current_params_it;
                        const std::vector<int> &current_params = *current_params_it;//u_knots_params_relation[u_min_index];
                        T alpha = (u_new_knots[r] - u_new_knots[r - k + 1]) / (u_new_knots[r - k + m_u_degree + 2] - u_new_knots[r - k + 1]);
                        for (const int index : current_params)
                        {
                            int temp_index = -1;
                            find_span<T>(params[index][1], m_v_degree, v_new_knots, temp_index);
                            for (int j = temp_index - m_v_degree; j <= temp_index; ++j)
                            {
                                if (temp_control_points(0, j) == -1.0)
                                {
                                    u_Brs[j][u_min_index] = -1.0;
                                    flag = false;
                                    break;
                                }
                            }
                            if (flag == false)
                                break;
                            Eigen::Vector<T, 1> right;
                            temp_nurbs.point_on_curve(params[index][1], right);
                            T basis;
                            one_basis_function<T>(i, params[index][0], m_u_degree, u_new_knots, basis);
                            temp_error[index] += basis * (1.0 - alpha) * right[0];
                        }
                    }
                    T max_error = *(std::max_element(temp_error.begin(), temp_error.end()));
                    if (max_error < E && flag == true)
                    {
                        errors = temp_error;
                        Eigen::VectorX<T> temp_u_new_knots;
                        for (int k = 0; k < new_v_points_count; ++k)
                        {
                            temp_u_new_knots = u_new_knots;
                            //可能有问题
                            remove_curve_knots_no_check_error<T, dim, false>(r, 1, m_u_degree, mul, temp_u_new_knots, new_control_points[k]);
                        }
                        u_new_knots = temp_u_new_knots;
                        u_points_count = new_control_points[0].cols();
                        auto u_knots_params_relation_it = u_knots_params_relation.begin();
                        for (int k = 0; k < r - m_u_degree - 1; ++k)
                            ++u_knots_params_relation_it;

                        for (int k = r - m_u_degree - 1; k <= r - mul; ++k)
                        {
                            T u_min = u_new_knots[k];
                            T u_max = u_new_knots[k + m_u_degree + 1];
                            std::vector<int> indexs;
                            for (int index = 0; index < params_count; ++index)
                            {
                                if (params[index][0] >= u_min && params[index][0] < u_max)
                                    indexs.push_back(index);
                            }
                            *u_knots_params_relation_it = indexs;
                            ++u_knots_params_relation_it;
                        }
                        u_knots_params_relation.erase(u_knots_params_relation_it);

                        int beg = std::max(r - m_u_degree, m_u_degree + 1);
                        int end = std::min(r + m_u_degree - mul + 1, new_u_points_count - 1);
                        //更新v_Brs
                        for (int u_index = r - m_u_degree; u_index < r - mul; ++u_index)
                        {
                            int current_index = m_v_degree;
                            std::vector<T> Brs(v_different_knots_count - 2, 0.0);
                            //计算内节点Br的值
                            for (int v_index = 1; v_index < v_different_knots_count - 1; ++v_index)
                            {
                                Eigen::Matrix<T, dim, Eigen::Dynamic> points(dim, v_points_count);
                                for (int i = 0; i < v_points_count; ++i)
                                    points.col(i) = new_control_points[i].col(u_index);
                                int r1 =  current_index + v_multiple[v_index];
                                T Br;
                                get_removal_bnd_curve<T, dim>(points, v_new_knots, m_v_degree, 
                                                        v_different_knots[v_index], r1, v_multiple[v_index], Br);
                                Brs[v_index - 1] = Br;
                                current_index = r1;
                            }
                            v_Brs[u_index] = Brs;
                        }
                        v_Brs.erase(v_Brs.begin() + (r - mul));
                        std::vector<int> need_eval_br_index;
                        int current_index = m_u_degree + 1;
                        for (int i = 1; i < u_different_knots_count - 1; ++i)
                        {
                            int next_index = current_index + u_multiple[i];
                            if (current_index >= beg && current_index <= end)
                            {
                                if (i == (u_min_index + 1) && mul == 1)
                                    continue;
                                if (i > (u_min_index + 1) && mul == 1)
                                    need_eval_br_index.push_back(i - 2);
                                else
                                    need_eval_br_index.push_back(i - 1);
                            }
                            current_index = next_index;
                        }
                        bool removed = false;
                        if (mul == 1)
                        {
                            removed = true;
                            u_different_knots.erase(u_different_knots.begin() + u_min_index + 1);
                            u_different_knots_count = u_different_knots.size();
                            u_multiple.erase(u_multiple.begin() + u_min_index + 1);
                        }
                        else
                        {
                            mul -= 1;
                        }

                        //计算Br
                        for (int v_index = 0; v_index < v_points_count; ++v_index)
                        {
                            // std::vector<T> Brs(u_different_knots_count - 2, 0.0);
                            //计算内节点Br的值
                            if (removed == true)
                                u_Brs[v_index].erase(u_Brs[v_index].begin() + u_min_index);
                            for (int u_need_eval_index : need_eval_br_index)
                            {
                                int r = std::accumulate(u_multiple.begin(), u_multiple.begin() + u_need_eval_index + 2, -1);
                                T Br;
                                get_removal_bnd_curve<T, dim>(new_control_points[v_index], u_new_knots, m_u_degree, 
                                                        u_different_knots[u_need_eval_index + 1], r, u_multiple[u_need_eval_index + 1], Br);
                                u_Brs[v_index][u_need_eval_index] = Br;
                            }
                        }
                    }
                    else
                    {
                        for (int index = 0; index < v_points_count; ++index)
                            u_Brs[index][u_min_index] = -1.0;
                    }
                }
                else
                {
                    if (v_Brs.size() == 0)
                        break;
                    Eigen::Matrix<T, 1, Eigen::Dynamic> temp_control_points(1, new_u_points_count);
                    for (int i = 0; i < new_u_points_count; ++i)
                        temp_control_points(0, i) = v_Brs[i][v_min_index];
                    nurbs_curve<T, 1, false, -1, -1> temp_nurbs(u_new_knots, temp_control_points);
                    int r = std::accumulate(v_multiple.begin(), v_multiple.begin() + v_min_index + 2, -1);

                    int &mul = v_multiple[v_min_index + 1];
                    bool flag = true;
                    if ((m_v_degree + mul) % 2 == 0)
                    {
                        int k = (m_v_degree + mul) / 2;
                        int i = r - k;
                        auto current_params_it = v_knots_params_relation.begin();
                        for (int index = 1; index <= i; ++index)
                            ++current_params_it;
                        const std::vector<int> &current_params = *current_params_it;//v_knots_params_relation[u_min_index];
                        for (const int index : current_params)
                        {
                            int temp_index = -1;
                            find_span<T>(params[index][0], m_u_degree, u_new_knots, temp_index);
                            for (int j = temp_index - m_u_degree; j <= temp_index; ++j)
                            {
                                if (temp_control_points(0, j) == -1.0)
                                {
                                    for (int i = 1; i < u_different_knots_count - 1; ++i)
                                        v_Brs[i - 1][v_min_index] = -1.0;
                                    flag = false;
                                    break;
                                }
                            }
                            if (flag == false)
                                break;
                            Eigen::Vector<T, 1> right;
                            temp_nurbs.point_on_curve(params[index][0], right);
                            T basis;
                            one_basis_function<T>(i, params[index][1], m_v_degree, v_new_knots, basis);
                            temp_error[index] += basis * right[0];
                        }
                    }
                    else
                    {
                        int k = (mul + m_v_degree + 1) / 2;
                        int i = r - k + 1;
                        auto current_params_it = v_knots_params_relation.begin();
                        for (int index = 1; index <= i; ++index)
                            ++current_params_it;
                        T alpha = (v_new_knots[r] - v_new_knots[r - k + 1]) / (v_new_knots[r - k + m_v_degree + 2] - v_new_knots[r - k + 1]);
                        const std::vector<int> &current_params = *current_params_it;//v_knots_params_relation[v_min_index];
                        for (const int index : current_params)
                        {
                            int temp_index = -1;
                            find_span<T>(params[index][0], m_u_degree, u_new_knots, temp_index);
                            for (int j = temp_index - m_u_degree; j <= temp_index; ++j)
                            {
                                if (temp_control_points(0, j) == -1.0)
                                {
                                    for (int l = 1; l < u_different_knots_count - 1; ++l)
                                        v_Brs[l - 1][v_min_index] = -1.0;
                                    flag = false;
                                    break;
                                }
                            }
                            if (flag == false)
                                break;
                            Eigen::Vector<T, 1> right;
                            temp_nurbs.point_on_curve(params[index][0], right);
                            T basis;
                            one_basis_function<T>(i, params[index][1], m_v_degree, v_new_knots, basis);
                            // T right_temp;
                            // one_basis_function<T>(i - 1, params[index][1], m_v_degree, v_new_knots, right_temp);
                            // T temp_test = right_temp * right[0] * (1.0 - alpha);
                            temp_error[index] += basis * right[0] * (1.0 - alpha);
                        }
                    }
                    T max_error = *(std::max_element(temp_error.begin(), temp_error.end()));
                    if (max_error < E && flag == true)
                    {
                        errors = temp_error;
                        Eigen::VectorX<T> temp_v_new_knots;
                        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> points(new_v_points_count - 1);
                        for (int l = 0; l < new_v_points_count - 1; ++l)
                            points[l].resize(dim, new_u_points_count);
                        for (int k = 0; k < new_u_points_count; ++k)
                        {
                            temp_v_new_knots = v_new_knots;
                            Eigen::Matrix<T, dim, Eigen::Dynamic> temp_control_points(dim, new_v_points_count);
                            for (int j = 0; j < new_v_points_count; ++j)
                            {
                                temp_control_points.col(j) = new_control_points[j].col(k);
                            }
                            remove_curve_knots_no_check_error<T, dim, false>(r, 1, m_v_degree, mul, temp_v_new_knots, temp_control_points);
                            for (int j = 0; j < new_v_points_count - 1; ++j)
                            {
                                points[j].col(k) = temp_control_points.col(j);
                            }
                        }
                        new_control_points = points;
                        v_new_knots = temp_v_new_knots;
                        v_points_count = new_control_points.rows();
                        auto v_knots_params_relation_it = v_knots_params_relation.begin();
                        for (int k = 0; k < r - m_v_degree - 1; ++k)
                            ++v_knots_params_relation_it;

                        for (int k = r - m_v_degree - 1; k <= r - mul; ++k)
                        {
                            T v_min = v_new_knots[k];
                            T v_max = v_new_knots[k + m_v_degree + 1];
                            std::vector<int> indexs;
                            for (int index = 0; index < params_count; ++index)
                            {
                                if (params[index][1] >= v_min && params[index][1] < v_max)
                                    indexs.push_back(index);
                            }
                            *v_knots_params_relation_it = indexs;
                            ++v_knots_params_relation_it;
                        }
                        v_knots_params_relation.erase(v_knots_params_relation_it);

                        int beg = std::max(r - m_v_degree, m_v_degree + 1);
                        int end = std::min(r + m_v_degree - mul + 1, new_v_points_count - 1);
                        //更新u_Brs
                        for (int v_index = r - m_v_degree; v_index < r - mul; ++v_index)
                        {
                            int current_index = m_u_degree;
                            std::vector<T> Brs(u_different_knots_count - 2, 0.0);
                            //计算内节点Br的值
                            for (int u_index = 1; u_index < u_different_knots_count - 1; ++u_index)
                            {
                                int r1 = current_index + u_multiple[u_index];
                                T Br;
                                get_removal_bnd_curve<T, dim>(new_control_points[v_index], u_new_knots, m_u_degree, 
                                                        u_different_knots[u_index], r1, u_multiple[u_index], Br);
                                Brs[u_index - 1] = Br;
                                current_index = r1;
                            }
                            u_Brs[v_index] = Brs;
                        }
                        u_Brs.erase(u_Brs.begin() + (r - mul));
                        std::vector<int> need_eval_br_index;
                        int current_index = m_v_degree + 1;
                        for (int i = 1; i < v_different_knots_count - 1; ++i)
                        {
                            int next_index = current_index + v_multiple[i];
                            if (current_index >= beg && current_index <= end)
                            {
                                if (i == (v_min_index + 1) && mul == 1)
                                    continue;
                                if (i > (v_min_index + 1) && mul == 1)
                                    need_eval_br_index.push_back(i - 2);
                                else
                                    need_eval_br_index.push_back(i - 1);
                            }
                            current_index = next_index;
                        }
                        bool removed = false;
                        if (mul == 1)
                        {
                            v_different_knots.erase(v_different_knots.begin() + v_min_index + 1);
                            v_different_knots_count = v_different_knots.size();
                            v_multiple.erase(v_multiple.begin() + v_min_index + 1);
                            removed = true;
                        }
                        else
                        {
                            mul -= 1;
                        }

                        //计算Br
                        for (int u_index = 0; u_index < u_points_count; ++u_index)
                        {
                            // std::vector<T> Brs(v_different_knots_count - 2, 0.0);
                            //计算内节点Br的值
                            if (removed == true)
                                v_Brs[u_index].erase(v_Brs[u_index].begin() + v_min_index);
                            Eigen::Matrix<T, dim, Eigen::Dynamic> temp_points(dim, v_points_count);
                            for (int i = 0; i < v_points_count; ++i)
                                temp_points.col(i) = new_control_points[i].col(u_index);
                            for (int v_need_eval_index : need_eval_br_index)
                            {
                                int r = std::accumulate(v_multiple.begin(), v_multiple.begin() + v_need_eval_index + 2, -1);
                                T Br;
                                get_removal_bnd_curve<T, dim>(temp_points, v_new_knots, m_v_degree, 
                                                        v_different_knots[v_need_eval_index + 1], r, v_multiple[v_need_eval_index + 1], Br);
                                v_Brs[u_index][v_need_eval_index] = Br;
                            }
                        }
                    }
                    else
                    {
                        for (int index = 0; index < u_points_count; ++index)
                            v_Brs[index][v_min_index] = -1.0;
                    }
                }
                
            }

            new_nurbs.set_control_points(new_control_points);
            new_nurbs.set_uv_degree(m_u_degree, m_v_degree);
            new_nurbs.set_uv_knots(u_new_knots, v_new_knots);
            return ENUM_NURBS::NURBS_SUCCESS;
        
        }


    };

    template<typename T, int dim, int rows, int cols, int u_degree, int v_degree, bool is_rational>
    struct geo_traits<nurbs_surface<T, dim, rows, cols, u_degree, v_degree, is_rational> >
    {
        static constexpr int point_size = is_rational ? dim + 1 : dim;
        using type = nurbs_surface<T, dim, rows, cols, u_degree, v_degree, is_rational>;
        using point_number_type = T;
        // using dimension = dim;
        using point_type = typename  Eigen::Vector<T, dim> ;
    };
}
