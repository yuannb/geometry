#pragma once
#include "nurbs_tool.h"
#include "bezier_surface.h"

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
    //TODO:
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
        Eigen::MatrixX<Eigen::Matrix<Eigen::Vector<T, point_size>, u_degree + 1, v_degree + 1>> bezier_surface_control_points(u_interval_count, v_interval_count);

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
                    bezier_surface_control_points(r, i)(br, k) = v_new_control_points[i].col(k);
                }
                if (r != 0 && br == 0)
                {
                    for (int k = 0; k <= v_degree; ++k)
                    {
                        bezier_surface_control_points(r - 1, i)(u_degree, k) = v_new_control_points[i].col(k);
                    }    
                }
            }        
        }
        decompose_curve_to_bezier<T, point_size>(v_degree, v_interval_count, m_v_knots_vector, temp_control_points[cols_count], new_knots_vector, v_new_control_points);
        for (int i = 0; i < v_interval_count; ++i)
        {
            for (int k = 0; k <= v_degree; ++k)
            {
                bezier_surface_control_points(u_interval_count - 1, i)(u_degree, k) = v_new_control_points[i].col(k);
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
class nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>
{
private:
    static constexpr int point_size = is_rational ? dim + 1 : dim;
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

    ENUM_NURBS set_control_points(Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> &points)
    {
        m_control_points = points;
        return ENUM_NURBS::NURBS_SUCCESS;
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
        // std::cout << nv << std::endl;
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


    /// @brief 将nurbs曲线分解成bezier曲面
    /// @param bezier_curves 分解的bezier曲线, 内存用户释放
    /// @return ENUM_NURBS错误码
    ENUM_NURBS decompose_to_bezier(Eigen::MatrixX<bezier_surface<T, dim,-1, -1, is_rational> *> &bezier_surfaces)
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
                bezier_surface_control_points(i,j).resize(m_u_degree + 1, m_v_degree + 1);
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
                    bezier_surface_control_points(r, i)(br, k) = v_new_control_points[i].col(k);
                }
                if (r != 0 && br == 0)
                {
                    for (int k = 0; k <= m_v_degree; ++k)
                    {
                        bezier_surface_control_points(r - 1, i)(m_u_degree, k) = v_new_control_points[i].col(k);
                    }    
                }
            }        
        }
        decompose_curve_to_bezier<T, point_size>(m_v_degree, v_interval_count, m_v_knots_vector, temp_control_points[cols_count], new_knots_vector, v_new_control_points);
        for (int i = 0; i < v_interval_count; ++i)
        {
            for (int k = 0; k <= m_v_degree; ++k)
            {
                bezier_surface_control_points(u_interval_count - 1, i)(m_u_degree, k) = v_new_control_points[i].col(k);
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


    bool is_u_closed()
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

    bool is_v_closed()
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
    ENUM_NURBS find_nearst_point_on_surface(const Eigen::Vector<T, dim> &point, T &u, T &v, Eigen::Vector<T, dim> &nearst_point)
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
                    Eigen::JacobiSVD<Eigen::Matrix<T, 2, 2>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(J);
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
        const Eigen::Vector<T, dim> &point ,Eigen::Vector<T, 2> &tangent_on_param_space)
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
        Eigen::JacobiSVD<Eigen::Matrix<T, 2, 2>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(M);
        tangent_on_param_space = matSvd.solve(vec);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

};

