#pragma once
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <array>
#include <vector>
#include <iostream>
#include <fstream>

// #define KNOTS_VECTOR_EPS    1e-8
#define KNOTS_VECTOR_ERROR     1e-8
#define DEFAULT_ERROR     1e-4
#define ANGLE_ERROR       1e-3
#define MAX_ITERATE_DEEP 1e5
#define MAX_SURFACE_ITERATE_DEEP 100
#define MAX_ITERATE_STEP    1e-3
#define INDEX_IS_OUTSIDE_OF_KNOTS_VECTOR    -1

enum ENUM_DIRECTION
{
    U_DIRECTION = 0,
    V_DIRECTION = 1
};

// #define U_DIRECTION = 0
// #define V_DIRECTION = 1

enum ENUM_NURBS
{
    NURBS_ERROR = 0,
    NURBS_SUCCESS = 1,
    NURBS_PARAM_IS_OUT_OF_DOMAIN = 2,
    NURBS_DERIV_DEGREE_IS_NOT_EXIST = 3,
    NUBRS_CANNOT_INSERT_KNOTS = 4,
    NUBRS_WIEGHT_IS_NONPOSITIVE = 5,
    NURBS_POINT_IS_ON_CURVE = 6,
    NURBS_POINT_IS_NOT_ON_CURVE = 7
};

/// @brief 查找节点矢量有几个不一样的数
/// @param T double float int...
/// @param knots_vector 节点矢量
/// @param count out_put_param 节点矢量中不同的数的个数
/// @return ENUM_NURBS错误码
template<typename T>
ENUM_NURBS find_knots_num(const Eigen::VectorX<T> &knots_vector, int &count)
{
    int knots_size = knots_vector.size();
    if (knots_size == 0)
        return ENUM_NURBS::NURBS_ERROR;
    count = 1;
    T current = knots_vector[0];
    for (int index = 1; index < knots_size; ++index)
    {
        if (current != knots_vector[index])
        {
            count += 1;
            current = knots_vector[index];
        }
    }
    return ENUM_NURBS::NURBS_SUCCESS;

}


/// @brief 如果是射影空间, 将射影点转换为参数点(投影最后一个坐标)
/// @tparam T double float int...
/// @tparam flag 是否是射影空间
/// @tparam rows 坐标点分量的个数
template<typename T, bool flag, int rows>
struct project_point
{
    static Eigen::Vector<T, rows - 1> project_point_to_euclidean_space(Eigen::Vector<T, rows> const &point)
    {
        return point.block(0, 0, rows - 1, 1) / point[rows - 1];
    }
};

template<typename T, int rows>
struct project_point<T, false, rows>
{
    static Eigen::Vector<T, rows> project_point_to_euclidean_space(Eigen::Vector<T, rows> const &point)
    {
        return point;
    }
};

template<int n>
Eigen::Matrix<int, n + 1, n + 1> binary_coeff()
{
    Eigen::Matrix<int, n + 1, n + 1 > result;
    result.col(0).setConstant(1);
    for (int i = 1; i <= n; ++i)
    {
        for (int j = 1; j <= i; ++j)
        {
            if (j == i)
                result(i, i) = 1;
            else
                result(i, j) = result(i - 1, j) + result(i - 1, j - 1);
        }
    }
    return result;
}


Eigen::MatrixX<int> binary_coeff(int d)
{
    Eigen::MatrixX<int> result(d, d);
    result.col(0).setConstant(1);
    for (int i = 1; i < d; ++i)
    {
        for (int j = 1; j <= i; ++j)
        {
            if (j == i)
                result(i, i) = 1;
            else
                result(i, j) = result(i - 1, j) + result(i - 1, j - 1);
        }
    }
    return result;
}


/// @brief 将射影空间中的导数转换为欧式空间中的导数
/// @tparam T double float int...
/// @tparam flag 是否是射影空间
/// @tparam rows 坐标点分量的个数
template<typename T, bool flag, int point_size, int n>
struct project_derivs_point
{
    static Eigen::Matrix<Eigen::Vector<T, point_size - 1>, n + 1, n + 1> project_point_to_euclidean_space(Eigen::Matrix<Eigen::Vector<T, point_size>, n + 1, n + 1> const &points)
    {
        Eigen::Matrix<Eigen::Vector<T, point_size - 1>, n + 1, n + 1> SKL;
        Eigen::Matrix<int, n + 1, n + 1> Bin = binary_coeff<n>();
        for (int k = 0; k <= n; ++k)
        {
            for (int l = 0; l <= n - k; ++l)
            {
                Eigen::Vector<T, point_size - 1> v = points(k, l).block(0, 0, point_size - 1, 1);
                for (int j = 1; j <= l; ++j)
                    v -= Bin(l, j) * points(0, j)[point_size - 1] * SKL(k, l - j);
                for (int i = 1; i <= k; ++i)
                {
                    v -= Bin(k, i) * points(i, 0)[point_size - 1] * SKL(k - i, l);
                    Eigen::Vector<T, point_size - 1>  v2;
                    v2.setConstant(0.0);
                    for (int j = 1; j <= l; ++j)
                        v2 = v2 + Bin(l, j) * points(i, j)[point_size - 1] * SKL(k - i, l - j);
                    v -= Bin(k, i) * v2;
                }
                SKL(k, l) = v / points(0, 0)[point_size - 1];
            }
        }
        return SKL;
    }
};

template<typename T, int point_size, int n>
struct project_derivs_point<T, false, point_size, n>
{
    static Eigen::Matrix<Eigen::Vector<T, point_size>, n + 1, n + 1> project_point_to_euclidean_space(Eigen::Matrix<Eigen::Vector<T, point_size>, n + 1, n + 1> const &points)
    {
        return points;
    }
};

template<typename T, int point_size>
struct project_derivs_point<T, true, point_size, -1>
{
    static Eigen::MatrixX<Eigen::Vector<T, point_size - 1>> project_point_to_euclidean_space(Eigen::MatrixX<Eigen::Vector<T, point_size>> const &points)
    {
        int d = points.cols();
        Eigen::MatrixX<Eigen::Vector<T, point_size - 1>> SKL(d, d);
        Eigen::MatrixX<int> Bin = binary_coeff(d);
        for (int k = 0; k < d; ++k)
        {
            for (int l = 0; l < d - k; ++l)
            {
                Eigen::Vector<T, point_size - 1> v = points(k, l).block(0, 0, point_size - 1, 1);
                for (int j = 1; j <= l; ++j)
                    v -= Bin(l, j) * points(0, j)[point_size - 1] * SKL(k, l - j);
                for (int i = 1; i <= k; ++i)
                {
                    v -= Bin(k, i) * points(i, 0)[point_size - 1] * SKL(k - i, l);
                    Eigen::Vector<T, point_size - 1>  v2;
                    v2.setConstant(0.0);
                    for (int j = 1; j <= l; ++j)
                        v2 = v2 + Bin(l, j) * points(i, j)[point_size - 1] * SKL(k - i, l - j);
                    v -= Bin(k, i) * v2;
                }
                SKL(k, l) = v / points(0, 0)[point_size - 1];
            }
        }
        return SKL;
    }
};

template<typename T, int point_size>
struct project_derivs_point<T, false, point_size, -1>
{
    static Eigen::MatrixX<Eigen::Vector<T, point_size>> project_point_to_euclidean_space(Eigen::MatrixX<Eigen::Vector<T, point_size>> const &points)
    {
        return points;
    }
};

//flag == true
template<typename T, bool flag, int point_size>
struct rat_curve_derivs_project
{
    static Eigen::VectorX<Eigen::Vector<T, point_size - 1>> project_point_to_euclidean_space(Eigen::VectorX<Eigen::Vector<T, point_size>> const &points)
    {
        int count = points.size();
        Eigen::MatrixX<int> Bin = binary_coeff(count);
        Eigen::VectorX<Eigen::Vector<T, point_size - 1>>  result(count);
        for (int k = 0; k < count; ++k)
        {
            Eigen::Vector<T, point_size - 1> v = points.block(0, k, point_size - 1, 1);
            for (int i = 1; i <= k; ++i)
                v = v - Bin(k, i) * points(point_size - 1, i) * result[k - i];
            result[k] = v / points(point_size - 1, 0);
        }
    }
};

//flag == false
template<typename T, int point_size>
struct rat_curve_derivs_project<T, false, point_size>
{
    static Eigen::VectorX<Eigen::Vector<T, point_size>> project_point_to_euclidean_space(Eigen::VectorX<Eigen::Vector<T, point_size>> const &points)
    {
        return points;
    }
};


/// @brief 求bezier曲线的基函数B_n^i在参数u处的值
/// @tparam T : double flaot int...
/// @param i : 第i个基函数
/// @param n : 基函数的阶数
/// @param u : 参数u
/// @param B : out_put_param 基函数的值
/// @return ENUM_NURBS错误码
template<typename T>
ENUM_NURBS Bernstein(int i, int n, T u, T &B)
{
    Eigen::Vector<T, Eigen::Dynamic> vec(n);
    vec.setConstant(0);
    vec[n - i] = 1.0;
    T u1 = 1.0 - u;
    for (int k = 1; k < n; ++k)
    {
        vec.block(0, 0, n - k, 1) = u *  vec.block(0, 0, n - k, 1) + u1 *  vec.block(0, 1, n - k, 1);
    }
    B = vec[0];
    return NURBS_SUCCESS;
}

/// @brief 求bezier曲线的在参数u处所有的非零基函数的值
/// @tparam T : float double int
/// @tparam points_count : 控制点的个数
/// @param u : 参数u
/// @param B : out_put_param 基函数的值
/// @return ENUM_NURBS错误码
template<typename T, int points_count>
ENUM_NURBS AllBernstein(T u, Eigen::Vector<T, points_count> &B)
{
    Eigen::Vector<T, points_count + 1> vec;
    vec.setConstant(0);
    vec[1] = 1.0;
    T u1 = 1.0 - u;
    for (int k = 1; k < points_count; ++k)
    {
        //bug : 下列形式的eigen求和不能写出下面注释的形式，只能写出下面三行的形式
        Eigen::Vector<T, Eigen::Dynamic> left = u1 * vec.block(1, 0, k + 1, 1);
        Eigen::Vector<T, Eigen::Dynamic> right = u * vec.block(0, 0, k + 1, 1);
        vec.block(1, 0, k + 1, 1) = left + right;
        // vec.block(1, 0, k + 1, 1) = ((1.0 - u) * vec.block(1, 0, k + 1, 1)) + (u * vec.block(0, 0, k + 1, 1));
    }
    B = vec.block(1, 0, points_count, 1);
    return NURBS_SUCCESS;
}

/// @brief 求bezier曲线的在参数u处所有的非零基函数的值
/// @tparam T : float double int
/// @param u : 参数u
/// @param B : out_put_param 基函数的值
/// @return ENUM_NURBS错误码
template<typename T>
ENUM_NURBS AllBernstein(int n, T u, Eigen::Vector<T, Eigen::Dynamic> &B)
{
    Eigen::Vector<T, Eigen::Dynamic> vec(n + 1);
    vec.setConstant(0);
    vec[1] = 1.0;
    T u1 = 1.0 - u;
    for (int k = 1; k < n; ++k)
    {
        //bug : 下列形式的eigen求和不能写出下面注释的形式，只能写出下面三行的形式
        Eigen::Vector<T, Eigen::Dynamic> left = u1 * vec.block(1, 0, k + 1, 1);
        Eigen::Vector<T, Eigen::Dynamic> right = u * vec.block(0, 0, k + 1, 1);
        vec.block(1, 0, k + 1, 1) = left + right;
        // vec.block(1, 0, k + 1, 1) = ((1.0 - u) * vec.block(1, 0, k + 1, 1)) + (u * vec.block(0, 0, k + 1, 1));
    }
    B.resize(n);
    B = vec.block(1, 0, n, 1);
    return NURBS_SUCCESS;
}

/// @brief bezier曲线的DeCasteljaul算法, 计算bezier曲线在参数u处的向量值
/// @tparam T : float double int
/// @tparam dim : bezier所在欧式空间的维数
/// @tparam cols : 控制点的个数
/// @tparam is_rational : 是否时有理函数
/// @param control_points : 控制点
/// @param u : 参数u
/// @param point : out_put_param bezier在参数u处的向量值
/// @return ENUM_NURBS错误码
template<typename T, int dim, int cols, bool is_rational>
ENUM_NURBS  DeCasteljaul(const Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic> &control_points,
                    T u, Eigen::Vector<T, is_rational ? dim + 1 : dim> &point)
{
    int constexpr rows = is_rational ? dim + 1 : dim;
    Eigen::Matrix<T, rows, cols> local_array = control_points;
    T u1 = 1.0 - u;
    for (int k = 1; k < cols; ++k)
    {
        local_array.block(0, 0, rows, cols - k) = u1 * local_array.block(0, 0, rows, cols - k)
                            + u * local_array.block(0, 1, rows, cols - k);
    }
    point = local_array.block(0, 0, rows, 1);
    return NURBS_SUCCESS;
}

/// @brief bezier曲线的DeCasteljaul算法, 计算bezier曲线在参数u处的向量值
/// @tparam T : float double int
/// @tparam dim : bezier所在欧式空间的维数
/// @tparam is_rational : 是否时有理函数
/// @param control_points : 控制点
/// @param cols : 控制点的个数
/// @param u : 参数u
/// @param point : out_put_param bezier曲线在参数u处的向量值
/// @return ENUM_NURBS错误码
template<typename T, int dim, bool is_rational>
ENUM_NURBS  DeCasteljaul(const Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic> &control_points, int n,
                    T u, Eigen::Vector<T, is_rational ? dim + 1 : dim> &point)
{
    Eigen::Matrix<T, is_rational ? dim + 1 : dim, Eigen::Dynamic> local_array = control_points;
    T u1 = 1.0 - u;
    for (int k = 1; k < n; ++k)
    {
        local_array.block(0, 0, is_rational ? dim + 1 : dim, n - k) = u1 * local_array.block(0, 0, is_rational ? dim + 1 : dim, n - k)
                            + u * local_array.block(0, 1, is_rational ? dim + 1 : dim, n - k);
    }
    point = local_array.block(0, 0, is_rational ? dim + 1 : dim, 1);
    return NURBS_SUCCESS;
}

/// @brief bezier曲面的DeCasteljaul算法, 计算bezier曲线在参数(u, v)处的向量值
/// @tparam T : float double int
/// @tparam dim : bezier所在欧式空间的维数
/// @tparam is_rational : 是否时有理函数
/// @param control_points : 控制点
/// @param rows : 控制点的U向个数(不严谨)
/// @param cols : 控制点的V向个数(不严谨)
/// @param u : 参数u
/// @param v : 参数v
/// @param point : out_put_param bezier曲面在参数(u, v)处的向量值
/// @return ENUM_NURBS错误码
template<typename T, int dim, bool is_rational>
ENUM_NURBS DeCasteljaul2(Eigen::MatrixX<Eigen::Vector<T, is_rational ? dim + 1 : dim>> const &control_points, int rows, int cols,
                    T u, T v, Eigen::Vector<T, is_rational ? dim + 1 : dim> &point)
{
    constexpr int col_size = is_rational ? dim + 1 : dim;
    Eigen::Matrix<T, col_size, Eigen::Dynamic> temp_points(col_size, cols);
    for (int col_index = 0; col_index < cols; ++col_index)
    {
        Eigen::Vector<T, col_size> vec;
        Eigen::Matrix<T, col_size, Eigen::Dynamic> points;
        points.resize(col_size, rows);
        for (int index = 0; index < rows; ++index)
        {
            points.col(index) = control_points(index, col_index);
        }
        DeCasteljaul<T, dim, is_rational>(points, rows, u, vec);
        temp_points.col(col_index) = vec;
    }
    DeCasteljaul<T, dim, is_rational>(temp_points, cols, v, point);
    return NURBS_SUCCESS;
}

/// @brief bezier曲面的DeCasteljaul算法, 计算bezier曲线在参数(u, v)处的向量值
/// @tparam T : float double int
/// @tparam dim : bezier所在欧式空间的维数
/// @tparam is_rational : 是否时有理函数
/// @tparam rows : 控制点的U向个数(不严谨)
/// @tparam cols : 控制点的V向个数(不严谨)
/// @param control_points : 控制点
/// @param u : 参数u
/// @param v : 参数v
/// @param point : out_put_param bezier曲面在参数(u, v)处的向量值
/// @return ENUM_NURBS错误码
template<typename T, int dim, bool is_rational, int rows, int cols>
ENUM_NURBS DeCasteljaul2(Eigen::Matrix<Eigen::Vector<T, is_rational ? dim + 1 : dim>, rows, cols> const &control_points,
                    T u, T v, Eigen::Vector<T, is_rational ? dim + 1 : dim> &point)
{
    constexpr int col_size = is_rational ? dim + 1 : dim;
    Eigen::Matrix<T, col_size, cols> temp_points;
    for (int col_index = 0; col_index < cols; ++col_index)
    {
        Eigen::Vector<T, col_size> vec;
        Eigen::Matrix<T, col_size, rows> points;
        for (int index = 0; index < rows; ++index)
        {
            points.col(index) = control_points(index, col_index);
        }
        DeCasteljaul<T, dim, rows, is_rational>(points, u, vec);
        temp_points.col(col_index) = vec;
    }
    DeCasteljaul<T, dim, cols, is_rational>(temp_points , v, point);
    return NURBS_SUCCESS;
}

/// @brief  计算所有在u处不为0的阶数为0 - degree的基函数的值
/// @tparam T double float int...
/// @tparam degree
/// @tparam points_count 控制点个数
/// @param i u \in [u_i, u_(i + 1))
/// @param u 参数u
/// @param knots_vector
/// @param result result(j, k)(j = 0, ... degree, k = i - j, ..., i) = N_j^k(u); 浪费一半内存
/// @return /// @return ENUM_NURBS错误码
template<typename T, int degree, int points_count>
ENUM_NURBS all_basis_functions(int i, T u, const Eigen::Vector<T, degree + 1 + points_count> &knots_vector, Eigen::Matrix<T, degree + 1, degree + 1> &result)
{
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> left;
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> right;
    // right[degree] = 0.0;
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, degree + 2> iterArray;
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    result.setConstant(0.0);
    result(0, degree) = 1.0;
    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;
        result.block(0, degree - iter_step, iter_step + 1, 1) = iterArray.block(degree - iter_step, 0, iter_step + 1, 1);
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief  计算所有在u处不为0的阶数为0 - degree的基函数的值
/// @tparam T double float int...
/// @tparam degree
/// @param i u \in [u_i, u_(i + 1))
/// @param u 参数u
/// @param knots_vector
/// @param result result(j, k)(j = 0, ... degree, k = i - j, ..., i) = N_j^k(u); 浪费一半内存
/// @return /// @return ENUM_NURBS错误码
template<typename T, int degree>
ENUM_NURBS all_basis_functions(int i, T u, const Eigen::VectorX<T> &knots_vector, Eigen::Matrix<T, degree + 1, degree + 1> &result)
{
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> left;
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> right;
    // right[degree] = 0.0;
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, degree + 2> iterArray;
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    result.setConstant(0.0);
    result(0, degree) = 1.0;
    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;
        result.block(0, degree - iter_step, iter_step + 1, 1) = iterArray.block(degree - iter_step, 0, iter_step + 1, 1);
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief  计算所有在u处不为0的阶数为0 - degree的基函数的值
/// @tparam T double float int...
/// @param i u \in [u_i, u_(i + 1))
/// @param degree
/// @param u 参数u
/// @param knots_vector
/// @param result result(j, k)(j = 0, ... degree, k = i - j, ..., i) = N_j^k(u); 浪费一半内存
/// @return /// @return ENUM_NURBS错误码
template<typename T>
ENUM_NURBS all_basis_functions(int i, int degree, T u, const Eigen::VectorX<T> &knots_vector, Eigen::MatrixX<T> &result)
{
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    result.resize(degree + 1, degree + 1);
    Eigen::VectorX<T> left(degree);
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::VectorX<T> right(degree);
    // right[degree] = 0.0;
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::VectorX<T> iterArray(degree + 2);
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    result.setConstant(0.0);
    result(0, degree) = 1.0;
    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;
        result.block(0, degree - iter_step, iter_step + 1, 1) = iterArray.block(degree - iter_step, 0, iter_step + 1, 1);
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief  计算所有在u处不为0的阶数为degree的基函数的值
/// @tparam T double float int...
/// @tparam degree
/// @tparam points_count 控制点个数
/// @param i u \in [u_i, u_(i + 1))
/// @param u 参数u
/// @param knots_vector
/// @param result result(j)(j = 0, ... degree) = N_degree^k (k = i - degree, ... i)
/// @return /// @return ENUM_NURBS错误码
template<typename T, int degree, int points_count>
ENUM_NURBS basis_functions(int i, T u, const Eigen::Vector<T, degree + 1 + points_count> &knots_vector, Eigen::Vector<T, degree + 1> &result)
{
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> left;
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> right;
    // right[degree] = 0.0;
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, degree + 2> iterArray;
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;
    }
    result = iterArray.block(0, 0, degree + 1, 1);
    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief  计算所有在u处不为0的阶数为degree的基函数的值
/// @tparam T double float int...
/// @tparam degree
/// @param i u \in [u_i, u_(i + 1))
/// @param u 参数u
/// @param knots_vector
/// @param result result(j)(j = 0, ... degree) = N_degree^k (k = i - degree, ... i)
/// @return ENUM_NURBS错误码
template<typename T, int degree>
ENUM_NURBS basis_functions(int i, T u, const Eigen::Vector<T, Eigen::Dynamic> &knots_vector, Eigen::Vector<T, degree + 1> &result)
{
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> left;
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> right;
    // right[degree] = 0.0;
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, degree + 2> iterArray;
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;
    }
    result = iterArray.block(0, 0, degree + 1, 1);
    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief  计算所有在u处不为0的阶数为degree的基函数的值
/// @tparam T double float int...
/// @param i u \in [u_i, u_(i + 1))
/// @param u 参数u
/// @param degree
/// @param knots_vector
/// @param result result(j)(j = 0, ... degree) = N_degree^k (k = i - degree, ... i)
/// @return /// @return ENUM_NURBS错误码
template<typename T>
ENUM_NURBS basis_functions(int i, T u, int degree, const Eigen::Vector<T, Eigen::Dynamic> &knots_vector, Eigen::Vector<T, Eigen::Dynamic> &result)
{
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, Eigen::Dynamic> left(degree);
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, Eigen::Dynamic> right(degree);
    // right[degree] = 0.0;
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, Eigen::Dynamic> iterArray(degree + 2);
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;
    }
    result = iterArray.block(0, 0, degree + 1, 1);
    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief 查找u \in [u_index, u_(index + 1))的i; remark: 当u == 节点矢量的右端点时, index = knots_vector.size() - degree - 2;
///         可以看出此函数优先计算右极限, 对于右端点其计算的为左极限
/// @tparam T double float int...
/// @tparam points_count 控制点个数
/// @tparam degree nurbs的阶数
/// @param u 参数
/// @param knots_vector nurbs的节点矢量
/// @param index out_put_param u \in [u_index, u_(index + 1))
/// @return ENUM_NURBS错误码
template<typename T, int points_count, int degree>
ENUM_NURBS find_span(T u, Eigen::Vector<T, degree + 1 + points_count> const  &knots_vector, int &index)
{
    if (u > knots_vector[points_count] || u < knots_vector[0])
        return ENUM_NURBS::NURBS_PARAM_IS_OUT_OF_DOMAIN;
    if (u == knots_vector[points_count])
    {
        index = points_count - 1;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    int mid = (degree + points_count) / 2;
    int low = degree;
    int high = points_count;
    while (u < knots_vector[mid] || u >= knots_vector[mid + 1])
    {
        if (u < knots_vector[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }
    index = mid;
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief 查找u \in [u_index, u_(index + 1))的i; remark: 当u == 节点矢量的右端点时, index = knots_vector.size() - degree - 2;
///         可以看出此函数优先计算右极限, 对于右端点其计算的为左极限
/// @tparam T double float int...
/// @tparam degree nurbs的阶数
/// @param u 参数
/// @param knots_vector nurbs的节点矢量
/// @param index out_put_param u \in [u_index, u_(index + 1))
/// @return ENUM_NURBS错误码
template<typename T, int degree>
ENUM_NURBS find_span(T u, Eigen::Vector<T, Eigen::Dynamic> const &knots_vector, int &index)
{
    int points_count = knots_vector.size() - degree - 1;
    if (u > knots_vector[points_count] || u < knots_vector[0])
        return ENUM_NURBS::NURBS_PARAM_IS_OUT_OF_DOMAIN;
    if (u == knots_vector[points_count])
    {
        index = points_count - 1;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    int mid = (degree + points_count) / 2;
    int low = degree;
    int high = points_count;
    while (u < knots_vector[mid] || u >= knots_vector[mid + 1])
    {
        if (u < knots_vector[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }
    index = mid;
    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief 查找u \in [u_index, u_(index + 1))的i; remark: 当u == 节点矢量的右端点时, index = knots_vector.size() - degree - 2;
///         可以看出此函数优先计算右极限, 对于右端点其计算的为左极限
/// @tparam T double float int...
/// @param u 参数
/// @param degree nurbs的阶数
/// @param knots_vector nurbs的节点矢量
/// @param index out_put_param u \in [u_index, u_(index + 1))
/// @return ENUM_NURBS错误码
template<typename T>
ENUM_NURBS find_span(T u, int degree, Eigen::Vector<T, Eigen::Dynamic> const &knots_vector, int &index)
{
    int points_count = knots_vector.size() - degree - 1;
    if (u > knots_vector[points_count] || u < knots_vector[0])
        return ENUM_NURBS::NURBS_PARAM_IS_OUT_OF_DOMAIN;
    if (u == knots_vector[points_count])
    {
        index = points_count - 1;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    int mid = (degree + points_count) / 2;
    int low = degree;
    int high = points_count;
    while (u < knots_vector[mid] || u >= knots_vector[mid + 1])
    {
        if (u < knots_vector[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }
    index = mid;
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief 元模板计算n!
/// @tparam n
template<int n>
struct my_factorial
{
    static constexpr int value = n * my_factorial<n - 1>::value;
};
template<>
struct my_factorial<1>
{
    static constexpr int value = 1;
};

/// @brief 求nurbs的基函数的导数
/// @tparam T : float double ...
/// @tparam degree : nurbs基函数的次数
/// @tparam points_count : nurbs控制点的个数
/// @tparam n : 求(0, 1, 2 ... n)阶导数
/// @param i : knots的index
/// @param u :曲线的参数u
/// @param knots_vector :nurbs的节点矢量
/// @param result out_put_param : result(k, l)为基函数N_(i - p + k)的第l次导数
/// @return ENUM_NURBS错误码
template<typename T, int degree, int points_count, int n>
ENUM_NURBS ders_basis_funs(int i, T u, const Eigen::Vector<T, degree + 1 + points_count> &knots_vector,
                Eigen::Matrix<T, degree + 1, n + 1> &result)
{
    result.setConstant(0.0);
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> left;
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> right;

    Eigen::Matrix<T, degree + 1, degree + 1> ndu;
    ndu.setConstant(0.0);
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, degree + 2> iterArray;
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    ndu(0, 0) = 1.0;

    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;

        ndu.block(0, iter_step, iter_step + 1, 1) = left_part + right_part;
        ndu.block(iter_step, 0, 1, iter_step) = coeff_denominator.transpose();
    }

    result.col(0) = ndu.col(degree);

    for (int r = i - degree; r <= i; ++r) //对基函数进行循环
    {
        // Eigen::Vector<T, degree + 1> current_array;
        Eigen::Vector<Eigen::Vector<T, degree + 1>, 2> arrays;
        int current_index = 0;
        int next_index = 1;
        Eigen::Vector<T, degree + 1> &current_array = arrays[0];

        arrays[0].setConstant(0.0);
        arrays[0][0] = 1.0;
        for (int k = 1; k <= n; ++k)
        {
            Eigen::Vector<T, degree + 1> &current_array = arrays[current_array];
            Eigen::Vector<T, degree + 1> &next_array = arrays[next_index];
            next_array.setConstant(0.0);

            int left_num = r - (i - degree) - k;
            int left_index = std::max(0, -left_num);
            int right_num = r - (i - degree);
            int right_index = right_num > (degree - k) ? degree - k - left_num : right_num - left_num;
            int next_array_length = right_index - left_index + 1;

            // int col_of_array = std::max(0, degree - i + r - 1);
            int col_of_array = std::max(0, degree - i + r - k);
            Eigen::Array<T, Eigen::Dynamic, 1> denominator = ndu.block(degree + 1 - k, col_of_array, 1, next_array_length).transpose();
            int current_left_index = left_index;
            int current_right_index = right_index;
            if (left_index == 0)
            {
                next_array[0] = current_array[0] / denominator(0, 0);
                current_left_index += 1;
            }
            if (right_index == k)
            {
                next_array[k] = -current_array[k - 1] / denominator(next_array_length - 1, 0);
                current_right_index -= 1;
            }
            if (current_right_index >= current_left_index)
            {
                int arrayLength = current_right_index - current_left_index + 1;
                next_array.block(current_left_index, 0, arrayLength, 1) = (current_array.block(current_left_index, 0, arrayLength, 1) -
                    current_array.block(current_left_index - 1, 0, arrayLength, 1)).array() / denominator.block(current_left_index - left_index, 0, arrayLength, 1).array();
            }

            left_num = std::max(0, left_num);
            Eigen::VectorX<T> temp_r = ndu.block(left_num, degree - k,  next_array_length, 1);
            Eigen::VectorX<T> temp_l = next_array.block(left_index, 0, next_array_length, 1);
            result(r - (i - degree), k) = my_factorial<degree>::value / std::tgamma<int>(degree - k + 1) *
                        temp_l.dot(temp_r);

            std::swap(current_index, next_index);
        }
    }

    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief 求nurbs的基函数的导数
/// @tparam T : float double ...
/// @tparam degree : nurbs基函数的次数
/// @tparam points_count : nurbs控制点的个数
/// @param i : knots的index
/// @param n : 求(0, 1, 2 ... n)阶导数
/// @param u :曲线的参数u
/// @param knots_vector :nurbs的节点矢量
/// @param result out_put_param : result(k, l)为基函数N_(i - p + k)的第l次导数
/// @return ENUM_NURBS错误码
template<typename T, int degree, int points_count>
ENUM_NURBS ders_basis_funs(int i, int n, T u, const Eigen::Vector<T, degree + 1 + points_count> &knots_vector,
                Eigen::Matrix<T, degree + 1, Eigen::Dynamic, Eigen::ColMajor, degree + 1, degree + 1> &result)
{
    result.resize(degree + 1, n + 1);
    result.setConstant(0.0);
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> left;
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> right;

    Eigen::Matrix<T, degree + 1, degree + 1> ndu;
    ndu.setConstant(0.0);
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, degree + 2> iterArray;
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    ndu(0, 0) = 1.0;

    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;

        ndu.block(0, iter_step, iter_step + 1, 1) = left_part + right_part;
        ndu.block(iter_step, 0, 1, iter_step) = coeff_denominator.transpose();
    }

    result.col(0) = ndu.col(degree);

    for (int r = i - degree; r <= i; ++r) //对基函数进行循环
    {
        Eigen::Vector<Eigen::Vector<T, degree + 1>, 2> arrays;
        int current_index = 0;
        int next_index = 1;

        arrays[0].setConstant(0.0);
        arrays[0][0] = 1.0;
        for (int k = 1; k <= n; ++k)
        {
            Eigen::Vector<T, degree + 1> &current_array = arrays[current_index];
            Eigen::Vector<T, degree + 1> &next_array = arrays[next_index];
            next_array.setConstant(0.0);

            int left_num = r - (i - degree) - k;
            int left_index = std::max(0, -left_num);
            int right_num = r - (i - degree);
            int right_index = right_num > (degree - k) ? degree - k - left_num : right_num - left_num;
            int next_array_length = right_index - left_index + 1;

            // int col_of_array = std::max(0, degree - i + r - 1);
            int col_of_array = std::max(0, degree - i + r - k);
            Eigen::Array<T, Eigen::Dynamic, 1> denominator = ndu.block(degree + 1 - k, col_of_array, 1, next_array_length).transpose();
            int current_left_index = left_index;
            int current_right_index = right_index;
            if (left_index == 0)
            {
                next_array[0] = current_array[0] / denominator(0, 0);
                current_left_index += 1;
            }
            if (right_index == k)
            {
                next_array[k] = -current_array[k - 1] / denominator(next_array_length - 1, 0);
                current_right_index -= 1;
            }
            if (current_right_index >= current_left_index)
            {
                int arrayLength = current_right_index - current_left_index + 1;
                next_array.block(current_left_index, 0, arrayLength, 1) = (current_array.block(current_left_index, 0, arrayLength, 1) -
                    current_array.block(current_left_index - 1, 0, arrayLength, 1)).array() / denominator.block(current_left_index - left_index, 0, arrayLength, 1).array();
            }

            left_num = std::max(0, left_num);
            Eigen::VectorX<T> temp_r = ndu.block(left_num, degree - k,  next_array_length, 1);
            Eigen::VectorX<T> temp_l = next_array.block(left_index, 0, next_array_length, 1);
            result(r - (i - degree), k) = my_factorial<degree>::value / std::tgamma<int>(degree - k + 1) *
                        temp_l.dot(temp_r);

            std::swap(current_index, next_index);
        }
    }

    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief 求nurbs的基函数的导数
/// @tparam T : float double ...
/// @tparam degree : nurbs基函数的次数
/// @param i : knots的index
/// @param n : 求(0, 1, 2 ... n)阶导数
/// @param u :曲线的参数u
/// @param knots_vector :nurbs的节点矢量
/// @param result out_put_param : result(k, l)为基函数N_(i - p + k)的第l次导数
/// @return ENUM_NURBS错误码
template<typename T, int degree>
ENUM_NURBS ders_basis_funs(int i, int n, T u, const Eigen::Vector<T, Eigen::Dynamic> &knots_vector,
                Eigen::Matrix<T, degree + 1, Eigen::Dynamic, Eigen::ColMajor, degree + 1, degree + 1> &result)
{
    result.resize(degree + 1, n + 1);
    result.setConstant(0.0);
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> left;
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> right;

    Eigen::Matrix<T, degree + 1, degree + 1> ndu;
    ndu.setConstant(0.0);
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, degree + 2> iterArray;
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    ndu(0, 0) = 1.0;

    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;

        ndu.block(0, iter_step, iter_step + 1, 1) = left_part + right_part;
        ndu.block(iter_step, 0, 1, iter_step) = coeff_denominator.transpose();
    }

    result.col(0) = ndu.col(degree);

    for (int r = i - degree; r <= i; ++r) //对基函数进行循环
    {
        // Eigen::Vector<T, degree + 1> current_array;
        Eigen::Vector<Eigen::Vector<T, degree + 1>, 2> arrays;
        int current_index = 0;
        int next_index = 1;

        arrays[0].setConstant(0.0);
        arrays[0][0] = 1.0;
        for (int k = 1; k <= n; ++k)
        {
            Eigen::Vector<T, degree + 1> &current_array = arrays[current_index];
            Eigen::Vector<T, degree + 1> &next_array = arrays[next_index];
            next_array.setConstant(0.0);

            int left_num = r - (i - degree) - k;
            int left_index = std::max(0, -left_num);
            int right_num = r - (i - degree);
            int right_index = right_num > (degree - k) ? degree - k - left_num : right_num - left_num;
            int next_array_length = right_index - left_index + 1;

            // int col_of_array = std::max(0, degree - i + r - 1);
            int col_of_array = std::max(0, degree - i + r - k);
            Eigen::Array<T, Eigen::Dynamic, 1> denominator = ndu.block(degree + 1 - k, col_of_array, 1, next_array_length).transpose();
            int current_left_index = left_index;
            int current_right_index = right_index;
            if (left_index == 0)
            {
                next_array[0] = current_array[0] / denominator(0, 0);
                current_left_index += 1;
            }
            if (right_index == k)
            {
                next_array[k] = -current_array[k - 1] / denominator(next_array_length - 1, 0);
                current_right_index -= 1;
            }
            if (current_right_index >= current_left_index)
            {
                int arrayLength = current_right_index - current_left_index + 1;
                next_array.block(current_left_index, 0, arrayLength, 1) = (current_array.block(current_left_index, 0, arrayLength, 1) -
                    current_array.block(current_left_index - 1, 0, arrayLength, 1)).array() / denominator.block(current_left_index - left_index, 0, arrayLength, 1).array();
            }

            left_num = std::max(0, left_num);
            Eigen::VectorX<T> temp_r = ndu.block(left_num, degree - k,  next_array_length, 1);
            Eigen::VectorX<T> temp_l = next_array.block(left_index, 0, next_array_length, 1);
            result(r - (i - degree), k) = my_factorial<degree>::value / std::tgamma<int>(degree - k + 1) *
                        temp_l.dot(temp_r);

            std::swap(current_index, next_index);
        }
    }

    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief 求nurbs的基函数的导数
/// @tparam T : float double ...
/// @tparam degree : nurbs基函数的次数
/// @tparam n : 求(0, 1, 2 ... n)阶导数
/// @param i : knots的index
/// @param u :曲线的参数u
/// @param knots_vector :nurbs的节点矢量
/// @param result out_put_param : result(k, l)为基函数N_(i - p + k)的第l次导数
/// @return ENUM_NURBS错误码
template<typename T, int degree, int n>
ENUM_NURBS ders_basis_funs(int i, T u, const Eigen::Vector<T, Eigen::Dynamic> &knots_vector,
                Eigen::Matrix<T, degree + 1, n + 1> &result)
{
    result.setConstant(0.0);
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> left;
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, degree> right;

    Eigen::Matrix<T, degree + 1, degree + 1> ndu;
    ndu.setConstant(0.0);
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, degree + 2> iterArray;
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    ndu(0, 0) = 1.0;

    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;

        ndu.block(0, iter_step, iter_step + 1, 1) = left_part + right_part;
        ndu.block(iter_step, 0, 1, iter_step) = coeff_denominator.transpose();
    }

    result.col(0) = ndu.col(degree);

    for (int r = i - degree; r <= i; ++r) //对基函数进行循环
    {
        // Eigen::Vector<T, degree + 1> current_array;
        Eigen::Vector<Eigen::Vector<T, degree + 1>, 2> arrays;
        int current_index = 0;
        int next_index = 1;

        arrays[0].setConstant(0.0);
        arrays[0][0] = 1.0;
        for (int k = 1; k <= n; ++k)
        {
            Eigen::Vector<T, degree + 1> &current_array = arrays[current_index];
            Eigen::Vector<T, degree + 1> &next_array = arrays[next_index];
            next_array.setConstant(0.0);

            int left_num = r - (i - degree) - k;
            int left_index = std::max(0, -left_num);
            int right_num = r - (i - degree);
            int right_index = right_num > (degree - k) ? degree - k - left_num : right_num - left_num;
            int next_array_length = right_index - left_index + 1;

            // int col_of_array = std::max(0, degree - i + r - 1);
            int col_of_array = std::max(0, degree - i + r - k);
            Eigen::Array<T, Eigen::Dynamic, 1> denominator = ndu.block(degree + 1 - k, col_of_array, 1, next_array_length).transpose();
            int current_left_index = left_index;
            int current_right_index = right_index;
            if (left_index == 0)
            {
                next_array[0] = current_array[0] / denominator(0, 0);
                current_left_index += 1;
            }
            if (right_index == k)
            {
                next_array[k] = -current_array[k - 1] / denominator(next_array_length - 1, 0);
                current_right_index -= 1;
            }
            if (current_right_index >= current_left_index)
            {
                int arrayLength = current_right_index - current_left_index + 1;
                next_array.block(current_left_index, 0, arrayLength, 1) = (current_array.block(current_left_index, 0, arrayLength, 1) -
                    current_array.block(current_left_index - 1, 0, arrayLength, 1)).array() / denominator.block(current_left_index - left_index, 0, arrayLength, 1).array();
            }

            left_num = std::max(0, left_num);
            Eigen::VectorX<T> temp_r = ndu.block(left_num, degree - k,  next_array_length, 1);
            Eigen::VectorX<T> temp_l = next_array.block(left_index, 0, next_array_length, 1);
            result(r - (i - degree), k) = my_factorial<degree>::value / std::tgamma<int>(degree - k + 1) *
                        temp_l.dot(temp_r);

            std::swap(current_index, next_index);
        }
    }

    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief 求nurbs的基函数的导数
/// @tparam T : float double ...
/// @param i : knots的index
/// @param n : 求(0, 1, 2 ... n)阶导数
/// @param degree : nurbs基函数的次数
/// @param u :曲线的参数u
/// @param knots_vector :nurbs的节点矢量
/// @param result out_put_param : result(k, l)为基函数N_(i - p + k)的第l次导数
/// @return ENUM_NURBS错误码
template<typename T>
ENUM_NURBS ders_basis_funs(int i, int n, int degree, T u, const Eigen::Vector<T, Eigen::Dynamic> &knots_vector,
                Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &result)
{
    result.resize(degree + 1, n + 1);
    result.setConstant(0.0);
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, Eigen::Dynamic> left(degree);
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, Eigen::Dynamic> right(degree);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ndu(degree + 1, degree + 1);
    ndu.setConstant(0.0);
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, Eigen::Dynamic> iterArray(degree + 2);
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    ndu(0, 0) = 1.0;

    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;

        ndu.block(0, iter_step, iter_step + 1, 1) = left_part + right_part;
        ndu.block(iter_step, 0, 1, iter_step) = coeff_denominator.transpose();
    }

    result.col(0) = ndu.col(degree);

    for (int r = i - degree; r <= i; ++r) //对基函数进行循环
    {
        Eigen::Vector<Eigen::Vector<T, Eigen::Dynamic>, 2> arrays;
        arrays[0].resize(degree + 1);
        arrays[0].setConstant(0.0);
        arrays[1].resize(degree + 1);
        arrays[1].setConstant(0.0);
        int current_index = 0;
        int next_index = 1;
        
        arrays[0].setConstant(0.0);
        arrays[0][0] = 1.0;
        for (int k = 1; k <= n; ++k)
        {
            Eigen::Vector<T, Eigen::Dynamic> &current_array = arrays[current_index];
            Eigen::Vector<T, Eigen::Dynamic> &next_array = arrays[next_index];
            next_array.setConstant(0.0);

            int left_num = r - (i - degree) - k;
            int left_index = std::max(0, -left_num);
            int right_num = r - (i - degree);
            int right_index = right_num > (degree - k) ? degree - k - left_num : right_num - left_num;
            int next_array_length = right_index - left_index + 1;

            // int col_of_array = std::max(0, degree - i + r - 1);
            int col_of_array = std::max(0, degree - i + r - k);
            Eigen::Array<T, Eigen::Dynamic, 1> denominator = ndu.block(degree + 1 - k, col_of_array, 1, next_array_length).transpose();
            int current_left_index = left_index;
            int current_right_index = right_index;
            if (left_index == 0)
            {
                next_array[0] = current_array[0] / denominator(0, 0);
                current_left_index += 1;
            }
            if (right_index == k)
            {
                next_array[k] = -current_array[k - 1] / denominator(next_array_length - 1, 0);
                current_right_index -= 1;
            }
            if (current_right_index >= current_left_index)
            {
                int arrayLength = current_right_index - current_left_index + 1;
                next_array.block(current_left_index, 0, arrayLength, 1) = (current_array.block(current_left_index, 0, arrayLength, 1) -
                    current_array.block(current_left_index - 1, 0, arrayLength, 1)).array() / denominator.block(current_left_index - left_index, 0, arrayLength, 1).array();
            }

            left_num = std::max(0, left_num);
            Eigen::VectorX<T> temp_r = ndu.block(left_num, degree - k,  next_array_length, 1);
            Eigen::VectorX<T> temp_l = next_array.block(left_index, 0, next_array_length, 1);
            result(r - (i - degree), k) = std::tgamma<int>(degree + 1) / std::tgamma<int>(degree - k + 1) *
                        temp_l.dot(temp_r);

            std::swap(current_index, next_index);
        }
    }

    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief 求nurbs的基函数的导数
/// @tparam T : float double ...
/// @tparam n : 求(0, 1, 2 ... n)阶导数
/// @param i : knots的index
/// @param degree : nurbs基函数的次数
/// @param u :曲线的参数u
/// @param knots_vector :nurbs的节点矢量
/// @param result out_put_param : result(k, l)为基函数N_(i - p + k)的第l次导数
/// @return ENUM_NURBS错误码
template<typename T, int n>
ENUM_NURBS ders_basis_funs(int i, int degree, T u, const Eigen::Vector<T, Eigen::Dynamic> &knots_vector,
                Eigen::Matrix<T, Eigen::Dynamic, n + 1> &result)
{
    result.resize(degree + 1, n + 1);
    result.setConstant(0.0);
    // 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
    //将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, Eigen::Dynamic> left(degree);
    // left[0] = 0.0;
    //将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
    Eigen::Vector<T, Eigen::Dynamic> right(degree);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ndu(degree + 1, degree + 1);
    ndu.setConstant(0.0);
    for (int index = 0; index < degree; ++index)
    {
        left[index] = u - knots_vector[i - degree + 1 + index];
        right[index] = knots_vector[i + 1 + index] - u;
    }
    Eigen::Vector<T, Eigen::Dynamic> iterArray(degree + 2);
    iterArray.setConstant(0.0);
    iterArray[degree] = 1.0;
    ndu(0, 0) = 1.0;

    for (int iter_step = 1; iter_step <= degree; ++iter_step)
    {
        //u - u_j 列
        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
        first_coeff[0] = 0.0;
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
        second_coeff[iter_step] = 0.0;

        Eigen::Array<T, Eigen::Dynamic, 1> first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

        Eigen::Array<T, Eigen::Dynamic, 1> coeff_denominator = first_coeff_numerator + second_coeff_numerator;

        first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
        second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
        Eigen::Array<T, Eigen::Dynamic, 1> left_part = first_coeff * iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
        Eigen::Array<T, Eigen::Dynamic, 1> right_part = second_coeff * iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
        iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = left_part + right_part;

        ndu.block(0, iter_step, iter_step + 1, 1) = left_part + right_part;
        ndu.block(iter_step, 0, 1, iter_step) = coeff_denominator.transpose();
    }

    result.col(0) = ndu.col(degree);

    for (int r = i - degree; r <= i; ++r) //对基函数进行循环
    {
        Eigen::Vector<Eigen::Vector<T, Eigen::Dynamic>, 2> arrays;
        arrays[0].resize(degree + 1);
        arrays[1].resize(degree + 1);
        int current_index = 0;
        int next_index = 1;

        arrays[0].setConstant(0.0);
        arrays[0][0] = 1.0;
        for (int k = 1; k <= n; ++k)
        {
            Eigen::Vector<T, Eigen::Dynamic> &current_array = arrays[current_index];
            Eigen::Vector<T, Eigen::Dynamic> &next_array = arrays[next_index];
            next_array.setConstant(0.0);

            int left_num = r - (i - degree) - k;
            int left_index = std::max(0, -left_num);
            int right_num = r - (i - degree);
            int right_index = right_num > (degree - k) ? degree - k - left_num : right_num - left_num;
            int next_array_length = right_index - left_index + 1;

            // int col_of_array = std::max(0, degree - i + r - 1);
            int col_of_array = std::max(0, degree - i + r - k);
            Eigen::Array<T, Eigen::Dynamic, 1> denominator = ndu.block(degree + 1 - k, col_of_array, 1, next_array_length).transpose();
            int current_left_index = left_index;
            int current_right_index = right_index;
            if (left_index == 0)
            {
                next_array[0] = current_array[0] / denominator(0, 0);
                current_left_index += 1;
            }
            if (right_index == k)
            {
                next_array[k] = -current_array[k - 1] / denominator(next_array_length - 1, 0);
                current_right_index -= 1;
            }
            if (current_right_index >= current_left_index)
            {
                int arrayLength = current_right_index - current_left_index + 1;
                next_array.block(current_left_index, 0, arrayLength, 1) = (current_array.block(current_left_index, 0, arrayLength, 1) -
                    current_array.block(current_left_index - 1, 0, arrayLength, 1)).array() / denominator.block(current_left_index - left_index, 0, arrayLength, 1).array();
            }

            left_num = std::max(0, left_num);
            Eigen::VectorX<T> temp_r = ndu.block(left_num, degree - k,  next_array_length, 1);
            Eigen::VectorX<T> temp_l = next_array.block(left_index, 0, next_array_length, 1);
            result(r - (i - degree), k) = std::tgamma<int>(degree + 1) / std::tgamma<int>(degree - k + 1) *
                        temp_l.dot(temp_r);

            std::swap(current_index, next_index);
        }
    }

    return ENUM_NURBS::NURBS_SUCCESS;
}


template<int r1, int r2>
struct r1r2
{
    static_assert(r1 <= r2, "r1 is greater than r2");
    static constexpr int len = r2 - r1 + 1;
};

/// @brief 计算nurbs曲线求导后的新控制点
/// @tparam T double float int...
/// @tparam control_points_count : 控制点个数
/// @tparam degree : nurbs曲线的阶数
/// @tparam n : 求到n次导
/// @tparam rows : 坐标点分量的个数
/// @tparam r1 : 开始控制点的下标
/// @tparam r2 : 结束控制点的下标
/// @param knots_vector : 节点矢量
/// @param control_points : 控制点
/// @param PK : out_put_param PK[i](*,k)表示i阶导曲线的第k + r1个控制点
/// @return ENUM_NURBS错误码
template<typename T, int control_points_count, int degree, int n, int rows, int r1, int r2>
ENUM_NURBS curve_deriv_cpts(Eigen::Vector<T, control_points_count + degree + 1> const &knots_vector,
                    Eigen::Matrix<T, rows, control_points_count> const &control_points,
                    Eigen::Vector<Eigen::Matrix<T, rows, r1r2<r1, r2>::len>, n + 1> &PK)
{

    constexpr int r = r2 - r1;
    static_assert(r <= degree, "r2 - r1 > degree");
    static_assert(n <= degree, "n > degree");
    PK[0] = control_points.template block<rows, r1r2<r1, r2>::len>(0, r1);
    Eigen::Vector3d x{1, 32, 3};
    auto y = x.block<r1r2<r1, r2>::len - 2, 1>(0, 0);
    for (int k = 1; k <= n; ++k) //求导阶数
    {
        int temp = degree - k + 1;
        Eigen::Array<T, Eigen::Dynamic, 1> denominator = knots_vector.block(r1 + degree + 1, 0, r - k + 1, 1) - knots_vector.block(r1 + k, 0, r - k + 1, 1);
        PK[k].resize(rows, r + 1);
        PK[k].block(0, 0, rows, r - k + 1) = (temp * (PK[k - 1].block(0, 1, rows, r - k + 1) - PK[k - 1].block(0, 0, rows, r - k + 1))).array().rowwise() / denominator.transpose();
    }
    return  NURBS_SUCCESS;
}

/// @brief 计算nurbs曲线求导后的新控制点
/// @tparam T double float int...
/// @tparam degree : nurbs曲线的阶数
/// @tparam n : 求到n次导
/// @tparam rows : 坐标点分量的个数
/// @tparam r1 : 开始控制点的下标
/// @tparam r2 : 结束控制点的下标
/// @param knots_vector : 节点矢量
/// @param control_points : 控制点
/// @param PK : out_put_param PK[i](*,k)表示i阶导曲线的第k + r1个控制点
/// @return ENUM_NURBS错误码
template<typename T, int degree, int n, int rows, int r1, int r2>
ENUM_NURBS curve_deriv_cpts(Eigen::VectorX<T> const &knots_vector,
                    Eigen::Matrix<T, rows, Eigen::Dynamic> const &control_points,
                    Eigen::Vector<Eigen::Matrix<T, rows, r1r2<r1, r2>::len>, n + 1> &PK)
{

    constexpr int r = r2 - r1;
    static_assert(r <= degree, "r2 - r1 > degree");
    static_assert(n <= degree, "n > degree");
    PK[0] = control_points.template block<rows, r1r2<r1, r2>::len>(0, r1);
    Eigen::Vector3d x{1, 32, 3};
    auto y = x.block<r1r2<r1, r2>::len - 2, 1>(0, 0);
    for (int k = 1; k <= n; ++k) //求导阶数
    {
        int temp = degree - k + 1;
        Eigen::Array<T, Eigen::Dynamic, 1> denominator = knots_vector.block(r1 + degree + 1, 0, r - k + 1, 1) - knots_vector.block(r1 + k, 0, r - k + 1, 1);
        PK[k].resize(rows, r + 1);
        PK[k].block(0, 0, rows, r - k + 1) = (temp * (PK[k - 1].block(0, 1, rows, r - k + 1) - PK[k - 1].block(0, 0, rows, r - k + 1))).array().rowwise() / denominator.transpose();
    }
    return  NURBS_SUCCESS;
}


/// @brief 计算nurbs曲线求导后的新控制点
/// @tparam T double float int...
/// @tparam degree : nurbs曲线的阶数
/// @tparam rows : 坐标点分量的个数
/// @tparam r1 : 开始控制点的下标
/// @tparam r2 : 结束控制点的下标
/// @param n : 求到n次导
/// @param knots_vector : 节点矢量
/// @param control_points : 控制点
/// @param PK : out_put_param PK[i](*,k)表示i阶导曲线的第k + r1个控制点
/// @return ENUM_NURBS错误码
template<typename T, int degree, int rows, int r1, int r2>
ENUM_NURBS curve_deriv_cpts(int n, Eigen::VectorX<T> const &knots_vector,
                    Eigen::Matrix<T, rows, Eigen::Dynamic> const &control_points,
                    Eigen::VectorX<Eigen::Matrix<T, rows, r1r2<r1, r2>::len>> &PK)
{

    constexpr int r = r2 - r1;
    static_assert(r <= degree, "r2 - r1 > degree");
    PK[0] = control_points.template block<rows, r1r2<r1, r2>::len>(0, r1);
    Eigen::Vector3d x{1, 32, 3};
    auto y = x.block<r1r2<r1, r2>::len - 2, 1>(0, 0);
    for (int k = 1; k <= n; ++k) //求导阶数
    {
        int temp = degree - k + 1;
        Eigen::Array<T, Eigen::Dynamic, 1> denominator = knots_vector.block(r1 + degree + 1, 0, r - k + 1, 1) - knots_vector.block(r1 + k, 0, r - k + 1, 1);
        PK[k].resize(rows, r + 1);
        PK[k].block(0, 0, rows, r - k + 1) = (temp * (PK[k - 1].block(0, 1, rows, r - k + 1) - PK[k - 1].block(0, 0, rows, r - k + 1))).array().rowwise() / denominator.transpose();
    }
    return  NURBS_SUCCESS;
}

/// @brief 计算nurbs曲线求导后的新控制点
/// @tparam T double float int...
/// @tparam control_points_count : 控制点个数
/// @tparam degree : nurbs曲线的阶数
/// @tparam n : 求到n次导
/// @tparam rows : 坐标点分量的个数
/// @param r1 : 开始控制点的下标
/// @param r2 : 结束控制点的下标
/// @param knots_vector : 节点矢量
/// @param control_points : 控制点
/// @param PK : out_put_param PK(i, k)表示k阶导曲线的第i + r1个控制点; remark:最差情况浪费一半内存
/// @return ENUM_NURBS错误码
template<typename T, int control_points_count, int degree, int n, int rows>
ENUM_NURBS curve_deriv_cpts(int r1, int r2, Eigen::Vector<T, control_points_count + degree + 1> const &knots_vector,
                    Eigen::Matrix<T, rows, control_points_count> const &control_points,
                    Eigen::Vector<Eigen::Matrix<T, rows, Eigen::Dynamic>, n + 1> &PK)
{
    int r = r2 - r1;
    // static_assert(r <= degree, "r2 - r1 > degree");
    // static_assert(n <= degree, "n > degree");
    PK[0] = control_points.block(0, r1, rows, r + 1);
    for (int k = 1; k <= n; ++k) //求导阶数
    {
        int temp = degree - k + 1;
        Eigen::Array<T, Eigen::Dynamic, 1> denominator = knots_vector.block(r1 + degree + 1, 0, r - k + 1, 1) - knots_vector.block(r1 + k, 0, r - k + 1, 1);
        PK[k].resize(rows, r + 1);
        PK[k].block(0, 0, rows, r - k + 1) = (temp * (PK[k - 1].block(0, 1, rows, r - k + 1) - PK[k - 1].block(0, 0, rows, r - k + 1))).array().rowwise() / denominator.transpose();
    }
    return  NURBS_SUCCESS;
}

/// @brief 计算nurbs曲线求导后的新控制点
/// @tparam T double float int...
/// @tparam control_points_count : 控制点个数
/// @tparam degree : nurbs曲线的阶数
/// @tparam rows : 坐标点分量的个数
/// @param n : 求到n次导
/// @param r1 : 开始控制点的下标
/// @param r2 : 结束控制点的下标
/// @param knots_vector : 节点矢量
/// @param control_points : 控制点
/// @param PK : out_put_param PK(i, k)表示k阶导曲线的第i + r1个控制点; remark:最差情况浪费一半内存
/// @return ENUM_NURBS错误码
template<typename T, int control_points_count, int degree, int rows>
ENUM_NURBS curve_deriv_cpts(int n, int r1, int r2, Eigen::Vector<T, control_points_count + degree + 1> const &knots_vector,
                    Eigen::Matrix<T, rows, control_points_count> const &control_points,
                    Eigen::VectorX<Eigen::Matrix<T, rows, Eigen::Dynamic>> &PK)
{
    int r = r2 - r1;
    // static_assert(r <= degree, "r2 - r1 > degree");
    // static_assert(n <= degree, "n > degree");
    PK[0] = control_points.block(0, r1, rows, r + 1);
    for (int k = 1; k <= n; ++k) //求导阶数
    {
        int temp = degree - k + 1;
        Eigen::Array<T, Eigen::Dynamic, 1> denominator = knots_vector.block(r1 + degree + 1, 0, r - k + 1, 1) - knots_vector.block(r1 + k, 0, r - k + 1, 1);
        PK[k].resize(rows, r + 1);
        PK[k].block(0, 0, rows, r - k + 1) = (temp * (PK[k - 1].block(0, 1, rows, r - k + 1) - PK[k - 1].block(0, 0, rows, r - k + 1))).array().rowwise() / denominator.transpose();
    }
    return  NURBS_SUCCESS;
}

/// @brief 计算nurbs曲线求导后的新控制点
/// @tparam T double float int...
/// @tparam degree : nurbs曲线的阶数
/// @tparam n : 求到n次导
/// @tparam rows : 坐标点分量的个数
/// @param r1 : 开始控制点的下标
/// @param r2 : 结束控制点的下标
/// @param knots_vector : 节点矢量
/// @param control_points : 控制点
/// @param PK : out_put_param PK(i, k)表示k阶导曲线的第i + r1个控制点; remark:最差情况浪费一半内存
/// @return ENUM_NURBS错误码
template<typename T, int degree, int n, int rows>
ENUM_NURBS curve_deriv_cpts(int r1, int r2, Eigen::VectorX<T> const &knots_vector,
                    Eigen::Matrix<T, rows, Eigen::Dynamic> const &control_points,
                    Eigen::Vector<Eigen::Matrix<T, rows, Eigen::Dynamic>, n + 1> &PK)
{
    int r = r2 - r1;
    // static_assert(r <= degree, "r2 - r1 > degree");
    // static_assert(n <= degree, "n > degree");
    PK[0] = control_points.block(0, r1, rows, r + 1);
    for (int k = 1; k <= n; ++k) //求导阶数
    {
        int temp = degree - k + 1;
        Eigen::Array<T, Eigen::Dynamic, 1> denominator = knots_vector.block(r1 + degree + 1, 0, r - k + 1, 1) - knots_vector.block(r1 + k, 0, r - k + 1, 1);
        PK[k].resize(rows, r + 1);
        PK[k].block(0, 0, rows, r - k + 1) = (temp * (PK[k - 1].block(0, 1, rows, r - k + 1) - PK[k - 1].block(0, 0, rows, r - k + 1))).array().rowwise() / denominator.transpose();
    }
    return  NURBS_SUCCESS;
}

/// @brief 计算nurbs曲线求导后的新控制点
/// @tparam T double float int...
/// @tparam degree : nurbs曲线的阶数
/// @tparam rows : 坐标点分量的个数
/// @param n : 求到n次导
/// @param r1 : 开始控制点的下标
/// @param r2 : 结束控制点的下标
/// @param knots_vector : 节点矢量
/// @param control_points : 控制点
/// @param PK : out_put_param PK(i, k)表示k阶导曲线的第i + r1个控制点; remark:最差情况浪费一半内存
/// @return ENUM_NURBS错误码
template<typename T, int degree, int rows>
ENUM_NURBS curve_deriv_cpts(int n, int r1, int r2, Eigen::Vector<T, Eigen::Dynamic> const &knots_vector,
                    Eigen::Matrix<T, rows, Eigen::Dynamic> const &control_points,
                    Eigen::VectorX<Eigen::Matrix<T, rows, Eigen::Dynamic>> &PK)
{
    int r = r2 - r1;
    // static_assert(r <= degree, "r2 - r1 > degree");
    // static_assert(n <= degree, "n > degree");
    PK[0] = control_points.block(0, r1, rows, r + 1);
    for (int k = 1; k <= n; ++k) //求导阶数
    {
        int temp = degree - k + 1;
        Eigen::Array<T, Eigen::Dynamic, 1> denominator = knots_vector.block(r1 + degree + 1, 0, r - k + 1, 1) - knots_vector.block(r1 + k, 0, r - k + 1, 1);
        PK[k].resize(rows, r + 1);
        PK[k].block(0, 0, rows, r - k + 1) = (temp * (PK[k - 1].block(0, 1, rows, r - k + 1) - PK[k - 1].block(0, 0, rows, r - k + 1))).array().rowwise() / denominator.transpose();
    }
    return  NURBS_SUCCESS;
}

/// @brief 计算nurbs曲线求导后的新控制点
/// @tparam T double float int...
/// @tparam n : 求到n次导
/// @tparam rows : 坐标点分量的个
/// @param degree : nurbs曲线的阶数数
/// @param r1 : 开始控制点的下标
/// @param r2 : 结束控制点的下标
/// @param knots_vector : 节点矢量
/// @param control_points : 控制点
/// @param PK : out_put_param PK(i, k)表示k阶导曲线的第i + r1个控制点; remark:最差情况浪费一半内存
/// @return ENUM_NURBS错误码
template<typename T, int n, int rows>
ENUM_NURBS curve_deriv_cpts(int degree, int r1, int r2, Eigen::VectorX<T> const &knots_vector,
                    Eigen::Matrix<T, rows, Eigen::Dynamic> const &control_points,
                    Eigen::Vector<Eigen::Matrix<T, rows, Eigen::Dynamic>, n + 1> &PK)
{
    int r = r2 - r1;
    // static_assert(r <= degree, "r2 - r1 > degree");
    // static_assert(n <= degree, "n > degree");
    PK[0] = control_points.block(0, r1, rows, r + 1);
    for (int k = 1; k <= n; ++k) //求导阶数
    {
        int temp = degree - k + 1;
        Eigen::Array<T, Eigen::Dynamic, 1> denominator = knots_vector.block(r1 + degree + 1, 0, r - k + 1, 1) - knots_vector.block(r1 + k, 0, r - k + 1, 1);
        PK[k].resize(rows, r + 1);
        PK[k].block(0, 0, rows, r - k + 1) = (temp * (PK[k - 1].block(0, 1, rows, r - k + 1) - PK[k - 1].block(0, 0, rows, r - k + 1))).array().rowwise() / denominator.transpose();
    }
    return  NURBS_SUCCESS;
}

/// @brief 计算nurbs曲线求导后的新控制点
/// @tparam T double float int...
/// @tparam rows : 坐标点分量的个数
/// @param degree : nurbs曲线的阶数
/// @param n : 求到n次导
/// @param r1 : 开始控制点的下标
/// @param r2 : 结束控制点的下标
/// @param knots_vector : 节点矢量
/// @param control_points : 控制点
/// @param PK : out_put_param PK(i, k)表示k阶导曲线的第i + r1个控制点; remark:最差情况浪费一半内存
/// @return ENUM_NURBS错误码
template<typename T, int rows>
ENUM_NURBS curve_deriv_cpts(int degree, int n, int r1, int r2, Eigen::Vector<T, Eigen::Dynamic> const &knots_vector,
                    Eigen::Matrix<T, rows, Eigen::Dynamic> const &control_points,
                    Eigen::VectorX<Eigen::Matrix<T, rows, Eigen::Dynamic>> &PK)
{
    int r = r2 - r1;
    // static_assert(r <= degree, "r2 - r1 > degree");
    // static_assert(n <= degree, "n > degree");
    PK[0] = control_points.block(0, r1, rows, r + 1);
    for (int k = 1; k <= n; ++k) //求导阶数
    {
        int temp = degree - k + 1;
        Eigen::Array<T, Eigen::Dynamic, 1> denominator = knots_vector.block(r1 + degree + 1, 0, r - k + 1, 1) - knots_vector.block(r1 + k, 0, r - k + 1, 1);
        PK[k].resize(rows, r + 1);
        PK[k].block(0, 0, rows, r - k + 1) = (temp * (PK[k - 1].block(0, 1, rows, r - k + 1) - PK[k - 1].block(0, 0, rows, r - k + 1))).array().rowwise() / denominator.transpose();
    }
    return  NURBS_SUCCESS;
}


// template<typename T, int point_size, int u_degree, int v_degree, int d>
// using surface_derv_cpts_matrix  =  Eigen::Matrix<
//                 Eigen::Matrix<
//                         Eigen::Matrix<T, point_size, Eigen::Dynamic, Eigen::ColMajor, point_size, u_degree + 1>,
//                 Eigen::Dynamic, 1, Eigen::ColMajor, v_degree + 1>,
//                 d + 1, d + 1>;

// template<typename T, int u_degree, int v_degree, int r1, int r2, int s1, int s2, int point_size, int d>
// ENUM_NURBS surface_deriv_cpts(int u_knots_size, int v_knots_size , Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> const & control_points,
//     Eigen::VectorX<T> const & u_knots_vector, Eigen::VectorX<T> const & v_knots_vector,
//     surface_derv_cpts_matrix<T, point_size, u_degree, v_degree, d> &PKL)
// {
//     int constexpr du = std::min(d, u_degree);
//     int constexpr dv = std::min(d, v_degree);
//     int constexpr r = r2 - r1;
//     int constexpr s = s2 - s1;
//     static_assert(du >= 0 && dv >= 0 && r >= 0 && s >= 0, "param is unvalid");
//     Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> PK;
//     for (int j = s1; j <= s2; ++j)
//     {
//         curve_deriv_cpts<T, u_degree, du, point_size, r1, r2>(u_knots_vector, control_points[j], PK);
//         for (int k = 0; k <= du; ++k)
//         {
//             PKL(k, 0)[j - s1] = PK[k];
//         }
//     }
//     auto const v_sub_knots_vector = v_knots_vector.block(s1, 0, v_knots_size - s1, 1);
//     for (int k = 0; k <= du; ++k)              //??????
//         for (int i = 0; i <= r - k; ++i)
//         {
//             Eigen::Matrix<T, point_size, s + 1> temp;
//             for (int index = 0; index <= s; ++index)
//                 temp.col(index) = PKL(k, 0)[index].col(i);
//             int dd = std::min(d - k, dv);
//             curve_deriv_cpts<T, v_degree, point_size, 0, s>(dd, v_sub_knots_vector, temp, PK);
//             for (int l = 1; l <= dd; ++l)
//                 for (int j = 0; j <= s - l; ++j)
//                     PKL(k, l)[j].col(i) = PK(l).col(j);
//         }
//     return ENUM_NURBS::NURBS_SUCCESS;
// }

/// @brief nurbs曲线的节点矢量插入一个几点r次
/// @tparam T double float int...
/// @tparam r 插入节点次数
/// @tparam point_size 控制点坐标分量的个数
/// @param control_points_count 控制点的个数
/// @param degree nurbs的阶数
/// @param knots_vector nurbs的节点矢量
/// @param old_control_points 插入节点之前的nurbs的控制点
/// @param u 插入节点
/// @param span u \in konts_vector[u_span, u_(span+1))
/// @param repeat_count 参数u在节点矢量中的重复度
/// @param new_control_points 插入节点之后的节点矢量
/// @return ENUM_NURBS错误码
template<typename T, int r, int point_size>
ENUM_NURBS curve_knots_insert(int control_points_count, int degree, const Eigen::VectorX<T> &knots_vector,
        const Eigen::Matrix<T, point_size, Eigen::Dynamic> &old_control_points, T u, int span, int repeat_count,
        Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points)
{
    new_control_points.resize(point_size, control_points_count + r);
    new_control_points.block(0, 0, point_size, span - degree + 1) = old_control_points.block(0, 0, point_size, span - degree + 1);
    new_control_points.block(0, span - repeat_count + r,point_size, control_points_count - span + repeat_count)
                    = old_control_points.block(0, span - repeat_count, point_size, control_points_count - span + repeat_count);
    Eigen::Matrix<T, point_size, Eigen::Dynamic> RW;
    RW.resize(point_size, degree - repeat_count + 1);
    RW = old_control_points.block(0, span - degree, point_size, degree - repeat_count + 1);
    //插入节点r次
    for (int j = 1; j <= r; ++j)
    {
        int L = span - degree + j;
        int count = degree - j - repeat_count;
        for (int i = 0; i <= count; ++i)
        {
            T alpha = (u - knots_vector[i + L]) / (knots_vector[i + span + 1] - knots_vector[i + L]);
            RW.col(i) = alpha * RW.col(i + 1) + (1.0 - alpha) * RW.col(i);
        }
        new_control_points.col(L) = RW.col(0);
        new_control_points.col(span + r - j - repeat_count) = RW.col(degree - j - repeat_count);
    }
    new_control_points.block(0, span - degree + 1 + r, point_size, degree - repeat_count - r) = RW.block(0, 1, point_size, degree - repeat_count - r);
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief nurbs曲线的节点矢量插入一个几点r次
/// @tparam T double float int...
/// @tparam point_size 控制点坐标分量的个数
/// @param r 插入节点次数
/// @param control_points_count 控制点的个数
/// @param degree nurbs的阶数
/// @param knots_vector nurbs的节点矢量
/// @param old_control_points 插入节点之前的nurbs的控制点
/// @param u 插入节点
/// @param span u \in konts_vector[u_span, u_(span+1))
/// @param repeat_count 参数u在节点矢量中的重复度
/// @param new_control_points 插入节点之后的节点矢量
/// @return ENUM_NURBS错误码
template<typename T, int point_size>
ENUM_NURBS curve_knots_insert(int r, int control_points_count, int degree, const Eigen::VectorX<T> &knots_vector,
        const Eigen::Matrix<T, point_size, Eigen::Dynamic> &old_control_points, T u, int span, int repeat_count,
        Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points)
{
    new_control_points.resize(point_size, control_points_count + r);
    new_control_points.block(0, 0, point_size, span - degree + 1) = old_control_points.block(0, 0, point_size, span - degree + 1);
    new_control_points.block(0, span - repeat_count + r,point_size, control_points_count - span + repeat_count)
                    = old_control_points.block(0, span - repeat_count, point_size, control_points_count - span + repeat_count);
    Eigen::Matrix<T, point_size, Eigen::Dynamic> RW;
    RW.resize(point_size, degree - repeat_count + 1);
    RW = old_control_points.block(0, span - degree, point_size, degree - repeat_count + 1);
    //插入节点r次
    for (int j = 1; j <= r; ++j)
    {
        int L = span - degree + j;
        int count = degree - j - repeat_count;
        for (int i = 0; i <= count; ++i)
        {
            T alpha = (u - knots_vector[i + L]) / (knots_vector[i + span + 1] - knots_vector[i + L]);
            RW.col(i) = alpha * RW.col(i + 1) + (1.0 - alpha) * RW.col(i);
        }
        new_control_points.col(L) = RW.col(0);
        new_control_points.col(span + r - j - repeat_count) = RW.col(degree - j - repeat_count);
    }
    new_control_points.block(0, span - degree + 1 + r, point_size, degree - repeat_count - r) = RW.block(0, 1, point_size, degree - repeat_count - r);
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief nurbs曲线的节点矢量插入一个几点r次
/// @tparam T double float int...
/// @tparam r 插入节点次数
/// @tparam point_size 控制点坐标分量的个数
/// @tparam degree nurbs的阶数
/// @tparam control_points_count 控制点的个数
/// @param knots_vector nurbs的节点矢量
/// @param old_control_points 插入节点之前的nurbs的控制点
/// @param u 插入节点
/// @param span u \in konts_vector[u_span, u_(span+1))
/// @param repeat_count 参数u在节点矢量中的重复度
/// @param new_control_points 插入节点之后的节点矢量
/// @return ENUM_NURBS错误码
template<typename T, int r, int point_size, int degree, int control_points_count>
ENUM_NURBS curve_knots_insert(const Eigen::Vector<T, degree + 1 + control_points_count> &knots_vector,
        const Eigen::Matrix<T, point_size, control_points_count> &old_control_points, T u, int span, int repeat_count,
        Eigen::Matrix<T, point_size, control_points_count + r> &new_control_points)
{
    // new_control_points.resize(point_size, control_points_count + r);
    new_control_points.block(0, 0, point_size, span - degree + 1) = old_control_points.block(0, 0, point_size, span - degree + 1);
    new_control_points.block(0, span - repeat_count + r,point_size, control_points_count - span + repeat_count)
                    = old_control_points.block(0, span - repeat_count, point_size, control_points_count - span + repeat_count);
    Eigen::Matrix<T, point_size, Eigen::Dynamic> RW;
    RW.resize(point_size, degree - repeat_count + 1);
    RW = old_control_points.block(0, span - degree, point_size, degree - repeat_count + 1);
    //插入节点r次
    for (int j = 1; j <= r; ++j)
    {
        int L = span - degree + j;
        int count = degree - j - repeat_count;
        for (int i = 0; i <= count; ++i)
        {
            T alpha = (u - knots_vector[i + L]) / (knots_vector[i + span + 1] - knots_vector[i + L]);
            RW.col(i) = alpha * RW.col(i + 1) + (1.0 - alpha) * RW.col(i);
        }
        new_control_points.col(L) = RW.col(0);
        new_control_points.col(span + r - j - repeat_count) = RW.col(degree - j - repeat_count);
    }
    new_control_points.block(0, span - degree + 1 + r, point_size, degree - repeat_count - r) = RW.block(0, 1, point_size, degree - repeat_count - r);
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief 通过将节点插入r次来计算此参数对应的曲线上的点
/// @tparam T double float int...
/// @tparam point_size 控制点坐标分量的个数
/// @tparam degree 曲线的阶数
/// @tparam control_point_count 控制点的个数
/// @param u 参数u
/// @param knots_vector 节点矢量
/// @param control_point 控制点
/// @param point out_put_param 参数u对应的点
/// @return ENUM_NURBS错误码
template<typename T, int point_size, int degree, int control_point_count>
ENUM_NURBS curve_pnt_by_corner_cut(T u,const Eigen::Vector<T, degree + control_point_count + 1> &knots_vector,
    const Eigen::Matrix<T, point_size, control_point_count> &control_point, Eigen::Vector<T, point_size> &point)
{
    // if (std::abs(u - knots_vector[0]) < KNOTS_VECTOR_EPS)
    if (u == knots_vector[0])
    {
        point = control_point.col(0);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    // if (std::abs(u - knots_vector[control_point_count + degree]) < KNOTS_VECTOR_EPS)
    if (u == knots_vector[control_point_count + degree])
    {
        point = control_point.col(control_point_count - 1);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    int index = -1;
    find_span<T, control_point_count, degree>(u, knots_vector, index);
    int repeat = 0;
    int current_index = index;
    while (current_index >= 0)
    {
        // if (std::abs(knots_vector[current_index--] - u) < KNOTS_VECTOR_EPS)
        if (knots_vector[current_index--] == u)
        {
            ++repeat;
            continue;
        }
        break;
    }
    int r = degree - repeat;
    Eigen::Matrix<T, point_size, Eigen::Dynamic> RW;
    RW.resize(point_size, r + 1);
    RW = control_point.block(0, index - degree, point_size, r + 1);
    for (int j = 1; j <= r; ++j)
    {
        for (int i = 0; i <= r - j; ++i)
        {
            T alpha = (u - knots_vector[index - degree + j + i]) / (knots_vector[u + index + 1] - knots_vector[index - degree + j + i]);
            RW.col(i) = alpha * RW.col(i + 1) + (1.0 - alpha) * RW.col(i);
        }
    }
    point = RW.col(0);
    return ENUM_NURBS::NURBS_SUCCESS;

}


/// @brief 通过将节点插入r次来计算此参数对应的曲线上的点
/// @tparam T double float int...
/// @tparam point_size 控制点坐标分量的个数
/// @tparam degree 曲线的阶数
/// @param u 参数u
/// @param control_point_count 控制点的个数
/// @param knots_vector 节点矢量
/// @param control_point 控制点
/// @param point out_put_param 参数u对应的点
/// @return ENUM_NURBS错误码
template<typename T, int point_size, int degree>
ENUM_NURBS curve_pnt_by_corner_cut(T u, int control_point_count, const Eigen::VectorX<T> &knots_vector,
    const Eigen::Matrix<T, point_size, Eigen::Dynamic> &control_point, Eigen::Vector<T, point_size> &point)
{
    // if (std::abs(u - knots_vector[0]) < KNOTS_VECTOR_EPS)
    if (u == knots_vector[0])
    {
        point = control_point.col(0);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    // if (std::abs(u - knots_vector[control_point_count + degree]) < KNOTS_VECTOR_EPS)
    if (u == knots_vector[control_point_count + degree])
    {
        point = control_point.col(control_point_count - 1);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    int index = -1;
    find_span<T, degree>(u, knots_vector, index);
    int repeat = 0;
    int current_index = index;
    while (current_index >= 0)
    {
        // if (std::abs(knots_vector[current_index--] - u) < KNOTS_VECTOR_EPS)
        if (knots_vector[current_index--] == u)
        {
            ++repeat;
            continue;
        }
        break;
    }
    int r = degree - repeat;
    Eigen::Matrix<T, point_size, Eigen::Dynamic> RW;
    RW.resize(point_size, r + 1);
    RW = control_point.block(0, index - degree, point_size, r + 1);
    for (int j = 1; j <= r; ++j)
    {
        for (int i = 0; i <= r - j; ++i)
        {
            T alpha = (u - knots_vector[index - degree + j + i]) / (knots_vector[u + index + 1] - knots_vector[index - degree + j + i]);
            RW.col(i) = alpha * RW.col(i + 1) + (1.0 - alpha) * RW.col(i);
        }
    }
    point = RW.col(0);
    return ENUM_NURBS::NURBS_SUCCESS;

}

/// @brief 通过将节点插入r次来计算此参数对应的曲线上的点
/// @tparam T double float int...
/// @tparam point_size 控制点坐标分量的个数
/// @param u 参数u
/// @param control_point_count 控制点的个数
/// @param degree 曲线的阶数
/// @param knots_vector 节点矢量
/// @param control_point 控制点
/// @param point out_put_param 参数u对应的点
/// @return ENUM_NURBS错误码
template<typename T, int point_size>
ENUM_NURBS curve_pnt_by_corner_cut(T u, int control_point_count, int degree, const Eigen::VectorX<T> &knots_vector,
    const Eigen::Matrix<T, point_size, Eigen::Dynamic> &control_point, Eigen::Vector<T, point_size> &point)
{
    // if (std::abs(u - knots_vector[0]) < KNOTS_VECTOR_EPS)
    if (u == knots_vector[0])
    {
        point = control_point.col(0);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    // if (std::abs(u - knots_vector[control_point_count + degree]) < KNOTS_VECTOR_EPS)
    if (u == knots_vector[control_point_count + degree])
    {
        point = control_point.col(control_point_count - 1);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    int index = -1;
    find_span<T>(u, degree, knots_vector, index);
    int repeat = 0;
    int current_index = index;
    while (current_index >= 0)
    {
        // if (std::abs(knots_vector[current_index--] - u) < KNOTS_VECTOR_EPS)
        if (knots_vector[current_index--] == u)
        {
            ++repeat;
            continue;
        }
        break;
    }
    int r = degree - repeat;
    Eigen::Matrix<T, point_size, Eigen::Dynamic> RW;
    RW.resize(point_size, r + 1);
    RW = control_point.block(0, index - degree, point_size, r + 1);
    for (int j = 1; j <= r; ++j)
    {
        for (int i = 0; i <= r - j; ++i)
        {
            T alpha = (u - knots_vector[index - degree + j + i]) / (knots_vector[u + index + 1] - knots_vector[index - degree + j + i]);
            RW.col(i) = alpha * RW.col(i + 1) + (1.0 - alpha) * RW.col(i);
        }
    }
    point = RW.col(0);
    return ENUM_NURBS::NURBS_SUCCESS;

}


/// @brief 节点加细
/// @tparam T double float
/// @tparam point_size 控制点坐标分量的个数
/// @param knots_vector_size 节点矢量的个数
/// @param degree 曲线阶数
/// @param knots_vector 节点矢量
/// @param control_points 控制点
/// @param insert_knots 插入节点
/// @param new_knots_vector 新的节点矢量
/// @param new_control_points 新的控制点
/// @return ENUM_NURBS错误码
template<typename T, int point_size>
ENUM_NURBS refine_knots_vector_curve(int knots_vector_size, int degree, const Eigen::VectorX<T> &knots_vector,
    const Eigen::Matrix<T, point_size, Eigen::Dynamic> &control_points, const Eigen::VectorX<T> &insert_knots,
    Eigen::VectorX<T> &new_knots_vector, Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points)
{
    int start_index = -1, end_index = -1;
    if (find_span<T>(insert_knots[0], degree, knots_vector, start_index) != ENUM_NURBS::NURBS_SUCCESS)
        return ENUM_NURBS::NURBS_ERROR;
    int insert_knots_size = insert_knots.size();
    if (find_span<T>(insert_knots[insert_knots_size - 1], degree, knots_vector, end_index) != ENUM_NURBS::NURBS_SUCCESS)
        return ENUM_NURBS::NURBS_ERROR;
    end_index += 1;

    new_knots_vector.resize(knots_vector_size + insert_knots_size);
    int control_points_size = control_points.cols();
    new_control_points.resize(point_size, control_points_size + insert_knots_size);
    new_control_points.block(0, 0, point_size, start_index - degree + 1) = control_points.block(0, 0, point_size, start_index - degree + 1);
    new_control_points.block(0, end_index + insert_knots_size - 1, point_size, control_points_size - end_index + 1) =
        control_points.block(0, end_index - 1, point_size, control_points_size - end_index + 1);

    new_knots_vector.block(0, 0, start_index + 1, 1) = knots_vector.block(0, 0, start_index + 1, 1);
    new_knots_vector.block(end_index + degree + insert_knots_size, 0, control_points_size - end_index + 1, 1) = knots_vector.block(end_index + degree, 0, control_points_size - end_index + 1, 1);

    int i = end_index + degree - 1, k = end_index + insert_knots_size + degree - 1;
    for (int j = insert_knots_size - 1; j >= 0; --j)
    {
        while (insert_knots[j] <= knots_vector[i] && i > start_index)
        {
            new_control_points.col(k - degree - 1) = control_points.col(i - degree - 1);
            new_knots_vector[k] = knots_vector[i];
            --k; --i;
        }
        new_control_points.col(k - degree - 1) = new_control_points.col(k - degree);
        for (int l = 1; l <= degree; ++l)
        {
            int ind = k - degree + l;
            T alpha = new_knots_vector[k + l] - insert_knots[j];
            // if (std::abs(alpha) < KNOTS_VECTOR_EPS)
            if (alpha == 0.0)
                new_control_points.col(ind - 1) = new_control_points.col(ind);
            else
            {
                alpha /= (new_knots_vector[k + l] - knots_vector[i - degree + l]);
                new_control_points.col(ind - 1) = alpha * new_control_points.col(ind - 1) + (1.0 - alpha) * new_control_points.col(ind);
            }
        }
        new_knots_vector[k] = insert_knots[j];
        --k;
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief 节点加细
/// @tparam T double float
/// @tparam point_size 控制点坐标分量的个数
/// @param knots_vector_size 节点矢量的个数
/// @param degree 曲线阶数
/// @param knots_vector 节点矢量
/// @param control_points 控制点
/// @param insert_knots 插入节点
/// @param new_knots_vector 新的节点矢量
/// @param new_control_points 新的控制点
/// @return ENUM_NURBS错误码
template<typename T, int point_size, int knots_vector_size, int degree>
ENUM_NURBS refine_knots_vector_curve(const Eigen::Vector<T, knots_vector_size> &knots_vector,
    const Eigen::Matrix<T, point_size, knots_vector_size - degree - 1> &control_points, const Eigen::VectorX<T> &insert_knots,
    Eigen::VectorX<T> &new_knots_vector, Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points)
{
    int start_index = -1, end_index = -1;
    if (find_span<T, knots_vector_size - degree - 1, degree>(insert_knots[0], knots_vector, start_index) != ENUM_NURBS::NURBS_SUCCESS)
        return ENUM_NURBS::NURBS_ERROR;
    int insert_knots_size = insert_knots.size();
    if (find_span<T, knots_vector_size - degree - 1, degree>(insert_knots[insert_knots_size - 1], knots_vector, end_index) != ENUM_NURBS::NURBS_SUCCESS)
        return ENUM_NURBS::NURBS_ERROR;
    end_index += 1;

    new_knots_vector.resize(knots_vector_size + insert_knots_size);
    constexpr int control_points_size = knots_vector_size - degree - 1;
    new_control_points.resize(point_size, control_points_size + insert_knots_size);
    new_control_points.block(0, 0, point_size, start_index - degree + 1) = control_points.block(0, 0, point_size, start_index - degree + 1);
    new_control_points.block(0, end_index + insert_knots_size - 1, point_size, control_points_size - end_index + 1) =
        control_points.block(0, end_index - 1, point_size, control_points_size - end_index + 1);

    new_knots_vector.block(0, 0, start_index + 1, 1) = knots_vector.block(0, 0, start_index + 1, 1);
    new_knots_vector.block(end_index + degree + insert_knots_size, 0, control_points_size - end_index + 1, 1) = knots_vector.block(end_index + degree, 0, control_points_size - end_index + 1, 1);

    int i = end_index + degree - 1, k = end_index + insert_knots_size + degree - 1;
    for (int j = insert_knots_size - 1; j >= 0; --j)
    {
        while (insert_knots[j] <= knots_vector[i] && i > start_index)
        {
            new_control_points.col(k - degree - 1) = control_points.col(i - degree - 1);
            new_knots_vector[k] = knots_vector[i];
            --k; --i;
        }
        new_control_points.col(k - degree - 1) = new_control_points.col(k - degree);
        for (int l = 1; l <= degree; ++l)
        {
            int ind = k - degree + l;
            T alpha = new_knots_vector[k + l] - insert_knots[j];
            // if (std::abs(alpha) < KNOTS_VECTOR_EPS)
            if (alpha == 0.0)
                new_control_points.col(ind - 1) = new_control_points.col(ind);
            else
            {
                alpha /= (new_knots_vector[k + l] - knots_vector[i - degree + l]);
                new_control_points.col(ind - 1) = alpha * new_control_points.col(ind - 1) + (1.0 - alpha) * new_control_points.col(ind);
            }
        }
        new_knots_vector[k] = insert_knots[j];
        --k;
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief 计算将nurbs曲线分解成bezier的控制点
/// @tparam T double float int...
/// @tparam point_size 控制点坐标分量的个数
/// @param degree 曲线阶数
/// @param interval_segment nurbs节点矢量区间的个数(不算0长度区间)
/// @param knots_vector 节点矢量
/// @param control_points 控制点
/// @param new_knots_vector out_put_param 新的节点矢量, 要确保new_control_points的行数是正确的
/// @param new_control_points out_put_param 新的控制点
/// @return ENUM_NURBS错误码
template<typename T, int point_size>
ENUM_NURBS decompose_curve_to_bezier(int degree, int interval_segment, const Eigen::VectorX<T> &knots_vector, const Eigen::Matrix<T, point_size, Eigen::Dynamic> &control_points,
    Eigen::VectorX<T> &new_knots_vector, Eigen::VectorX<Eigen::Matrix<T, point_size, Eigen::Dynamic>> &new_control_points)
{
    for (int index = 0; index < interval_segment; ++index)
        new_control_points[index].resize(point_size, degree + 1);

    // new_control_points[0].block(0, 0, point_size, degree + 1) = control_points.block(0, 0, point_size, degree + 1);
    new_control_points[0] = control_points.block(0, 0, point_size, degree + 1);
    int b = degree + 1;
    int a = degree;
    int nb = 0;
    // int knots_vector_size = knots_vector.size();
    int konts_end_index = knots_vector.size() - 1;
    while (b < konts_end_index)
    {
        int i = b;
        while (b < konts_end_index && knots_vector[b + 1] == knots_vector[b])
        {
            ++b;
        }
        // new_control_points[nb + 1].resize(point_size, degree + 1);
        int mult = b - i + 1;
        if (mult < degree)
        {
            T numer = knots_vector[b] - knots_vector[a];
            Eigen::VectorX<T> alphas(degree - mult);
            for (int j = degree; j > mult; --j)
                alphas[j - mult - 1] = numer / (knots_vector[a + j] - knots_vector[a]);
            // alphas.setConstant(numer);
            // alphas = alphas.array() / (knots_vector.block(a + mult + 1, 0, degree - mult, 1).rowwise() - knots_vector[a]).array();
            int r = degree - mult;
            for (int j = 1; j <= r; ++j)
            {
                int s = mult + j;
                for (int k = degree; k >= s; --k)
                {
                    T alpha = alphas[k - s];
                    new_control_points[nb].col(k) = alpha * new_control_points[nb].col(k) + (1.0 - alpha) * new_control_points[nb].col(k - 1);
                    // std::cout << new_control_points[nb].col(k) << std::endl;
                }
                if (b < konts_end_index)
                {
                    new_control_points[nb + 1].col(r - j) = new_control_points[nb].col(degree);
                    // std::cout << new_control_points[nb + 1].col(r - j) << std::endl;
                }
            }
        }
        nb += 1;

        if (b < konts_end_index)
        {
            new_control_points[nb].block(0, degree - mult, point_size, mult + 1) = control_points.block(0, b - mult, point_size, mult + 1);
            // std::cout << new_control_points[nb].block(0, degree - mult, point_size, mult + 1) << std::endl;
            a = b;
            b += 1;
        }
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief 计算将nurbs曲线分解成bezier的控制点
/// @tparam T double float int...
/// @tparam point_size 控制点坐标分量的个数
/// @tparam degree 曲线阶数
/// @param interval_segment nurbs节点矢量区间的个数(不算0长度区间)
/// @param knots_vector 节点矢量
/// @param control_points 控制点
/// @param new_knots_vector out_put_param 新的节点矢量, 要确保new_control_points的行数是正确的
/// @param new_control_points out_put_param 新的控制点
/// @return ENUM_NURBS错误码
template<typename T, int point_size, int degree>
ENUM_NURBS decompose_curve_to_bezier(int interval_segment, const Eigen::VectorX<T> &knots_vector, const Eigen::Matrix<T, point_size, Eigen::Dynamic> &control_points,
    Eigen::VectorX<T> &new_knots_vector, Eigen::VectorX<Eigen::Matrix<T, point_size, degree + 1>> &new_control_points)
{
    // new_control_points[0].block(0, 0, point_size, degree + 1) = control_points.block(0, 0, point_size, degree + 1);
    new_control_points[0] = control_points.block(0, 0, point_size, degree + 1);
    int b = degree + 1;
    int a = degree;
    int nb = 0;
    int konts_end_index = knots_vector.size() - 1;
    while (b < konts_end_index)
    {
        int i = b;
        while (b < konts_end_index && knots_vector[b + 1] == knots_vector[b])
        {
            ++b;
        }
        // new_control_points[nb + 1].resize(point_size, degree + 1);
        int mult = b - i + 1;
        if (mult < degree)
        {
            T numer = knots_vector[b] - knots_vector[a];
            Eigen::VectorX<T> alphas(degree - mult);
            for (int j = degree; j > mult; --j)
                alphas[j - mult - 1] = numer / (knots_vector[a + j] - knots_vector[a]);
            // alphas.setConstant(numer);
            // alphas = alphas.array() / (knots_vector.block(a + mult + 1, 0, degree - mult, 1) - knots_vector.block(a, 0, degree - mult, 1)).array();
            int r = degree - mult;
            for (int j = 1; j <= r; ++j)
            {
                int s = mult + j;
                for (int k = degree; k >= s; --k)
                {
                    T alpha = alphas[k - s];
                    new_control_points[nb].col(k) = alpha * new_control_points[nb].col(k) + (1.0 - alpha) * new_control_points[nb].col(k - 1);
                }
                if (b < konts_end_index)
                {
                    new_control_points[nb + 1].col(r - j) = new_control_points[nb].col(degree);
                }
            }
        }
        nb += 1;

        if (b < konts_end_index)
        {
            new_control_points[nb].block(0, degree - mult, point_size, mult + 1) = control_points.block(0, b - mult, point_size, mult + 1);
            a = b;
            b += 1;
        }
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}

/// @brief 计算将nurbs曲线分解成bezier的控制点
/// @tparam T double float int...
/// @tparam point_size 控制点坐标分量的个数
/// @tparam point_count 控制点个数
/// @tparam degree 曲线阶数
/// @param interval_segment nurbs节点矢量区间的个数(不算0长度区间)
/// @param knots_vector 节点矢量
/// @param control_points 控制点
/// @param new_knots_vector out_put_param 新的节点矢量, 要确保new_control_points的行数是正确的
/// @param new_control_points out_put_param 新的控制点
/// @return ENUM_NURBS错误码
template<typename T, int point_size, int point_count, int degree>
ENUM_NURBS decompose_curve_to_bezier(int interval_segment, const Eigen::Vector<T, degree + point_count + 1> &knots_vector, const Eigen::Matrix<T, point_size, point_count> &control_points,
    Eigen::VectorX<T> &new_knots_vector, Eigen::VectorX<Eigen::Matrix<T, point_size, degree + 1>> &new_control_points)
{
    new_control_points[0] = control_points.block(0, 0, point_size, degree + 1);
    int b = degree + 1;
    int a = degree;
    int nb = 0;
    int konts_end_index = knots_vector.size() - 1;
    while (b < konts_end_index)
    {
        int i = b;
        while (b < konts_end_index && knots_vector[b + 1] == knots_vector[b])
        {
            ++b;
        }
        int mult = b - i + 1;
        if (mult < degree)
        {
            T numer = knots_vector[b] - knots_vector[a];
            Eigen::VectorX<T> alphas(degree - mult);
            for (int j = degree; j > mult; --j)
                alphas[j - mult - 1] = numer / (knots_vector[a + j] - knots_vector[a]);
            // alphas.setConstant(numer);
            // alphas = alphas.array() / (knots_vector.block(a + mult + 1, 0, degree - mult, 1) - knots_vector.block(a, 0, degree - mult, 1)).array();
            int r = degree - mult;
            for (int j = 1; j <= r; ++j)
            {
                int s = mult + j;
                for (int k = degree; k >= s; --k)
                {
                    T alpha = alphas[k - s];
                    new_control_points[nb].col(k) = alpha * new_control_points[nb].col(k) + (1.0 - alpha) * new_control_points[nb].col(k - 1);
                }
                if (b < konts_end_index)
                {
                    new_control_points[nb + 1].col(r - j) = new_control_points[nb].col(degree);
                }
            }
        }
        nb += 1;

        if (b < konts_end_index)
        {
            new_control_points[nb].block(0, degree - mult, point_size, mult + 1) = control_points.block(0, b - mult, point_size, mult + 1);
            a = b;
            b += 1;
        }
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}


template<typename T>
ENUM_NURBS find_interval_segment_count(int degree, const Eigen::VectorX<T> &knots_vector, int &count)
{
    int knots_count = knots_vector.size();
    count = 0;
    knots_count -= degree;
    T current_knots = knots_vector[0];
    for (int index = degree + 1; index < knots_count; ++index)
    {
        if (current_knots != knots_vector[index])
        {
            current_knots = knots_vector[index];
            count += 1;
        }
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}

template<typename T, int point_count, int degree>
ENUM_NURBS find_interval_segment_count(const Eigen::Vector<T, point_count + degree + 1> &knots_vector, int &count)
{
    constexpr int knots_count = point_count + 1;
    count = 0;
    T current_knots = knots_vector[0];
    for (int index = degree + 1; index < knots_count; ++index)
    {
        if (current_knots != knots_vector[index])
        {
            current_knots = knots_vector[index];
            count += 1;
        }
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief 去除nurbs曲线的节点
/// @tparam T double float int...
/// @tparam point_size 控制点坐标分量的个数
/// @tparam is_rational 是否是有理曲线
/// @param index 去除节点在节点矢量列表中的下指标
/// @param count 期望去除的次数
/// @param degree nurbs曲线的阶数
/// @param repeat 去除节点的重复度
/// @param knots_vector 节点矢量
/// @param control_points 控制点
/// @param time 实际去除的次数
/// @param error 节点出去误差
/// @return ENUM_NURBS错误码
template<typename T, int point_size, bool is_rational>
ENUM_NURBS remove_curve_knots(int index, int count, int degree, int repeat, Eigen::VectorX<T> &knots_vector,
     Eigen::Matrix<T, point_size, Eigen::Dynamic> &control_points , int &time, T error = DEFAULT_ERROR)
{
    T TOL = -1;
    int cols = control_points.cols();
    if (is_rational == false)
    {
        TOL = error;
    }
    else
    {
        T w_min = control_points(3, 0);
        if (w_min <= 0)
            return ENUM_NURBS::NUBRS_WIEGHT_IS_NONPOSITIVE;
        T p_max = (control_points.block(0, 0, 3, 1) / control_points(3, 0)).norm();
        for (int col = 1; col < cols; ++col)
        {
            if (w_min > control_points(3, col))
                w_min = control_points(3, col);
            if (w_min <= 0)
                return ENUM_NURBS::NUBRS_WIEGHT_IS_NONPOSITIVE;
            T temp = (control_points.block(0, col, 3, 1) / control_points(3, col)).norm();
            if (p_max < temp)
                p_max = temp;
        }
        TOL = (error * w_min) / (1 + p_max);
    }
    T u = knots_vector[index];
    int order = degree + 1;
    int first = index - degree;
    int last = index - repeat;
    std::vector<Eigen::Vector<T, point_size>> temp(2 * degree + 1);
    time = 0;
    for (; time < count; ++time)
    {
        int offset = first - 1;
        temp[0] = control_points.col(offset);
        temp[last + 1 - offset] = control_points.col(last + 1);
        int i = first, j = last;
        int ii = 1, jj = last - offset;
        bool remflag = false;
        while (j - i > time)
        {
            T alfi = (u - knots_vector[i]) / (knots_vector[i + order + time] - knots_vector[i]);
            T alfj = (u - knots_vector[j - time]) / (knots_vector[j + order] - knots_vector[j - time]);
            temp[ii] = (control_points.col(i) - (1.0 - alfi) * temp[ii - 1]) / alfi;
            temp[jj] = (control_points.col(j) - alfj * temp[jj + 1]) / (1.0 - alfj);
            i += 1;
            ii += 1;
            j -= 1;
            jj -= 1;
        }
        if (j - i < time)
        {
            if ((temp[ii - 1] - temp[jj + 1]).norm() <= TOL)
                remflag = true;
        }
        else
        {
            T alfi = (u - knots_vector[i]) / (knots_vector[i + order + time] - knots_vector[i]);
            Eigen::Vector<T, point_size> temp_vec = control_points.col(i) - alfi * temp[ii + time + 1] - (1.0 - alfi) * temp[ii - 1];
            if (temp_vec.norm() <= TOL)
                remflag = true;
        }
        if (remflag == false)
            break;
        else
        {
            int i = first, j = last;
            while (j - i > time)
            {
                control_points.col(i) = temp[i - offset];
                control_points.col(j) = temp[j - offset];
                i += 1;
                j -= 1;
            }
        }
        first -= 1; last += 1;
    }

    if (time == 0)
        return ENUM_NURBS::NURBS_SUCCESS;

    int knots_size = knots_vector.size();
    Eigen::VectorX<T> new_knots_vector(knots_size - time);
    new_knots_vector.block(0, 0, index + 1 - time, 1) = knots_vector.block(0, 0, index + 1 - time, 1);
    new_knots_vector.block(index + 1 - time, 0, knots_size - index - 1, 1) = knots_vector.block(index + 1, 0, knots_size - index - 1, 1);
    knots_vector = new_knots_vector;

    int j = (2 * index - repeat - degree) / 2, i = j;
    for (int k = 1; k < time; ++k)
    {
        if (k % 2 == 1)
            i += 1;
        else
            j -= 1;
    }
    int ij = i + 1 - j;
    for (int k = i + 1; k < cols; ++k)
    {
        control_points.col(j) = control_points.col(k);
        j += 1;
    }

    Eigen::Matrix<T, point_size, Eigen::Dynamic> temp_points = control_points.block(0, 0, point_size, cols - ij);
    control_points = temp_points;
    return ENUM_NURBS::NURBS_SUCCESS;
}


/// @brief 去除nurbs曲线的节点
/// @tparam T double float int...
/// @tparam degree 曲线阶数
/// @tparam point_count 控制点个数
/// @tparam point_size 控制点坐标分量的个数
/// @tparam is_rational 是否是有理曲线
/// @param index 去除节点在节点矢量列表中的下指标
/// @param count 期望去除的次数
/// @param repeat 去除节点的重复度
/// @param knots_vector 节点矢量
/// @param control_points 控制点
/// @param new_knots_vector out_put_param 新的节点矢量
/// @param new_control_points out_put_param 新的控制点
/// @param time 实际去除的次数
/// @param error 节点出去误差
/// @return ENUM_NURBS错误码
template<typename T, int degree, int point_count, int point_size, bool is_rational>
ENUM_NURBS remove_curve_knots(int index, int count, int repeat, const Eigen::Vector<T, point_count + degree + 1> &knots_vector,
     Eigen::Matrix<T, point_size, point_count> control_points , Eigen::VectorX<T> &new_knots_vector,
     Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points ,int &time, T error = DEFAULT_ERROR)
{
    T TOL = -1;
    int cols = control_points.cols();
    if (is_rational == false)
    {
        TOL = error;
    }
    else
    {
        T w_min = control_points(3, 0);
        if (w_min <= 0)
            return ENUM_NURBS::NUBRS_WIEGHT_IS_NONPOSITIVE;
        T p_max = (control_points.block(0, 0, 3, 1) / control_points(3, 0)).norm();
        for (int col = 1; col < cols; ++col)
        {
            if (w_min > control_points(3, col))
                w_min = control_points(3, col);
            if (w_min <= 0)
                return ENUM_NURBS::NUBRS_WIEGHT_IS_NONPOSITIVE;
            T temp = (control_points.block(0, col, 3, 1) / control_points(3, col)).norm();
            if (p_max < temp)
                p_max = temp;
        }
        TOL = (error * w_min) / (1 + p_max);
    }
    T u = knots_vector[index];
    int order = degree + 1;
    int first = index - degree;
    int last = index - repeat;
    std::vector<Eigen::Vector<T, point_size>> temp(2 * degree + 1);
    time = 0;
    for (; time < count; ++time)
    {
        int offset = first - 1;
        temp[0] = control_points.col(offset);
        temp[last + 1 - offset] = control_points.col(last + 1);
        int i = first, j = last;
        int ii = 1, jj = last - offset;
        bool remflag = false;
        while (j - i > time)
        {
            T alfi = (u - knots_vector[i]) / (knots_vector[i + order + time] - knots_vector[i]);
            T alfj = (u - knots_vector[j - time]) / (knots_vector[j + order] - knots_vector[j - time]);
            temp[ii] = (control_points.col(i) - (1.0 - alfi) * temp[ii - 1]) / alfi;
            temp[jj] = (control_points.col(j) - alfj * temp[jj + 1]) / (1.0 - alfj);
            i += 1;
            ii += 1;
            j -= 1;
            jj -= 1;
        }
        if (j - i < time)
        {
            if ((temp[ii - 1] - temp[jj + 1]).norm() <= TOL)
                remflag = true;
        }
        else
        {
            T alfi = (u - knots_vector[i]) / (knots_vector[i + order + time] - knots_vector[i]);
            Eigen::Vector<T, point_size> temp_vec = control_points.col(i) - alfi * temp[ii + time + 1] - (1.0 - alfi) * temp[ii - 1];
            if (temp_vec.norm() <= TOL)
                remflag = true;
        }
        if (remflag == false)
            break;
        else
        {
            int i = first, j = last;
            while (j - i > time)
            {
                control_points.col(i) = temp[i - offset];
                control_points.col(j) = temp[j - offset];
                i += 1;
                j -= 1;
            }
        }
        first -= 1; last += 1;
    }

    if (time == 0)
    {
        new_control_points = control_points;
        new_knots_vector = knots_vector;
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    int knots_size = knots_vector.size();
    // Eigen::VectorX<T> new_knots_vector(knots_size - time);
    new_knots_vector.resize(knots_size - time);
    new_knots_vector.block(0, 0, index + 1 - time, 1) = knots_vector.block(0, 0, index + 1 - time, 1);
    new_knots_vector.block(index + 1 - time, 0, knots_size - index - 1, 1) = knots_vector.block(index + 1, 0, knots_size - index - 1, 1);
    // knots_vector = new_knots_vector;

    int j = (2 * index - repeat - degree) / 2, i = j;
    for (int k = 1; k < time; ++k)
    {
        if (k % 2 == 1)
            i += 1;
        else
            j -= 1;
    }
    int ij = i + 1 - j;
    for (int k = i + 1; k < cols; ++k)
    {
        control_points.col(j) = control_points.col(k);
        j += 1;
    }

    new_control_points = control_points.block(0, 0, point_size, cols - ij);
    return ENUM_NURBS::NURBS_SUCCESS;
}



/// @brief 升阶
/// @tparam T double float int...
/// @tparam point_size 控制点坐标分量的个数
/// @param t 升阶次数
/// @param degree nurbs曲线阶数
/// @param knots_vector 节点向量
/// @param control_points 控制点
/// @param new_knots_vector 新的节点向量
/// @param new_control_points 新的控制点
/// @return ENUM_NURBS错误码
template<typename T, int point_size>
ENUM_NURBS degree_elevate_curve(int t, int degree, const Eigen::VectorX<T> &knots_vector, const Eigen::Matrix<T, point_size, Eigen::Dynamic> &control_points,
        Eigen::VectorX<T> &new_knots_vector, Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points)
{
    int knots_size = knots_vector.size();
    int ph = degree + t;
    int ph2 = ph / 2;
    Eigen::MatrixX<T> bezalfs(ph + 1, degree + 1);
    bezalfs(0, 0) = bezalfs(ph, degree) = 1.0;
    Eigen::MatrixX<int> Bin = binary_coeff(degree + t + 1);
    Eigen::Matrix<T, point_size, Eigen::Dynamic> next_bpts(point_size, degree - 1);
    for (int i = 1; i <= ph2; ++i)
    {
        T inv = 1.0 / Bin(ph, i);
        int mpi = std::min(degree, i);
        for (int j = std::max(0, i - t); j <= mpi; ++j)
            bezalfs(i, j) = inv * Bin(degree, j) * Bin(t, i - j);
    }
    for (int i = ph2 + 1; i < ph; ++i)
    {
        int mpi = std::min(degree, i);
        for (int j = std::max(0, i - t); j <= mpi; ++j)
            bezalfs(i, j) = bezalfs(ph - i, degree - j);
    }
    int knots_num;
    find_knots_num(knots_vector, knots_num);
    new_knots_vector.resize(knots_size + knots_num * t);
    new_control_points.resize(point_size, knots_size - degree - 1 + (knots_num - 1) * t);
    T ua = knots_vector[0];
    for (int i = 0; i <= ph; ++i)
        new_knots_vector[i] = ua;

    Eigen::Matrix<T, point_size, Eigen::Dynamic> bpts = control_points.block(0, 0, point_size, degree + 1);
    int a = degree;
    int b = degree + 1;
    int m = knots_size - 1;
    int r = -1;
    int cind = 1;
    int kind = ph + 1;
    new_control_points.col(0) = control_points.col(0);
    Eigen::VectorX<T> alphas(degree - 1);
    Eigen::Matrix<T, point_size, Eigen::Dynamic> ebpts(point_size, degree + t + 1);
    while (b < m)
    {
        int i = b;
        while (b < m && knots_vector[b] == knots_vector[b + 1])
            b += 1;
        int mul = b - i + 1;
        T ub = knots_vector[b];
        int older = r;
        r = degree - mul;
        int lbz = older > 0 ? (older + 2) / 2 : 1;
        int rbz = r > 0 ? ph - (r + 1) / 2 : ph;
        if (r > 0)
        {
            T number = ub - ua;
            for (int k = degree; k > mul; --k)
                alphas[k - mul - 1] = number / (knots_vector[a + k] - ua);
            for (int j = 1; j <= r; ++j)
            {
                int save = r - j;
                int s = mul + j;
                for (int k = degree; k >= s; --k)
                {
                    bpts.col(k) = alphas[k - s] * bpts.col(k) + (1.0 - alphas[k - s]) * bpts.col(k - 1);
                }
                next_bpts.col(save) = bpts.col(degree);
            }
        }

        Eigen::Vector<T, point_size> zero_vector;
        zero_vector.setConstant(0.0);
        for (int i = lbz; i <= ph; ++i)
        {
            ebpts.col(i) = zero_vector;
            int mpi = std::min(degree, i);
            for (int j = std::max(0, i - t); j <= mpi; ++j)
                ebpts.col(i) = ebpts.col(i) + bezalfs(i, j) * bpts.col(j);
        }
        if (older > 1)
        {
            int first = kind - 2;
            int last = kind;
            T den = ub - ua;
            T bet = (ub - knots_vector[kind - 1]) / den;
            for (int tr = 1; tr < older; ++tr)
            {
                int i = first, j = last, kj = j - kind + 1;
                while (j - i > tr)
                {
                    if (i < cind)
                    {
                        T alf = (ub - new_knots_vector[i]) / (ua - new_knots_vector[i]);
                        new_control_points.col(i) = alf * new_control_points.col(i) + (1.0 - alf) * new_control_points.col(i - 1);
                    }
                    if (j >= lbz)
                    {
                        if (j - tr <= kind - ph + older)
                        {
                            T gam = (ub - new_knots_vector[j - tr]) / den;
                            ebpts.col(kj) = gam * ebpts.col(kj) + (1.0 - gam) * ebpts.col(kj + 1);
                        }
                        else
                        {
                            ebpts.col(kj) = bet * ebpts.col(kj) + (1.0 - bet) * ebpts.col(kj + 1);
                        }
                    }
                    i += 1;
                    j -= 1;
                    kj -= 1;
                }
                first -= 1;
                last += 1;
            }
        }

        if (a != degree)
            for (int i = 0; i < ph - older; ++i)
            {
                new_knots_vector[kind] = ua;
                kind += 1;
            }
        for (int j = lbz; j <= rbz; ++j)
        {
            new_control_points.col(cind) = ebpts.col(j);
            cind += 1;
        }
        if (b < m)
        {
            bpts.block(0, 0, point_size, r) = next_bpts.block(0, 0, point_size, r);
            bpts.block(0, r, point_size, degree - r + 1) = control_points.block(0, b - degree + r, point_size, degree - r + 1);
            a = b; b += 1; ua = ub;
        }
        else
            for (int i = 0; i <= ph; ++i)
                new_knots_vector[kind + i] = ub;
    }
    return ENUM_NURBS::NURBS_SUCCESS;
}


template<typename T, int point_size>
ENUM_NURBS bez_degree_reduce(const Eigen::Matrix<T, point_size, Eigen::Dynamic> &control_points, Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points, T &max_error)
{
    int degree = control_points.cols() - 1;
    new_control_points.resize(point_size, degree);
    new_control_points.col(0) = control_points.col(0);
    new_control_points.col(degree - 1) = control_points.col(degree);
    int r = (degree - 1) / 2;

    for (int i = 1; i <= r; ++i)
    {
        T alph = static_cast<T>(i) / degree;
        new_control_points.col(i) = (control_points.col(i) - alph * new_control_points.col(i - 1)) / (1.0 - alph);
    }
    for (int i = degree - 2; i > r; --i)
    {
        T alph = static_cast<T>(i + 1) / degree;
        new_control_points.col(i) = (control_points.col(i + 1) - (1.0 - alph) * new_control_points.col(i + 1)) / alph;
    }
    if (degree % 2 == 0)
    {
        max_error = (control_points.col(r + 1) - (new_control_points.col(r) + new_control_points.col(r + 1)) / 2.0).norm();
    }
    else
    {
        T alph_r = static_cast<T> (r) / degree;
        T next_alph = static_cast<T> (r +1) / degree;
        Eigen::Vector<T, point_size> p_rl = (control_points.col(r) - alph_r * new_control_points.col(r - 1)) / (1.0 - alph_r);
        Eigen::Vector<T, point_size> p_rr = (control_points.col(r + 1) - (1.0 - next_alph) * new_control_points.col(r + 1)) / next_alph;
        new_control_points.col(r) = (p_rl + p_rr) / 2.0;
        max_error = (p_rl - p_rr).norm();
    }

    return ENUM_NURBS::NURBS_SUCCESS;
}

template<typename T, int point_size, bool is_rational>
ENUM_NURBS degree_reduce_curve(int degree, const Eigen::VectorX<T> &knots_vector, const Eigen::Matrix<T, point_size, Eigen::Dynamic> &control_points,
        Eigen::VectorX<T> &new_knots_points, Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points, T error = DEFAULT_ERROR)
{
    T TOL = -1;
    int cols = control_points.cols();
    if (is_rational == false)
    {
        TOL = error;
    }
    else
    {
        T w_min = control_points(3, 0);
        if (w_min <= 0)
            return ENUM_NURBS::NUBRS_WIEGHT_IS_NONPOSITIVE;
        T p_max = (control_points.block(0, 0, 3, 1) / control_points(3, 0)).norm();
        for (int col = 1; col < cols; ++col)
        {
            if (w_min > control_points(3, col))
                w_min = control_points(3, col);
            if (w_min <= 0)
                return ENUM_NURBS::NUBRS_WIEGHT_IS_NONPOSITIVE;
            T temp = (control_points.block(0, col, 3, 1) / control_points(3, col)).norm();
            if (p_max < temp)
                p_max = temp;
        }
        TOL = (error * w_min) / (1 + p_max);
    }

    int n = control_points.cols() - 1;
    int ph = degree - 1;// mh = ph;
    int kind = ph + 1;
    int r = -1;
    int a = degree;
    int b = degree + 1;
    int cind = 1;
    int mult = degree;
    int m = degree + n + 1;
    int knots_num;
    find_knots_num(knots_vector, knots_num);
    new_control_points.resize(point_size, cols - knots_num + 1);
    new_knots_points.resize(m + 1 - knots_num);
    new_control_points.col(0) = control_points.col(0);
    new_knots_points.block(0, 0, degree, 1) = knots_vector.block(0, 0, degree, 1);
    Eigen::Matrix<T, point_size, Eigen::Dynamic> bpts(point_size, degree + 1);
    bpts = control_points.block(0, 0, point_size, degree + 1);
    Eigen::VectorX<T> e(m);
    e.setConstant(0.0);
    Eigen::VectorX<T> alphas(degree - 1);
    Eigen::Matrix<T, point_size, Eigen::Dynamic> next_bpts(point_size, degree - 1);
    Eigen::Matrix<T, point_size, Eigen::Dynamic> rbpts(point_size, degree);
    while (b < m)
    {
        int i = b;
        while (b < m && knots_vector[b] == knots_vector[b + 1])
        {
            b += 1;
        }
        mult = b - i + 1;
        // mh += mult - 1;
        int oldr = r;
        r = degree - mult;
        int lbz = oldr > 0 ? (oldr + 2) / 2 : 1;
        if (r > 0)
        {
            T numer = knots_vector[b] - knots_vector[a];
            for (int k = degree; k > mult; --k) ///////？？？？？？？？？？？
                alphas[k - mult -1] = numer / (knots_vector[a + k] - knots_vector[a]);
            for (int j = 1; j <= r; ++j)
            {
                int save = r - j, s = mult + j;
                for (int k = degree; k >= s; --k)
                    bpts.col(k) = alphas[k-s] * bpts.col(k) + (1.0 - alphas[k - s]) * bpts.col(k - 1);
                next_bpts.col(save) = bpts.col(degree);
            }
        }
        T max_error;
        bez_degree_reduce<T, point_size>(bpts, rbpts, max_error);
        e[a] += max_error;
        if (e[a] > TOL)
            return ENUM_NURBS::NURBS_ERROR;

        if (oldr > 0)
        {
            int first = kind, last = kind;
            for (int k = 0; k < oldr; ++k)
            {
                i = first;
                int j = last;
                int kj = j - kind;
                while (j - i > k)
                {
                    T alfa = (knots_vector[a] - new_knots_points[i - 1]) / (knots_vector[b] - new_knots_points[i - 1]);
                    T beta = (knots_vector[a] - new_knots_points[j - k - 1]) / (knots_vector[b] - new_knots_points[j - k - 1]);
                    new_control_points.col(i - 1) = (new_control_points.col(i - 1) - (1.0 - alfa) * new_control_points.col(i - 2)) / alfa;
                    rbpts.col(kj) = (rbpts.col(kj) - beta * rbpts.col(kj + 1)) / (1.0 - beta);
                    i += 1;
                    j -= 1;
                    kj -= 1;
                }
                T Br = -1;
                if (j - i < k) // 偶数
                {
                    Br = (new_control_points.col(i - 2) - rbpts.col(kj + 1)).norm();
                }
                else    //奇数
                {
                    T delta = (knots_vector[a] - new_knots_points[i - 1]) / (knots_vector[b] - new_knots_points[i - 1]);
                    Eigen::Vector<T, point_size> A = delta * rbpts.col(kj + 1) + (1.0 - delta) * new_control_points.col(i - 2);
                    Br = (new_control_points.col(i - 1) - A).norm();
                }
                int K = a + oldr - k;
                int q = (2 * degree - k + 1) / 2;
                int L = K - q;
                for (int ii = L; ii <= a; ++ii)
                {
                    e[ii] += Br;
                    if (e[ii] > TOL)
                        return ENUM_NURBS::NURBS_ERROR;
                }
                first -= 1;
                last += 1;
            }
            cind = i - 1;
        }

        if (a != degree)
            for (int i = 0; i < ph - oldr; ++i)
            {
                new_knots_points[kind] = knots_vector[a];
                kind += 1;
            }
        for (int i = lbz; i <= ph; ++i)
        {
            new_control_points.col(cind) = rbpts.col(i);
            cind += 1;
        }

        if (b < m)
        {
            bpts.block(0, 0, point_size, r) = next_bpts.block(0, 0, point_size, r);
            bpts.block(0, r, point_size, degree - r + 1) = control_points.block(0, b - degree + r, point_size, degree - r + 1);
            a = b;
            b += 1;
        }
        else
        {
            new_knots_points.block(kind, 0, ph + 1, 1).setConstant(knots_vector[b]);
        }

    }

    return ENUM_NURBS::NURBS_SUCCESS;

}



// is_rational == false
template<typename T, int dim, bool is_rational>
struct prallel_projection_curve
{
    static ENUM_NURBS parallel_projection(const Eigen::Vector<T, dim> &reference_point, const Eigen::Vector<T, dim> &normal, const Eigen::Vector<T, dim> &projection_direction, 
        const Eigen::Matrix<T, dim, Eigen::Dynamic> &control_points, Eigen::Matrix<T, dim, Eigen::Dynamic> &new_control_points)
    {
        int cols = control_points.cols();
        T np = normal.dot(projection_direction);
        if (np < DEFAULT_ERROR)
            return ENUM_NURBS::NURBS_ERROR;
        T nr = normal.dot(reference_point);
        T np_nr = nr / np;

        new_control_points.resize(dim, cols);
        for (int index = 0; index < cols; ++index)
        {
            Eigen::Vector<T, dim> point = control_points.col(index);
            new_control_points.col(index) = point + (np_nr - normal.dot(point) / np) * projection_direction;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }
};

template<typename T, int dim>
struct prallel_projection_curve<T, dim, true>
{
    static ENUM_NURBS parallel_projection(const Eigen::Vector<T, dim> &reference_point, const Eigen::Vector<T, dim> &normal, const Eigen::Vector<T, dim> &projection_direction, 
        const Eigen::Matrix<T, dim + 1, Eigen::Dynamic> &control_points, Eigen::Matrix<T, dim + 1, Eigen::Dynamic> &new_control_points)
    {
        int cols = control_points.cols();
        T np = normal.dot(projection_direction);
        if (np < DEFAULT_ERROR)
            return ENUM_NURBS::NURBS_ERROR;
        T nr = normal.dot(reference_point);
        T np_nr = nr / np;
        new_control_points.resize(dim + 1, cols);
        for (int index = 0; index < cols; ++index)
        {
            T weight = control_points(dim, index);
            Eigen::Vector<T, dim> pi = control_points.block(0, index, dim, 1) / weight;
            new_control_points.block(0, index, dim, 1) = (pi - (np_nr + normal.dot(pi) / np) * projection_direction) * weight;
            new_control_points(dim, index) = weight;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }
};

// is_rational == false
template<typename T, int dim, bool is_rational>
struct perspective_projection_curve
{
    static ENUM_NURBS perspective_projection(const Eigen::Vector<T, dim> &reference_point, const Eigen::Vector<T, dim> &normal, const Eigen::Vector<T, dim> &eye, 
        const Eigen::Matrix<T, dim, Eigen::Dynamic> &control_points, Eigen::Matrix<T, dim + 1, Eigen::Dynamic> &new_control_points)
    {
        int cols = control_points.cols();
        T ne = normal.dot(eye);
        T nr = normal.dot(reference_point);
        T ne_nr = ne - nr;

        new_control_points.resize(dim + 1, cols);
        for (int index = 0; index < cols; ++index)
        {
            Eigen::Vector<T, dim> point = control_points.col(index);
            T np = normal.dot(point);
            T weight = ne - np;
            if (weight < DEFAULT_ERROR)
                return ENUM_NURBS::NURBS_ERROR;
            new_control_points.block(0, index, dim, 1) = ((ne_nr / weight) * point + ((nr - normal.dot(point)) / weight) * eye) * weight;
            new_control_points(dim, index) = weight;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }
};

template<typename T, int dim>
struct perspective_projection_curve<T, dim, true>
{
    static ENUM_NURBS perspective_projection(const Eigen::Vector<T, dim> &reference_point, const Eigen::Vector<T, dim> &normal, const Eigen::Vector<T, dim> &eye, 
        const Eigen::Matrix<T, dim + 1, Eigen::Dynamic> &control_points, Eigen::Matrix<T, dim + 1, Eigen::Dynamic> &new_control_points)
    {
        int cols = control_points.cols();
        T ne = normal.dot(eye);
        T nr = normal.dot(reference_point);
        T ne_nr = ne - nr;

        new_control_points.resize(dim + 1, cols);
        for (int index = 0; index < cols; ++index)
        {
            T old_weight = control_points(dim, index);
            Eigen::Vector<T, dim> point = control_points.block(0, index, dim, 1) / old_weight;
            T np = normal.dot(point);
            T denominator = ne - np;
            T weight = (ne - np) * old_weight;
            if (denominator < DEFAULT_ERROR)
                return ENUM_NURBS::NURBS_ERROR;
            new_control_points.block(0, index, dim, 1) = ((ne_nr / denominator) * point + ((nr - normal.dot(point)) / denominator) * eye) * weight;
            new_control_points(dim, index) = weight;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }
};


std::vector<std::vector<int>> find_sum_n(int k, int n, int add_num, bool flag)
{
    std::vector<std::vector<int>> reuslt;
    if (k == 1)
    {
        std::vector<int> temp{n};
        if (flag)
            temp.push_back(add_num);
        reuslt.push_back(temp);
        return reuslt;
    }
    if (n == 0)
    {
        std::vector<int> temp(k, 0);
        if (flag)
            temp.push_back(add_num);
        reuslt.push_back(temp);
        return reuslt;
    }
    for (int last = 0; last <= n; ++last)
    {
        std::vector<std::vector<int>> vec = find_sum_n(k - 1, n -last, last, true);
        if (flag)
        {
            for (auto &v : vec)
            {
                v.push_back(add_num);
            }
        }

        reuslt.insert(reuslt.end(), vec.begin(), vec.end());
    }
    return reuslt;
}



template<typename T, int point_size>
ENUM_NURBS reparameter_bezier_curve(int old_degree, const Eigen::VectorX<T> &old_knots_vector, const Eigen::Matrix<T, point_size, Eigen::Dynamic> &old_control_points,
    int function_degree, const Eigen::VectorX<Eigen::Vector<T, 1>> &function_low_ders, const Eigen::VectorX<Eigen::Vector<T, 1>> &function_high_ders, T delta_s,
    Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points)
{
    int new_degree = function_degree * old_degree;
    int pqf =  std::tgamma<int>(new_degree + 1);
    int ml = function_low_ders.size() - 1;
    int mr = function_high_ders.size() - 1;
    
    Eigen::MatrixX<int> Bin = binary_coeff(ml + 1);

    Eigen::VectorX<Eigen::Vector<T, point_size>> ders_low(ml + 1);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ders_temp;
    int new_n = std::min(ml, old_degree);
    ders_basis_funs<T>(old_degree, new_n, old_degree, old_knots_vector[0], old_knots_vector, ders_temp);
    Eigen::Vector<T, point_size> zero;
    zero.setConstant(0.0);
    for (int idx = 0; idx <= new_n; ++idx)
    {
        Eigen::Vector<T, point_size> vec = old_control_points.block(0, 0, point_size, old_degree + 1) * ders_temp.col(idx);
        ders_low[idx] = vec;
    }
    for (int idx = new_n + 1; idx <= ml; ++idx)
    {
        ders_low[idx] = zero;
    }

    new_n = std::min(mr, old_degree);
    ders_basis_funs<T>(old_degree, new_n, old_degree, old_knots_vector[old_degree + 1], old_knots_vector, ders_temp);
    Eigen::VectorX<Eigen::Vector<T, point_size>> ders_high(mr + 1);
    int knots_index = old_knots_vector.size() - old_degree * 2 - 2;
    for (int idx = 0; idx <= new_n; ++idx)
    {
        Eigen::Vector<T, point_size> vec = old_control_points.block(0, knots_index, point_size, old_degree + 1) * ders_temp.col(idx);
        ders_high[idx] = vec;
    }
    for (int idx = new_n + 1; idx <= mr; ++idx)
    {
        ders_high[idx] = zero;
    }

    // mr <= ml

    Eigen::VectorX<Eigen::Vector<T, point_size>> ders_s_low(ml + 1);
    for (int index = 0; index <= ml; ++index)
    {
        ders_s_low[index].setConstant(0.0);
    }
    ders_s_low[0] = ders_low[0];

    Eigen::VectorX<Eigen::Vector<T, point_size>> ders_s_high(mr + 1);
    ders_s_high[0] = ders_high[0];
    for (int index = 1; index <= mr; ++index)
        ders_s_high[index] = zero;

    for (int n = 1; n <= mr; ++n)
    {
        for (int j = 0; j <= n; ++j)
        {
            std::vector<std::vector<int>> temp = find_sum_n(n, j, 0, false);
            for (auto vec : temp)
            {
                int vec_sum = 0;
                for (int index = 0; index < n; ++index)
                    vec_sum += (index + 1) * vec[index];
                if (vec_sum == n)
                {
                    T coeff1 = std::tgamma<int>(n + 1) * std::pow(function_low_ders[1][0], vec[0]);
                    T coeff2 = std::tgamma<int>(n + 1) * std::pow(function_high_ders[1][0], vec[0]);
                    int denominator = std::tgamma<int>(vec[0] + 1);
                    for (int index = 1; index < n; ++index)
                    {
                        denominator *= std::tgamma<int>(vec[index] + 1) * std::pow(std::tgamma<int>(index + 2), vec[index]);
                        coeff1 *= std::pow(function_low_ders[index + 1][0], vec[index]);
                        coeff2 *= std::pow(function_high_ders[index + 1][0], vec[index]);
                    }
                    ders_s_low[n] += ders_low[j] * (coeff1 / denominator);
                    ders_s_high[n] += ders_high[j] * (coeff2 / denominator);
                }
            }   
        }
    }

    if (new_degree % 2 == 0)
    {
        for (int j = 0; j <= ml; ++j)
        {
            std::vector<std::vector<int>> temp = find_sum_n(ml, j, 0, false);
            for (auto vec : temp)
            {
                int vec_sum = 0;
                for (int index = 0; index < ml; ++index)
                    vec_sum += (index + 1) * vec[index];
                if (vec_sum == ml)
                {
                    T coeff1 = std::tgamma<int>(ml + 1) * std::pow(function_low_ders[1][0], vec[0]);
                    int denominator = std::tgamma<int>(vec[0] + 1);
                    for (int index = 1; index < ml; ++index)
                    {
                        denominator *= std::tgamma<int>(vec[index] + 1) * std::pow(std::tgamma<int>(index + 2), vec[index]);
                        coeff1 *= std::pow(function_low_ders[index + 1][0], vec[index]);
                    }
                    ders_s_low[ml] += ders_low[j] * (coeff1 / denominator);
                }
            }   
        }
    }


    for (int index = 1; index <= mr; ++index)
    {
        int num = std::tgamma<int>(new_degree + 1 - index);
        T s = std::pow(delta_s, index);
        int coeff = index % 2 == 0 ? 1 : -1;
        new_control_points.col(index) = ((num * s) / pqf) * ders_s_low[index];
        new_control_points.col(new_degree - index) = ((coeff * num * s) / pqf) * ders_s_high[index];
        for (int j = 0; j < index; ++j)
        {
            int coeff1 = (index + j - 1) % 2 == 0 ? 1 : -1;
            new_control_points.col(index) = new_control_points.col(index) + coeff1 * Bin(index, j) * new_control_points.col(j);
            new_control_points.col(new_degree - index) = new_control_points.col(new_degree - index) + coeff1 * Bin(index, j) * new_control_points.col(new_degree - j);
        }
    }

    if (new_degree % 2 == 0)
    {
        int num = std::tgamma<int>(new_degree + 1 - ml);
        T s = std::pow(delta_s, ml);
        new_control_points.col(ml) = ((num * s) / pqf) * ders_s_low[ml];
        for (int j = 0; j < ml; ++j)
        {
            int coeff1 = (ml + j - 1) % 2 == 0 ? 1 : -1;
            new_control_points.col(ml) = new_control_points.col(ml) + coeff1 * Bin(ml, j) * new_control_points.col(j);
        }
    }

    return ENUM_NURBS::NURBS_SUCCESS; 
}

//误差可以为0
template<typename T, int point_size>
ENUM_NURBS merge_two_curve(int left_degree, int right_degree, const Eigen::VectorX<T> &left_knots, 
    const Eigen::VectorX<T> &right_knots, const Eigen::Matrix<T, point_size, Eigen::Dynamic> &left_control_points,
    const Eigen::Matrix<T, point_size, Eigen::Dynamic> &right_control_points, Eigen::VectorX<T> &new_knots_vector,
    Eigen::Matrix<T, point_size, Eigen::Dynamic> &new_control_points, T eps = DEFAULT_ERROR)
{
    if (left_degree != right_degree)
        return ENUM_NURBS::NURBS_ERROR;
    int left_knots_count = left_knots.size();
    if (std::abs(left_knots[left_knots_count - 1] - right_knots[0]) > eps)
        return ENUM_NURBS::NURBS_ERROR;
    Eigen::Vector<T, point_size> vec = left_control_points.col(left_knots_count - left_degree - 2) - right_control_points.col(0);
    if (vec.squaredNorm() > eps * eps)
        return ENUM_NURBS::NURBS_ERROR;
    
    int right_knots_count = right_knots.size();
    new_knots_vector.resize(left_knots_count + right_knots_count - left_degree - 2);
    for (int index = 0; index < left_knots_count - 1; ++index)
    {
        new_knots_vector[index] = left_knots[index];
    }

    for (int index = right_degree + 1; index < right_knots_count; ++index)
    {
        new_knots_vector[index + left_knots_count - 2 - right_degree] = right_knots[index];
    }

    int left_control_points_count = left_knots_count - left_degree - 1;
    int right_control_points_count = right_knots_count - right_degree - 1;
    new_control_points.resize(point_size, right_control_points_count + left_control_points_count - 1);
    new_control_points.block(0, 0, point_size, left_control_points_count) = left_control_points;
    new_control_points.block(0, left_control_points_count - 1, point_size, right_control_points_count) = right_control_points;
    return ENUM_NURBS::NURBS_SUCCESS;

}
