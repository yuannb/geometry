#include "declare.h"
#include "nurbs_tool.h"
#include <array>


namespace tnurbs
{
    //两区间相加
    template<typename T>
    Interval<T> plus(const Interval<T> &left, const Interval<T> &right, T eps = PRECISION<T>::value)
    {
        Interval<T> result;
        result.m_interval[0] = left.m_interval[0] + right.m_interval[0] - eps;
        result.m_interval[1] = left.m_interval[1] + right.m_interval[1] + eps;
        return result;
    }

    template<typename T, unsigned dim>
    std::array<Interval<T>, (unsigned)dim> plus(const std::array<Interval<T>, (unsigned)dim>& left, const std::array<Interval<T>, (unsigned)dim>& right, T eps = PRECISION<T>::value)
    {
        std::array<Interval<T>, (unsigned)dim> result;
        for (unsigned index = 0; index < dim; ++index)
            result[index] = plus(left[index], right[index]);
        return result;
    }


    //两区间相减
    template<typename T>
    Interval<T> minus(const Interval<T> &left, const Interval<T> &right, T eps = PRECISION<T>::value)
    {
        Interval<T> result;
        result.m_interval[0] = left.m_interval[0] - right.m_interval[1] - eps;
        result.m_interval[1] = left.m_interval[1] - right.m_interval[0] + eps;
        return result;
    }

    template<typename T, unsigned dim>
    std::array<Interval<T>, (unsigned)dim> minus(const std::array<Interval<T>, (unsigned)dim>& left, const std::array<Interval<T>, (unsigned)dim>& right, T eps = PRECISION<T>::value)
    {
        std::array<Interval<T>, (unsigned)dim> result;
        for (unsigned index = 0; index < dim; ++index)
            result[index] = minus(left[index], right[index]);
        return result;
    }

    //两区间相乘
    template<typename T>
    Interval<T> mutilply(const Interval<T> &left, const Interval<T> &right, T eps = PRECISION<T>::value)
    {
        Interval<T> result;
        T x1 = left.m_interval[0] * right.m_interval[0];
        T x2 = left.m_interval[0] * right.m_interval[1];
        T x3 = left.m_interval[1] * right.m_interval[0];
        T x4 = left.m_interval[1] * right.m_interval[1];

        result.m_interval[0] = std::min({x1, x2, x3, x4}) - eps;
        result.m_interval[1] = std::max({x1, x2, x3, x4}) + eps;
        return result;
    }

    template<typename T, unsigned dim>
    std::array<Interval<T>, (unsigned)dim> mutilply(const std::array<Interval<T>, (unsigned)dim>& left, const Interval<T>& right, T eps = PRECISION<T>::value)
    {
        std::array<Interval<T>, (unsigned)dim> result;
        for (unsigned index = 0; index < dim; ++index)
            result[index] = mutilply(left[index], right);
        return result;
    }

    //两区间相除
    template<typename T>
    bool divide(const Interval<T> &left, const Interval<T> &right, Interval<T> &result, T eps = PRECISION<T>::value)
    {
        if (0 >= right.m_interval[0] && 0 <= right.m_interval[1])
            return false;
        T x1 = left.m_interval[0] / right.m_interval[0];
        T x2 = left.m_interval[0] / right.m_interval[1];
        T x3 = left.m_interval[1] / right.m_interval[0];
        T x4 = left.m_interval[1] / right.m_interval[1];

        result.m_interval[0] = std::min({x1, x2, x3, x4}) - eps;
        result.m_interval[1] = std::max({x1, x2, x3, x4}) + eps;
        return true;
    }


    template<typename T, unsigned dim>
    bool divide(const std::array<Interval<T>, (unsigned)dim> &left, const Interval<T> &right, std::array<Interval<T>, (unsigned)dim> &result, T eps = PRECISION<T>::value)
    {
        for (unsigned index = 0; index < dim; ++index)
            if (divide(left[index], right, result[index]) == false)
                return false;
        return true;
    }



    //两个向量作点积
    template<typename T, int dim = 3>
    Interval<T> dot_product(const std::array<Interval<T>, (unsigned)dim> &left, const std::array<Interval<T>, (unsigned)dim> &right, T eps = PRECISION<T>::value)
    {
        Interval<T> result;
        for (int index = 0; index < dim; ++index)
        {
            result = plus(result, mutilply(left[index], right[index]));
        }
        return result;
    }

    //两个三维向量作叉乘
    template<typename T>
    std::array<Interval<T>, 3> vector_product(const std::array<Interval<T>, 3> &left, const std::array<Interval<T>, 3> &right, T eps = PRECISION<T>::value)
    {
        std::array<Interval<T>, 3> result;
        result[0] = minus(mutilply(left[1], right[2]), mutilply(left[2], right[1]));
        result[1] = minus(mutilply(left[2], right[0]), mutilply(left[0], right[2]));
        result[2] = minus(mutilply(left[0], right[1]), mutilply(left[1], right[0]));  
        
        return result;
    }

}
