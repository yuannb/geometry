#pragma once
#include "declare.h"
#include "nurbs_tool.h"
#include <array>


namespace tnurbs
{
    template<typename T>
    Interval<T> plus(const Interval<T>& left, const Interval<T>& right, T eps = PRECISION<T>::value)
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


    template<typename T>
    Interval<T> minus(const Interval<T>& left, const Interval<T>& right, T eps = PRECISION<T>::value)
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

    template<typename T>
    Interval<T> mutilply(const Interval<T>& left, const Interval<T>& right, T eps = PRECISION<T>::value)
    {
        Interval<T> result;
        T x1 = left.m_interval[0] * right.m_interval[0];
        T x2 = left.m_interval[0] * right.m_interval[1];
        T x3 = left.m_interval[1] * right.m_interval[0];
        T x4 = left.m_interval[1] * right.m_interval[1];

        result.m_interval[0] = std::min<T>({ x1, x2, x3, x4 }) - eps;
        result.m_interval[1] = std::max<T>({ x1, x2, x3, x4 }) + eps;
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

    template<typename T>
    bool divide(const Interval<T>& left, const Interval<T>& right, Interval<T>& result, T eps = PRECISION<T>::value)
    {
        if (0 >= right.m_interval[0] && 0 <= right.m_interval[1])
            return false;
        T x1 = left.m_interval[0] / right.m_interval[0];
        T x2 = left.m_interval[0] / right.m_interval[1];
        T x3 = left.m_interval[1] / right.m_interval[0];
        T x4 = left.m_interval[1] / right.m_interval[1];

        result.m_interval[0] = std::min<T>({ x1, x2, x3, x4 }) - eps;
        result.m_interval[1] = std::max<T>({ x1, x2, x3, x4 }) + eps;
        return true;
    }


    template<typename T, unsigned dim>
    bool divide(const std::array<Interval<T>, (unsigned)dim>& left, const Interval<T>& right, std::array<Interval<T>, (unsigned)dim>& result, T eps = PRECISION<T>::value)
    {
        for (unsigned index = 0; index < dim; ++index)
            if (divide(left[index], right, result[index]) == false)
                return false;
        return true;
    }



    template<typename T, int dim = 3>
    Interval<T> dot_product(const std::array<Interval<T>, (unsigned)dim>& left, const std::array<Interval<T>, (unsigned)dim>& right, T eps = PRECISION<T>::value)
    {
        Interval<T> result;
        for (int index = 0; index < dim; ++index)
        {
            result = plus(result, mutilply(left[index], right[index]));
        }
        return result;
    }

    template<typename T>
    std::array<Interval<T>, 3> vector_product(const std::array<Interval<T>, 3>& left, const std::array<Interval<T>, 3>& right, T eps = PRECISION<T>::value)
    {
        std::array<Interval<T>, 3> result;
        result[0] = minus(mutilply(left[1], right[2]), mutilply(left[2], right[1]));
        result[1] = minus(mutilply(left[2], right[0]), mutilply(left[0], right[2]));
        result[2] = minus(mutilply(left[0], right[1]), mutilply(left[1], right[0]));

        return result;
    }


    namespace interval_algorithm
    {

        template<typename T>
        Box<T, 3> cross(const Box<T, 3>& left, const Box<T, 3>& right)
        {
            Box<T, 3> result;
            const std::array<Box<T, 1>, 3>& lhs = left.sperate_box();
            const std::array<Box<T, 1>, 3>& rhs = right.sperate_box();
            result.set_index_interval(0, lhs[1] * rhs[2] - lhs[2] * rhs[1]);
            result.set_index_interval(1, lhs[2] * rhs[0] - lhs[0] * rhs[2]);
            result.set_index_interval(2, lhs[0] * rhs[1] - lhs[1] * rhs[0]);
            return result;
        }

        template<typename T, int dim>
        Box<T,1> dot(const Box<T, dim>& left, const Box<T, dim>& right)
        {
            Box<T, 1> result;
            result.Min[0] = 0;
            result.Max[0] = 0;
            const std::array<Box<T, 1>, dim>& lhs = left.sperate_box();
            const std::array<Box<T, 1>, dim>& rhs = right.sperate_box();
            
            for (int index = 0; index < dim; ++index)
            {
                result = result + (lhs[index] * rhs[index]);
            }
            return result;
        }

        template<typename T, int dim>
        ENUM_NURBS normalized(const Box<T, dim>& vec, Box<T, dim>& unit_vec)
        {
            T min2 = 0;
            T max2 = 0;
            for (int index = 0; index < dim; ++index)
            {
                T max = std::max<int>(std::abs(vec.Min[index]), std::abs(vec.Max[index]));
                T min = vec.Min[index] * vec.Max[index] > 0 ? std::min<int>(std::abs(vec.Min[index]), std::abs(vec.Max[index])) : 0;
                min2 += min * min;
                max2 += max * max;
            }
            if (min2 < KNOTS_EPS<T>::value * KNOTS_EPS<T>::value)
            {
                return ENUM_NURBS::NURBS_ERROR;
            }
            for (int index = 0; index < dim; ++index)
            {
                T max = std::max<int>(std::abs(vec.Min[index]), std::abs(vec.Max[index]));
                T min = vec.Min[index] * vec.Max[index] > 0 ? std::min<int>(std::abs(vec.Min[index]), std::abs(vec.Max[index])) : 0;

                T max2_t = max2 - max * max;
                T min2_t = min2 - min * min;
                min2_t += vec.Max[index] * vec.Max[index];
                max2_t += vec.Min[index] * vec.Min[index];
                unit_vec.Min[index] = vec.Min[index] / std::sqrt(max2_t);
                unit_vec.Max[index] = vec.Max[index] / std::sqrt(min2_t);
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        template<typename T, int dim1, int dim2>
        ENUM_NURBS divide(const Box<T, dim1>& left, const Box<T, dim2>& right, Box<T, dim1>& result)
        {
            static_assert(dim2 == 1 || dim1 == dim2, "dim2 != 1 && dim1 != dim2");
            if constexpr (dim1 == dim2)
            {
                for (int index = 0; index < dim1; ++index)
                {
                    if (KNOTS_EPS<T>::value > right.Min[index] && -KNOTS_EPS<T>::value < right.Max[index])
                    {
                        return ENUM_NURBS::NURBS_ERROR;
                    }
                    T x1 = left.Min[index] / right.Min[index];
                    T x2 = left.Min[index] / right.Max[index];
                    T x3 = left.Max[index] / right.Min[index];
                    T x4 = left.Max[index] / right.Max[index];
                    result.Min[index] = std::min<T>({ x1, x2, x3, x4 });
                    result.Max[index] = std::max<T>({ x1, x2, x3, x4 });
                }
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            else
            {
                if (KNOTS_EPS<T>::value > right.Min[0] && -KNOTS_EPS<T>::value < right.Max[0])
                {
                    return ENUM_NURBS::NURBS_ERROR;
                }
                for (int index = 0; index < dim1; ++index)
                {
                    T x1 = left.Min[index] / right.Min[0];
                    T x2 = left.Min[index] / right.Max[0];
                    T x3 = left.Max[index] / right.Min[0];
                    T x4 = left.Max[index] / right.Max[0];
                    result.Min[index] = std::min<T>({ x1, x2, x3, x4 });
                    result.Max[index] = std::max<T>({ x1, x2, x3, x4 });
                }
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            
        }

    };

};


