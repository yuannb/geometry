#pragma once
#include "array"
namespace tnurbs
{
    template<typename T> struct geo_traits;

    template<typename T> struct geo_traits<const T> : geo_traits<T> {};

    template<typename T>
    struct MAX_WEIGHT
    {
    };

    template<>
    struct MAX_WEIGHT<double>
    {
        constexpr static double value = 1e5;
    };

    template<>
    struct MAX_WEIGHT<float>
    {
        constexpr static float value = 1e3;
    };

    template<typename T>
    struct MIN_WEIGHT
    {
    };

    template<>
    struct MIN_WEIGHT<double>
    {
        constexpr static double value = 1e-4;
    };

    template<>
    struct MIN_WEIGHT<float>
    {
        constexpr static float value = 1e-2;
    };


template<typename T>
    struct KNOTS_EPS
    {
    };

    template<>
    struct KNOTS_EPS<double>
    {
        constexpr static double value = 1e-8;
    };

    template<>
    struct KNOTS_EPS<float>
    {
        constexpr static float value = 1e-4;
    };


    template<typename T>
    struct TDEFAULT_ERROR
    {

    };

    template<>
    struct TDEFAULT_ERROR<double>
    {
        constexpr static double value = 1e-4;
    };

    template<>
    struct TDEFAULT_ERROR<float>
    {
        constexpr static double value = 1e-2;
    };

    template<typename T, int dim>
    struct frame
    {
        Eigen::Matrix<T, dim, dim> basis;
        Eigen::Vector<T, dim> origin;
        frame() = default;
        ~frame () { }
    };
    
    enum ENUM_DIRECTION
    {
        U_DIRECTION = 0,
        V_DIRECTION = 1
    };

    enum ENUM_LIMITDIRECTION
    {
        LEFT = 0,
        RIGHT = 1
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
        NURBS_POINT_IS_NOT_ON_CURVE = 7,
        NURBS_CHORD_IS_ZERO = 8,
        NURBS_PARAM_IS_INVALID = 9
    };


    template<typename T, int dim>
    struct Box
    {
        Eigen::Vector<T, dim> Min;
        Eigen::Vector<T, dim> Max;
        Box() = default;

        Box& operator=(const Box<T, dim> &other)
        {
            if (this != &other)
            {
                Min = other.Min;
                Max = other.Max;
            }
            return *this;
        }

        bool is_contain_point(const Eigen::Vector<T, dim> &point)
        {
            for (int index = 0; index < dim; ++index)
            {
                if (point[index] > Max[index] || point[index] < Min[index])
                    return false;
            }
            return true;
        }

        T eval_minimal_distance(const Eigen::Vector<T, dim> &point) const
        {
            struct Index
            {
                Eigen::Vector<int, dim> index;
                Index() { index.setConstant(0); }
                ENUM_NURBS add_one()
                {
                    index[dim - 1] += 1;
                    for (int i = dim - 1; i > 0; --i)
                    {
                        if (index[i] != 2)
                            break;
                        else
                        {
                            index[i] -= 1;
                            index[i - 1] += 1;
                        }
                    }
                    if (index[0] == 2)
                        return ENUM_NURBS::NURBS_ERROR;
                    return ENUM_NURBS::NURBS_SUCCESS;
                }
                ENUM_NURBS get_point(const Box<T, dim> &box, Eigen::Vector<T, dim> &point) const
                {
                    for (int i = 0; i < dim; ++i)
                    {
                        if (index[i] == 0)
                            point[i] = box.Min[i];
                        else
                            point[i] = box.Max[i];
                    }             
                    return ENUM_NURBS::NURBS_SUCCESS;
                }
            };
            Index idx;
            T dis = -1.0;
            constexpr int point_count = std::pow(2, dim);
            for (int i = 0; i < point_count; ++i)
            {
                Eigen::Vector<T, dim> vec;
                idx.get_point(*this, vec);
                T current_dis = (vec - point).norm();
                if (current_dis < dis || current_dis < 0.0)
                    dis = current_dis;
            }
            //计算point到box的平面的距离
            std::vector<int> out_index;
            for (int i = 0; i < dim; ++i)
            {
                if (point[i] > Max[i] || point[i] < Min[i])
                    out_index.push_back(i);
            }
            int out_index_count = out_index.size();
            if (out_index_count == 0)
            {
                for (int i = 0; i < dim; ++i)
                {
                    T d1 = point[i] - Min[i];
                    T d2 = Max[i] - point[i];
                    dis = std::min({d1, d2, dis});
                }
            }
            else if (out_index_count == 1)
            {
                T d1 = std::abs(point[out_index[0]] - Min[out_index[0]]);
                T d2 = std::abs(point[out_index[0]] - Max[out_index[0]]);
                dis = std::min({d1, d2, dis}); 
            }
            return dis;
        }
    };



}

