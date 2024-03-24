#pragma once
#include "array"
#include <Eigen/Dense>
#include <Eigen/Core>
namespace tnurbs
{
    constexpr int MAXINTERSETORPOINTNUMBER = 100000;

    // enum ENGEOMETRYTYPE
    // {

    // }

    enum ENGEOMETRYTYPE
    {
        UNKOWN = -1,
        //无理nurbs曲线
        NUNBS_CURVE = 0,
        //有理nurbs曲线
        NURBS_CURVE = 1
    };



    template<typename T> struct geo_traits;

    template<typename T> struct geo_traits<const T> : geo_traits<T> {};

    template<typename T>
    struct PRECISION
    {
    };

    template<>
    struct PRECISION<double>
    {
        constexpr static double value = 1e-10;
    };

    template<>
    struct PRECISION<float>
    {
        constexpr static float value = 1e-5;
    };




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
        Box(const Box &box): Min(box.Min), Max(box.Max) { }
        // Box(T min, T max) : Min(min), Max(max) { };
        Box(const Eigen::Vector<T, dim> &min, const Eigen::Vector<T, dim> &max): Min(min), Max(max) { }
        Box& operator=(const Box<T, dim> &other)
        {
            if (this != &other)
            {
                Min = other.Min;
                Max = other.Max;
            }
            return *this;
        }

        bool is_contain_point(const Eigen::Vector<T, dim> &point) const
        {
            for (int index = 0; index < dim; ++index)
            {
                if (point[index] > Max[index] || point[index] < Min[index])
                    return false;
            }
            return true;
        }

        bool is_contain_box(const Box &bx) const
        {
            for (int index = 0; index < dim; ++index)
            {
                if (bx.Min[index] < Min[index] || bx.Max[index] > Max[index])
                    return false;
            }
            return true;
        }

        T eval_minimal_distance(const Eigen::Vector<T, dim> &point) const
        {        
            T dis = 0.0;
            for (int i = 0; i < dim; ++i)
            {
                if (point[i] > Max[i] || point[i] < Min[i])
                {
                    T min_dis = std::min(std::abs(point[i] - Max[i]), std::abs(point[i] - Min[i]));
                    dis += (min_dis * min_dis);
                }
            }
            return std::sqrt(dis);
        }

        T eval_maximal_distance(const Eigen::Vector<T, dim> &point) const
        {        
            T dis = 0.0;
            for (int i = 0; i < dim; ++i)
            {
                T max_dis = std::max(std::abs(point[i] - Max[i]), std::abs(point[i] - Min[i]));
                dis += max_dis * max_dis;
            }
            return std::sqrt(dis);
        }
    
        std::array<Box<T, dim>, 2> split_at_middle(int dimension) const
        {
            Box<T, dim> box1(*this), box2(*this);
            box1.Max[dimension] = Min[dimension] + (Max[dimension] - Min[dimension]) / 2.0;
            box2.Min[dimension] = box1.Max[dimension];
            return { box1, box2 };
        }

        bool intersect(const Box &box, Box &intersect_box)
        {
            for (int index = 0; index < dim; ++index)
            {
                intersect_box.Min[index] = std::max(Min[index], box.Min[index]);
                intersect_box.Max[index] = std::min(Max[index], box.Max[index]);
                if (intersect_box.Max[index] < intersect_box.Min[index])
                {
                    return false;
                }
            }
            return true;
        }
    
        Box &scale(T s)
        {
            Min *= s;
            Max *= s;
            return *this;
        }

        //谨慎使用
        Box<T, dim> plus(const Box<T, dim> &box, T eps = 0.0) const
        {
            Box<T, dim> new_box;
            Eigen::Vector<T, dim> eps_vector;
            eps_vector.setConstant(eps);
            new_box.Min = box.Min + Min - eps_vector;
            new_box.Max = box.Max + Max + eps_vector;
            return new_box;
        }       

        template<int dim2>
        Box<T, dim + dim2> product_box(const Box<T, dim2> &left_box) const
        {
            Eigen::Vector<T, dim + dim2> new_min;
            new_min.template block<dim, 1>(0, 0) = Min;
            new_min.template block<dim2, 1>(dim, 0) = left_box.Min;
            
            Eigen::Vector<T, dim + dim2> new_max;
            new_max.template block<dim, 1>(0, 0) = Max;
            new_max.template block<dim2, 1>(dim, 0) = left_box.Max;
            return Box<T, dim + dim2>(new_min, new_max);
        }

        Eigen::Vector<T, dim> get_middle_point() const
        {
            Eigen::Vector<T, dim> middle_point = (Min + Max) / 2.0;
            return middle_point;
        }

    };



}

