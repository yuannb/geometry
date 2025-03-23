#pragma once
#include "array"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <math.h>
#include <algorithm>
// using M_PI
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
    struct TDEFAULTANGLETOL
    {
    };

    template<>
    struct TDEFAULTANGLETOL<double>
    {
        constexpr static double value = 0.1;
    };

    template<>
    struct TDEFAULTANGLETOL<float>
    {
        constexpr static float value = 0.3;
    };


    template<typename T>
    struct SPACE_PARAM_TIMES
    {
    };

    template<>
    struct SPACE_PARAM_TIMES<double>
    {
        constexpr static double value = 1e-4;
    };

    template<>
    struct SPACE_PARAM_TIMES<float>
    {
        constexpr static float value = 1e-3;
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
    struct TINFINITE
    {

    };

    template<>
    struct TINFINITE<double>
    {
        constexpr static double value = 1e4;
    };

    template<>
    struct TINFINITE<float>
    {
        constexpr static double value = 1e3;
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
        NURBS_PARAM_IS_INVALID = 9,
        NURBS_ISOLATED_TANGENTIAL_POINT = 10,
        NURBS_HIGH_ORDER_TANGENTIAL = 11
    };


    template<typename T>
    struct Interval
    {
        Eigen::Vector2<T> m_interval;
        Interval() : m_interval{ 0, 0 } { };
        Interval(T low, T high) { m_interval[0] = low; m_interval[1] = high; }
        Interval(const Interval& interval) { m_interval = interval.m_interval; }
        T get_low() const { return m_interval[0]; }
        T get_high() const { return m_interval[1]; }
        Interval& operator=(const Interval& rhs)
        {
            m_interval = rhs.m_interval;
            return *this;
        }
        Interval& operator+=(const Interval& rhs)
        {
            m_interval += rhs.m_interval;
            return *this;
        }
        ENUM_NURBS set_interval(T low, T high)
        {
            if (low > high)
                return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
            m_interval[0] = low;
            m_interval[1] = high;
        }
        bool contain(T number, T tol) const
        {
            if (m_interval[0] - tol > number || m_interval[1] + tol < number)
            {
                return false;
            }
            return true;
        }
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

        bool is_contain_point(const Eigen::Vector<T, dim> &point, T eps = 0) const
        {
            for (int index = 0; index < dim; ++index)
            {
                if (point[index] > Max[index] + eps || point[index] < Min[index] - eps)
                    return false;
            }
            return true;
        }

        bool is_contain_box(const Box &bx, const T &tol = PRECISION<T>::value) const
        {
            for (int index = 0; index < dim; ++index)
            {
                if (bx.Min[index] < Min[index] - tol || bx.Max[index] > Max[index] + tol)
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
        std::vector<Box<T, dim>> split_at_middle() const
        {
            std::vector<Box<T, dim>> sub_boxes;
            sub_boxes.push_back(*this);
            std::vector<Box<T, dim>> next_boxes;
            for (int index = 0; index < dim; ++index)
            {
                for (int box_index = 0; box_index < sub_boxes.size(); ++box_index)
                {
                    std::array<Box<T, dim>, 2> temp_result = sub_boxes[box_index].split_at_middle(index);
                    next_boxes.insert(next_boxes.end(), temp_result.begin(), temp_result.end());
                }
                std::swap(next_boxes, sub_boxes);
                next_boxes.clear();
            }
            return sub_boxes;
        }
        
        Eigen::Matrix2<Box<T, 2>> split_at_param(const Eigen::Vector<T, dim>& param) const
        {
            static_assert(dim == 2, "dim != 2");
            Eigen::Matrix2<Box<T, 2>> sub_boxes;
            Eigen::Vector<T, 2> dir = (Max - Min) * 0.5;
            sub_boxes(0, 0).Min = Min;
            sub_boxes(0, 0).Max = param;

            sub_boxes(0, 1).Min = Eigen::Vector2<T>(param[0], Min[1]);
            sub_boxes(0, 1).Max = Eigen::Vector2<T>(Max[0], param[1]);
            
            
            sub_boxes(1, 0).Min = Eigen::Vector2<T>(Min[0], param[1]);
            sub_boxes(1, 0).Max = Eigen::Vector2<T>(param[0], Max[1]);
            
            
            sub_boxes(1, 1).Min = param;
            sub_boxes(1, 1).Max = Max;
            return sub_boxes;
        }
        bool intersect(const Box &box, Box &intersect_box) const
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
            for (int index = 0; index < dim; ++index)
            {
                T l1 = Min[index] * s;
                T l2 = Max[index] * s;
                if (l1 > l2)
                {
                    Min[index] = l2;
                    Max[index] = l1;
                }
                else
                {
                    Min[index] = l1;
                    Max[index] = l2;
                }
            }
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

        T get_min_lenght() const
        {
            T min_length = Max[0] - Min[0];
            for (int index = 1; index < dim; ++index)
            {
                min_length = std::min(min_length, Max[index] - Min[index]);
            }
            return min_length;
        }

        Box<T, 1> get_index_interval(int index) const
        {
            Box<T, 1> result;
            result.Min[0] = Min[index];
            result.Max[0] = Max[index];
            return result;
        }

        std::array<Box<T, 1>, dim> sperate_box() const
        {
            std::array<Box<T, 1>, dim> result;
            for (int index = 0; index < dim; ++index)
            {
                result[index].Min[0] = Min[index];
                result[index].Max[0] = Max[index];
            }
            return result;
        }



        bool set_index_interval(int index, const Box<T, 1> &interval)
        {
            Min[index] = interval.Min[0];
            Max[index] = interval.Max[0];
            return true;
        }
    
        Box<T, dim> operator+(const Box<T, dim>& right) const
        {
            Box<T, dim> result;
            result.Min = Min + right.Min;
            result.Max = Max + right.Max;
            return result;
        }

        Box<T, dim> operator-(const Box<T, dim>& right) const
        {
            Box<T, dim> result;
            result.Min = Min - right.Max;
            result.Max = Max - right.Min;
            return result;
        }


        template<int n>
        Box<T, dim> operator*(const Box<T, n>& right) const
        {
            static_assert(n == 1 || n == dim, "n != 1 && n != dim");
            Box<T, dim> result;
            if constexpr (n == dim)
            {
                for (int index = 0; index < dim; ++index)
                {
                    if (Min[index] > 0 && right.Min[index] > 0)
                    {
						T x1 = Min[index] * right.Min[index];
						T x4 = Max[index] * right.Max[index];
						result.Min[index] = x1;
						result.Max[index] = x4;
                    }
                    else if (Max[index] < 0 && right.Max[index] < 0)
                    {
						T x1 = Min[index] * right.Min[index];
						T x4 = Max[index] * right.Max[index];
						result.Min[index] = x4;
						result.Max[index] = x1;
                    }
                    else
                    {
						T x1 = Min[index] * right.Min[index];
						T x2 = Min[index] * right.Max[index];
						T x3 = Max[index] * right.Min[index];
						T x4 = Max[index] * right.Max[index];
						auto min_max_element = std::minmax({ x1, x2, x3, x4 });
						result.Min[index] = min_max_element.first;
						result.Max[index] = min_max_element.second;
                    }
                    // result.Min[index] = std::min({ x1, x2, x3, x4 });
                    // result.Max[index] = std::max({ x1, x2, x3, x4 });
                }
                return result;
            }
            else
            {
                for (int index = 0; index < dim; ++index)
                {
                    T x1 = Min[index] * right.Min[0];
                    T x2 = Min[index] * right.Max[0];
                    T x3 = Max[index] * right.Min[0];
                    T x4 = Max[index] * right.Max[0];
                    auto min_max_element = std::minmax({ x1, x2, x3, x4 });
                    result.Min[index] = min_max_element.first;
                    result.Max[index] = min_max_element.second;
                    // result.Min[index] = std::min({ x1, x2, x3, x4 });
                    // result.Max[index] = std::max({ x1, x2, x3, x4 });
                }
                return result;
            }
            
        }

        Box<T, dim> operator*(const T scale) const
        {
            Box<T, dim> result;
            for (int index = 0; index < dim; ++index)
            {
                T x1 = scale * Min[index];
                T x2 = scale * Max[index];

                result.Min[index] = std::min({ x1, x2 });
                result.Max[index] = std::max({ x1, x2 });
            }
            return result;
        }

        Box<T, dim> operator*(const Interval<T> &scale) const
        {
            Box<T, dim> result;
            for (int index = 0; index < dim; ++index)
            {
                T x1 = scale.m_interval[0] * Min[index];
                T x2 = scale.m_interval[0] * Max[index];
                T x3 = scale.m_interval[1] * Min[index];
                T x4 = scale.m_interval[1] * Max[index];

                result.Min[index] = std::min({ x1, x2, x3, x4 });
                result.Max[index] = std::max({ x1, x2, x3, x4 });
            }
            return result;
        }

        void enlarge(const Eigen::Vector<T, dim>& point)
        {
            for (int index = 0; index < dim; ++index)
            {
                Min[index] = std::min(Min[index], point[index]);
                Max[index] = std::max(Max[index], point[index]);
            }
            return;
        }


        bool unit(Box<T, dim>& point)
        {
            enlarge(point.Min);
            enlarge(point.Max);
            return true;
        }

        void move_vector(const Eigen::Vector<T, dim>& vec)
        {
            Min += vec;
            Max += vec;
            return;
        }



};
       template<typename T, int dim> 
       Box<T, dim> envelop_points(std::vector<Eigen::Vector<T, dim>>& points)
       {
           Box<T, dim> box;
           T max = std::numeric_limits<T>::max();
           box.Min.setConstant(max);
           box.Max.setConstant(-max);
           for (const auto point : points)
           {
               box.enlarge(point);
           }
           return box;
       };
    // template<typename T, int dim>
    // struct Eigen::NumTraits<Box<T, dim>> : Eigen::NumTraits<T>
    // {

    // };

    template<typename T, int dim>
    ENUM_NURBS create_box(const Eigen::Vector<T, dim>& center, T len, Box<T, dim> &box)
    {
        if (len < 0)
            return ENUM_NURBS::NURBS_ERROR;
        Eigen::Vector<T, dim> len_vec;
        len_vec.setConstant(len);
        box.Min = center - len_vec;
        box.Max = center + len_vec;
        return ENUM_NURBS::NURBS_SUCCESS;
    }

}
namespace Eigen
{
    template<typename Scalar, int dim>
    struct NumTraits<tnurbs::Box<Scalar, dim>> 
    {
        typedef tnurbs::Box<Scalar, dim> BoxType;
        typedef typename NumTraits<Scalar>::Real RealScalar;
        typedef tnurbs::Box<RealScalar, dim> Real;
        typedef typename NumTraits<Scalar>::NonInteger NonIntegerScalar;
        typedef tnurbs::Box<NonIntegerScalar, dim> NonInteger;
        typedef BoxType & Nested;
        typedef typename NumTraits<Scalar>::Literal Literal;

        enum {
          IsComplex = NumTraits<Scalar>::IsComplex,
          IsInteger = NumTraits<Scalar>::IsInteger,
          IsSigned  = NumTraits<Scalar>::IsSigned,
          RequireInitialization = 1,
          ReadCost = HugeCost,
          AddCost  = HugeCost,
          MulCost  = HugeCost
        };

         EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
         static inline RealScalar epsilon() { return NumTraits<RealScalar>::epsilon(); }
         EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
         static inline RealScalar dummy_precision() { return NumTraits<RealScalar>::dummy_precision(); }
    
         EIGEN_CONSTEXPR
         static inline int digits10() { return NumTraits<Scalar>::digits10(); }
    };   
}
