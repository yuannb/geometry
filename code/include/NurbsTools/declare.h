#pragma once
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

}

