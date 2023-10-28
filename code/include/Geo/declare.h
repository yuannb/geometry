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

}

