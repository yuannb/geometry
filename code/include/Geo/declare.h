#pragma once
namespace tnurbs
{
    using namespace tnurbs;
    template<typename T> struct geo_traits;

    template<typename T> struct geo_traits<const T> : geo_traits<T> {};
}

