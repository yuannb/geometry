#pragma once
// #include "declare.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
namespace tnurbs
{
    //find the point with minimum y coordinate in case of tie choose the point with minimun x-coordinate
    template<typename T = double>
    size_t find_min_ycoord(const std::vector<Eigen::Vector2<T>> &points)
    {
        size_t points_count = points.size();
        size_t min_ycoord_index = 0;
        for (size_t index = 1; index < points_count; ++index)
        {
            if (points[index][1] < points[min_ycoord_index][1])
                min_ycoord_index = index;
            else if (points[index][1] == points[min_ycoord_index][1])
                if (points[index][0] < points[min_ycoord_index][0])
                    min_ycoord_index = index;
        }
        return min_ycoord_index;
    }

    template<typename T = double>
    std::vector<Eigen::Vector2<T>> graham_scan(std::vector<Eigen::Vector2<T>> &points)
    {
        size_t points_count = points.size();
        if (points_count <= 2)
            return { };
        size_t index = find_min_ycoord(points);
        std::swap(points[index], points[0]);
        Eigen::Vector2<T> start_point = points[0];
        std::sort(points.begin() + 1, points.end(), [&start_point](Eigen::Vector2<T> first, Eigen::Vector2<T> second) -> bool
        {
            Eigen::Vector2<T> v1 = first - start_point;
            Eigen::Vector2<T> v2 = second - start_point;
            T d = v1[0] * v2[1] - v1[1] * v2[0];
            if (d < -PRECISION<T>::value)
                return false;
            else if (d > PRECISION<T>::value)
                return true;
            else
            {
                T d1 = v1.squaredNorm();
                T d2 = v2.squaredNorm();
                if (d1 < d2)
                    return true;
                return false;
            }
            return true;
        });
        //remove repeat elem
        std::vector<Eigen::Vector2<T>> new_points;
        new_points.reserve(points_count);
        new_points.push_back(points[0]);
        new_points.push_back(points[1]);
        for (auto it = points.begin() + 2; it != points.end(); ++it)
        {
            Eigen::Vector2<T> dir1 = new_points.back() - new_points[0];
            Eigen::Vector2<T> dir2 = *it - new_points.back();
            T d = dir1[0] * dir2[1] - dir1[1] * dir2[0];
            if (d != 0)
                new_points.push_back(*it);
            else if (d == 0)
            {
                new_points.back() = *it;
            }
        }

        size_t new_points_count = new_points.size();
        if (new_points_count == 1)
            return std::vector<Eigen::Vector2<T>>{ new_points[0] };
        if (new_points_count == 2)
            return std::vector<Eigen::Vector2<T>> { new_points[0], new_points[1] };
        if (new_points_count == 3)
            return std::vector<Eigen::Vector2<T>> { new_points[0], new_points[1], new_points[2] };

        std::vector<Eigen::Vector2<T>> convex_hell{ new_points[0], new_points[1], new_points[2] };
        for (auto it = new_points.begin() + 3; it != new_points.end(); ++it)
        {
            size_t convex_hell_count = convex_hell.size();
            while (convex_hell_count > 2)
            {
                Eigen::Vector2<T> dir1 = convex_hell[convex_hell_count - 1] - convex_hell[convex_hell_count - 2];
                Eigen::Vector2<T> dir2 = *it - convex_hell.back();
                T d = dir1[0] * dir2[1] - dir1[1] * dir2[0];
                if (d > 0)
                    break;
                else
                {
                    convex_hell_count -= 1;
                    convex_hell.erase(convex_hell.end() - 1);
                }
            }
            convex_hell.push_back(*it);
        }
        return convex_hell;
    }
}

