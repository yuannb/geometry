#pragma once
// #include "declare.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
// #include <chrono>
// using namespace std::chrono;
namespace tnurbs
{

    inline void printx(std::vector<Eigen::Vector2<double>>& points)
    {
        for (auto& p : points)
        {
            std::cout << p << std::endl;
            std::cout << "EEEEEE" << std::endl;
        }
		std::cout << "************************" << std::endl;
    }

    inline void printx(std::vector<Eigen::Vector2<double>>& points, std::vector<int>& indices)
    {
        for (auto& i : indices)
        {
            std::cout << points[i] << std::endl;
            std::cout << "EEEEEE" << std::endl;
        }
		std::cout << "************************" << std::endl;
    }

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
    bool can_add_one_element(std::vector<Eigen::Vector2<T>>& convex_hell, const Eigen::Vector2<T>& point)
    {
        size_t count = convex_hell.size();
        if (convex_hell.size() == 0)
        {
            return true;
        }
        for (size_t index = 0; index < count; ++index)
        {
            if ((convex_hell[index] - point).norm() < 1e-11)
            {
                return false;
            }
        }
        return true;
    }



    template<typename T = double>
    void graham_scan(std::vector<Eigen::Vector2<T>> &points, std::vector<int>& indices, std::vector<Eigen::Vector2<T>>& convex_hell)
    {
        size_t points_count = points.size();
        if (points_count <= 2)
            return;

        size_t index = find_min_ycoord(points);
        std::swap(points[index], points[0]);
        Eigen::Vector2<T> start_point = points[0];
        Eigen::Vector2<T> v1;
		Eigen::Vector2<T> v2;
		
        // std::vector<int> indices(points_count);
		// for (int i = 0; i < points_count; ++i)
		// {
		// 	indices[i] = i;
		// }

        std::vector<double> vec_dots(points_count);
        vec_dots[0] = 1.0;
        std::vector<double> vec_norms(points_count);
        vec_norms[0] = start_point.norm();


        auto compareFun = [&start_point, &v1, &v2](const Eigen::Vector2<T>& first, const Eigen::Vector2<T>& second)
            {
                v1 = first - start_point;
                if (v1.norm() < 1e-8)
                {
                    return false;
                }
                
                v2 = second - start_point;
                if (v2.norm() < 1e-8)
                {
                    return true;
                }
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

            };


        Eigen::Vector2<T> dir{ 1, 0 };
        for (int index = 1; index < points_count; ++index)
        {
            v1 = points[index] - start_point;
            vec_norms[index] = v1.norm();
            if (vec_norms[index] < 1e-11)
            {
                vec_dots[index] = 0;
                vec_norms[index] = 0;
            }
            double is_positive = v1[0] * dir[0] + v1[1] * dir[1];
            vec_dots[index] = is_positive / vec_norms[index];
        }


        std::sort(indices.begin() + 1, indices.end(), [&vec_dots, &vec_norms](int i, int j) -> bool
        {
                if (vec_dots[i] > vec_dots[j])
                {
                    return true;
                }
                else if (vec_dots[i] < vec_dots[j])
                {
                    return false;
                }

                if (vec_norms[i] < vec_norms[j])
                {
                    return true;
                }
                return false;

        });

        // convex_hell.reserve(points_count);
        
        convex_hell.push_back(points[0]);
        size_t convex_hell_count = 1;
        size_t point_index = 1;
        for (point_index = 1; point_index < points_count; ++point_index)
        {
            if (can_add_one_element(convex_hell, points[indices[point_index]]) == true)
            {
                convex_hell_count += 1;
                convex_hell.push_back(points[indices[point_index]]);
            }
            if (convex_hell_count > 2)
            {
                break;
            }
        }

        if (convex_hell_count <= 2)
        {
            return;
        }

        for (auto it = indices.begin() + point_index; it != indices.end(); ++it)
        {
            convex_hell_count = convex_hell.size();
            bool flag = true;
            while (convex_hell_count > 2)
            {
                Eigen::Vector2<T> dir1 = convex_hell[convex_hell_count - 1] - convex_hell[convex_hell_count - 2];
                Eigen::Vector2<T> dir2 = points[*it] - convex_hell.back();
                if (dir2.norm() < 1e-11)
                {
                    flag = false;
                    // convex_hell.pop_back();
                    convex_hell.back() == points[*it];
                    break;
                }
                dir1.normalize();
                dir2.normalize();
                T d = dir1[0] * dir2[1] - dir1[1] * dir2[0];
                if (std::abs(d) < 1e-11)
                {
                    break;
                    // if (dir1.dot(dir2) > 0)
                    // {
                    //     break;
                    // }
                    // else
                    // {
                    //     convex_hell_count -= 1;
                    //     convex_hell.pop_back();
                    // }
                }
                else if (d > 0)
                    break;
                else
                {
                    convex_hell_count -= 1;
                    convex_hell.pop_back();
                    // convex_hell.erase(convex_hell.end() - 1);
                }
            }
            if (flag == true)
            {
                convex_hell.push_back(points[*it]);
            }
        }
        return;
    }

    template<typename T = double>
    void graham_scan(std::vector<Eigen::Vector2<T>>& points, std::vector<Eigen::Vector2<T>>& convex_hell)
    {
        convex_hell.clear();
        size_t points_count = points.size();
        if (points_count <= 2)
            return;
        size_t index = find_min_ycoord(points);
        std::swap(points[index], points[0]);
        Eigen::Vector2<T> start_point = points[0];
        Eigen::Vector2<T> v1;
        Eigen::Vector2<T> v2;
        std::sort(points.begin() + 1, points.end(), [&start_point, &v1, &v2](const Eigen::Vector2<T>& first, const Eigen::Vector2<T>& second) -> bool
        {
            // Eigen::Vector2<T> v1 = first - start_point;
            // Eigen::Vector2<T> v2 = second - start_point;
                v1 = first - start_point;
                v2 = second - start_point;
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
        new_points.push_back(std::move(points[0]));
        new_points.push_back(std::move(points[1]));
        for (auto it = points.begin() + 2; it != points.end(); ++it)
        {
            Eigen::Vector2<T> dir1 = new_points.back() - new_points[0];
            Eigen::Vector2<T> dir2 = *it - new_points.back();
            T d = dir1[0] * dir2[1] - dir1[1] * dir2[0];

            if (std::abs(d) < PRECISION<T>::value)
            {
                new_points.back() = std::move(*it);
            }
            else
                new_points.push_back(std::move(*it));
        }

        size_t new_points_count = new_points.size();
        convex_hell.reserve(new_points_count);
        if (new_points_count >= 1)
            convex_hell.push_back(std::move(new_points[0]));
        if (new_points_count >= 2)
            convex_hell.push_back(std::move(new_points[1]));
        if (new_points_count >= 3)
            convex_hell.push_back(std::move(new_points[2]));

        if (new_points_count <= 3)
        {
            return;
        }

        // std::vector<Eigen::Vector2<T>> convex_hell{ new_points[0], new_points[1], new_points[2] };
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
        // return convex_hell;
        return;
    }



}

