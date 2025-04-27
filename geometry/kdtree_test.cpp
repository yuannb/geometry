#include "kdtree.h"
#include "gtest/gtest.h"
#include <stdlib.h>
#include <random>
#include <iostream>
#include "debug_used.h"

using namespace tnurbs;

class kdtree_test : public testing::Test
{
};

void generate_points(int count, double max)
{
    std::random_device rd;
    unsigned int seed = rd();
    std::minstd_rand0 gen(seed);
    std::vector<Eigen::Vector2d> points;
    for (int index = 0; index < count; ++index)
    {
        Eigen::Vector2d p;
        for (int j = 0; j < 2; ++j)
        {
			double num = (double)gen();
			num -= int(num / max) * max;
            p[j] = num;
        }
        points.push_back(p);
    }
    const std::string path("D:\\geometry\\intersectData\\points1.json");
    save_points_to_json(points, path.data());
}


TEST_F(kdtree_test, test1)
{
    // generate_points(10, 100);
    const std::string path("D:\\geometry\\intersectData\\points1.json");
    std::vector<Eigen::Vector2d> points;
    read_points_from_json(points, path);
    Eigen::Vector2d test_point{ 1.313, 34.4234 };
    double dis1 = (points[0] - test_point).norm();
    int min_index = 0;
    for (int index = 1; index < 10; ++index)
    {
		double dis2 = (points[index] - test_point).norm();
        if (dis2 < dis1)
        {
            dis1 = dis2;
            min_index = index;
        }
    }

    KDTree<double, 10> test_kdtree;
    for (int index = 0; index < 10; ++index)
    {
        test_kdtree.insert(index, points);
    }
    KDTree<double, 10>::point_index j = test_kdtree.nearest(test_point, points);
    EXPECT_EQ(int(j), min_index);
    std::cout << "j : " << j << std::endl;
    std::cout << "min_index : " << min_index << std::endl;

}
TEST_F(kdtree_test, test2)
{
    // generate_points(10, 100);
    const std::string path("D:\\geometry\\intersectData\\points1.json");
    std::vector<Eigen::Vector2d> points;
    read_points_from_json(points, path);
    Eigen::Vector2d test_point{ 131.3, 74.4234 };
    double dis1 = (points[0] - test_point).norm();
    int min_index = 0;
    for (int index = 1; index < 10; ++index)
    {
		double dis2 = (points[index] - test_point).norm();
        if (dis2 < dis1)
        {
            dis1 = dis2;
            min_index = index;
        }
    }

    KDTree<double, 10> test_kdtree;
    for (int index = 0; index < 10; ++index)
    {
        test_kdtree.insert(index, points);
    }
    KDTree<double, 10>::point_index j = test_kdtree.nearest(test_point, points);
    EXPECT_EQ(int(j), min_index);
    std::cout << "j : " << j << std::endl;
    std::cout << "min_index : " << min_index << std::endl;

}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    
    ::testing::FLAGS_gtest_filter = "kdtree_test.test2";
    return RUN_ALL_TESTS();
}