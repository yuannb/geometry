
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "convert_nubrs_with_polynomial.h" 
#include "create_nurbs_arc.h"
#include "contruct_primitive_nurbs_surface.h"
#include "fit_nurbs.h"
#include "gtest/gtest.h"
#include <cmath>
#include "debug_used.h"
#include "nurbs_build.h"
#include "discret.h"
#include "nearest.h"
// #include <memory>
#include "convex_hell.h"
#include "bezier_curve_int.h"

#include <chrono>
using namespace tnurbs;

//vs��assert��ֹ���Բ鿴��ջ�ڴ�������������ĺ���
//_set_error_mode(_OUT_TO_MSGBOX);
const std::string src_path = "../intersectData/";

using nunbs_curve3d = nurbs_curve<double, 3, false, -1, -1>;
using nurbs_curve3d = nurbs_curve<double, 3, true, -1, -1>;
using nunbs_surface3d = nurbs_surface<double, 3, -1, -1, -1, -1, false>;
using nurbs_surface3d = nurbs_surface<double, 3, -1, -1, -1, -1, true>;
TEST(BEZIER_INT, test2)
{
    _set_error_mode(_OUT_TO_MSGBOX);
    Eigen::VectorX<double> u_knots_vector2(8);
    u_knots_vector2<< 0, 0, 0, 0, 1, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector2(8);
    v_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;

    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector<< 0, 0, 0, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 1, 1, 1;


    Eigen::Matrix<double, 3, Eigen::Dynamic> points1(3, 3);
    points1 << -10, 0, 10,
               -30, -30, -30,
               3, -3, 3;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points2(3, 3);
    points2 << -10, 0, 10,
               -10, -10, -10,
               3, -3, 3;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points3(3, 3);
    points3 << -10, 0, 10,
               10, 10, 10,
               3, -3, 3;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points4(3, 3);
    points4 << -10, 0, 10,
               30,  30, 15,
               3,  -3,  3;

    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points[0] = points1;
    control_points[1] = points2;
    control_points[2] = points3;
    control_points[3] = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false> *test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);
    
    Eigen::Matrix<double, 3, Eigen::Dynamic> points11(3, 4);
    points11 << 20,  10, -10, -20,
                -10, -10, -10, -10,
                3,  3, 3, 3;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points12(3, 4);
    points12 << 20,  10, -10, -20,
                -10, -10, -10, -10,
                -1,  -1, -1, -1;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points13(3, 4);
    points13 << 20,  10, -10, -20,
                10, 10, 10, 10,
                -1,  -1, -1, -1;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points14(3, 4);
    points14 << 20,  10, -10, -20,
                10, 10,  10,  10,
                3,  3, 3, 3;

    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2[0] = points11;
    control_points2[1] = points12;
    control_points2[2] = points13;
    control_points2[3] = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false> *test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    // save_obj2(*test_surface, "test_surface.obj");
    // save_obj2(*test_surface2, "test_surface2.obj");

    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface, test_surface2);
	using namespace std::chrono;
	auto start = steady_clock::now();
    ts.surafces_intersection2();
	auto end = steady_clock::now();
	auto last_time = duration_cast<microseconds>(end - start);

	std::cout << last_time.count() << "um" << std::endl;
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest2.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    // save_chat_points_file(intersection, "intersectcurve.obj");
    surf_surf_int<double, 3> test_data;
    const std::string path("D:\\geometry\\intersectData\\BEZIER_INTtest2.json");
    // read_chat_points(test_data, src_path + "BEZIER_INTtest2.json");
    read_chat_points(test_data, path);
    EXPECT_EQ(curve_count, test_data.m_int_chats.size());
    EXPECT_EQ(isolate_point_count, test_data.m_isolate_points.size());

    for (size_t index = 0; index < curve_count; ++index)
    {
        surfs_int_points_chat<double, 3>& int_chat = test_data.m_int_chats[index];
        size_t curve_index;
        bool flag = ts.is_point_in_intcurve(int_chat.m_inter_points[1].m_uv, int_chat.m_is_transversal, curve_index);
        EXPECT_TRUE(flag);
        size_t points_count = int_chat.m_inter_points.size();
        EXPECT_TRUE(points_count < 1.2 * intersection.m_int_chats[curve_index].m_inter_points.size());
        for (size_t i = 2; i < points_count; ++i)
        {
            size_t current_index;
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[i].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}  

TEST(BEZIER_INT, test3)
{
    Eigen::VectorX<double> u_knots_vector2(8);
    u_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector2(8);
    v_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;

    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector << 0, 0, 0, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 1, 1, 1;


    Eigen::Matrix<double, 3, Eigen::Dynamic> points1(3, 3);
    points1 << -150, 0, 150,
        10, -10, 10,
        0, 0, 0;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points2(3, 3);
    points2 << -10, 20, 50,
        40, 0, 40,
        50, 50, 50;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points3(3, 3);
    points3 << -10, 20, 50,
        40, 0, 40,
        110, 110, 110;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points4(3, 3);
    points4 << -55, -25, 5,
        -5, -45, -5,
        200, 200, 200;

    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points[0] = points1;
    control_points[1] = points2;
    control_points[2] = points3;
    control_points[3] = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);
    Eigen::Matrix<double, 3, Eigen::Dynamic> points11(3, 3);
    points11 << -50, 0, 50,
        -50, 50, -50,
        0, 0, 0;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points12(3, 3);
    points12 << -30, 20, 70,
        -30, 70, -30,
        50, 50, 50;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points13(3, 3);
    points13 << -30, 20, 70,
        -30, 70, -30,
        110, 110, 110;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points14(3, 3);
    points14 << -75, -25, 25,
        -75, 25, -75,
        200, 200, 200;

    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2[0] = points11;
    control_points2[1] = points12;
    control_points2[2] = points13;
    control_points2[3] = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points2);
    save_obj2(*test_surface, "test_surface.obj");
    save_obj2(*test_surface2, "test_surface2.obj");
    _set_error_mode(_OUT_TO_MSGBOX);
    // trace_nurbs_surface<nurbs_surface<double, 3, -1, -1, -1, -1, false>> ts;
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface, test_surface2);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest3.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    const std::string path("D:\\geometry\\intersectData\\BEZIER_INTtest3.json");
    read_chat_points(test_data, path);
    // read_chat_points(test_data, src_path + "BEZIER_INTtest3.json");
    EXPECT_EQ(curve_count, test_data.m_int_chats.size());
    EXPECT_EQ(isolate_point_count, test_data.m_isolate_points.size());

    for (size_t index = 0; index < curve_count; ++index)
    {
        surfs_int_points_chat<double, 3>& int_chat = test_data.m_int_chats[index];
        size_t curve_index;
        bool flag = ts.is_point_in_intcurve(int_chat.m_inter_points[1].m_uv, int_chat.m_is_transversal, curve_index);
        EXPECT_TRUE(flag);
        size_t points_count = int_chat.m_inter_points.size();
        EXPECT_TRUE(points_count < 1.2 * intersection.m_int_chats[curve_index].m_inter_points.size());
        for (size_t i = 2; i < points_count; ++i)
        {
            size_t current_index;
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[i].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}


TEST(BEZIER_INT, test4)
{
    Eigen::VectorX<double> u_knots_vector2(8);
    u_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector2(8);
    v_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;

    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector << 0, 0, 0, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 1, 1, 1;


    Eigen::Matrix<double, 3, Eigen::Dynamic> points1(3, 3);
    points1 << -10, 0, 10,
        -30, -30, -30,
        2.5, -3.5, 2.5;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points2(3, 3);
    points2 << -10, 0, 10,
        -10, -10, -10,
        2.5, -3.5, 2.5;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points3(3, 3);
    points3 << -10, 0, 10,
        10, 10, 10,
        2.5, -3.5, 2.5;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points4(3, 3);
    points4 << -10, 0, 10,
        30, 30, 15,
        2.5, -3.5, 2.5;

    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points[0] = points1;
    control_points[1] = points2;
    control_points[2] = points3;
    control_points[3] = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);

    Eigen::Matrix<double, 3, Eigen::Dynamic> points11(3, 4);
    points11 << 20, 10, -10, -20,
        -10, -10, -10, -10,
        3, 3, 3, 3;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points12(3, 4);
    points12 << 20, 10, -10, -20,
        -10, -10, -10, -10,
        -1, -1, -1, -1;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points13(3, 4);
    points13 << 20, 10, -10, -20,
        10, 10, 10, 10,
        -1, -1, -1, -1;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points14(3, 4);
    points14 << 20, 10, -10, -20,
        10, 10, 10, 10,
        3, 3, 3, 3;

    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2[0] = points11;
    control_points2[1] = points12;
    control_points2[2] = points13;
    control_points2[3] = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    _set_error_mode(_OUT_TO_MSGBOX);

    // save_obj2(*test_surface, "test_surface.obj");
    // save_obj2(*test_surface2, "test_surface2.obj");
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface2, test_surface);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    //save_chat_points(intersection, "BEZIER_INTtest4.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    const std::string path("D:\\geometry\\intersectData\\BEZIER_INTtest4.json");
    read_chat_points(test_data, path);
    EXPECT_EQ(curve_count, test_data.m_int_chats.size());
    EXPECT_EQ(isolate_point_count, test_data.m_isolate_points.size());

    for (size_t index = 0; index < curve_count; ++index)
    {
        surfs_int_points_chat<double, 3>& int_chat = test_data.m_int_chats[index];
        size_t curve_index;
        bool flag = ts.is_point_in_intcurve(int_chat.m_inter_points[1].m_uv, int_chat.m_is_transversal, curve_index);
        EXPECT_TRUE(flag);
        size_t points_count = int_chat.m_inter_points.size();
        EXPECT_TRUE(points_count < 1.2 * intersection.m_int_chats[curve_index].m_inter_points.size());
        for (size_t i = 2; i < points_count; ++i)
        {
            size_t current_index;
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[i].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}
TEST(BEZIER_INT, test5)
{
    Eigen::VectorX<double> u_knots_vector2(8);
    u_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector2(8);
    v_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;

    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector << 0, 0, 0, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 1, 1, 1;


    Eigen::Matrix<double, 3, Eigen::Dynamic> points1(3, 3);
    points1 << 0, 10, 20,
        0, 0, 0,
        -46.45, -46.45, -46.45;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points2(3, 3);
    points2 << 0, 10, 20,
        10, 10, 10,
        -46.45, -46.45, -46.45;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points3(3, 3);
    points3 << 0, 10, 20,
        20, 20, 20,
        -46.45, -46.45, -46.45;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points4(3, 3);
    points4 << 0, 10, 20,
        30, 30, 30,
        -46.45, -46.45, -46.45;
    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points[0] = points1;
    control_points[1] = points2;
    control_points[2] = points3;
    control_points[3] = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);

    Eigen::Matrix<double, 3, Eigen::Dynamic> points11(3, 4);
    points11 << 0, 2, 4, 6,
        0, 10, 10, 0,
        -200, -200, -200, -200;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points12(3, 4);
    points12 << 0, 2, 4, 6,
        13, 17, 17, 13,
        5, 10, 15, 20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points13(3, 4);
    points13 << 0, 2, 4, 6,
        17, 13, 13, 17,
        5, 10, 15, 20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points14(3, 4);
    points14 << 0, 2, 4, 6,
        15, 10, 10, 15,
        -200, -200, -200, -200;

    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2[0] = points11;
    control_points2[1] = points12;
    control_points2[2] = points13;
    control_points2[3] = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    int i = 0;
    _set_error_mode(_OUT_TO_MSGBOX);
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    // save_obj2(*test_surface, "test_surface.obj");
    // save_obj2(*test_surface2, "test_surface2.obj");
    ts.init(test_surface, test_surface2);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest5.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    const std::string path("D:\\geometry\\intersectData\\BEZIER_INTtest5.json");
    read_chat_points(test_data, path);
    EXPECT_EQ(curve_count, test_data.m_int_chats.size());
    EXPECT_EQ(isolate_point_count, test_data.m_isolate_points.size());

    for (size_t index = 0; index < curve_count; ++index)
    {
        surfs_int_points_chat<double, 3>& int_chat = test_data.m_int_chats[index];
        size_t curve_index;
        bool flag = ts.is_point_in_intcurve(int_chat.m_inter_points[1].m_uv, int_chat.m_is_transversal, curve_index);
        EXPECT_TRUE(flag);
        size_t points_count = int_chat.m_inter_points.size();
        EXPECT_TRUE(points_count < 1.2 * intersection.m_int_chats[curve_index].m_inter_points.size());
        for (size_t i = 2; i < points_count; ++i)
        {
            size_t current_index;
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[i].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}

TEST(BEZIER_INT, test6)
{
    Eigen::VectorX<double> u_knots_vector2(8);
    u_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector2(8);
    v_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;

    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector << 0, 0, 0, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 1, 1, 1;


    Eigen::Matrix<double, 3, Eigen::Dynamic> points1(3, 3);
    points1 << 0, 10, 20,
        0, 0, 0,
        1.57, 1.57, 1.57;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points2(3, 3);
    points2 << 0, 10, 20,
        10, 10, 10,
        1.57, 1.57, 1.57;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points3(3, 3);
    points3 << 0, 10, 20,
        20, 20, 20,
        1.57, 1.57, 1.57;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points4(3, 3);
    points4 << 0, 10, 20,
        30, 30, 30,
        1.57, 1.57, 1.57;
    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points[0] = points1;
    control_points[1] = points2;
    control_points[2] = points3;
    control_points[3] = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);

    Eigen::Matrix<double, 3, Eigen::Dynamic> points11(3, 4);
    points11 << 0, 4, 8, 12,
        0, 10, 10, 0,
        -20, -20, -20, -20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points12(3, 4);
    points12 << 0, 0, 0, 12,
        -50, 12, 12, -50,
        20, 5, 5, 20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points13(3, 4);
    points13 << 0, 12, 12, 12,
        70, 8, 8, 70,
        20, 5, 5, 20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points14(3, 4);
    points14 << 0, 4, 8, 12,
        20, 10, 10, 20,
        -20, -20, -20, -20;

    std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2[0] = points11;
    control_points2[1] = points12;
    control_points2[2] = points13;
    control_points2[3] = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    _set_error_mode(_OUT_TO_MSGBOX);
    // save_obj2(*test_surface, "test_surface.obj");
    // save_obj2(*test_surface2, "test_surface2.obj");
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface2, test_surface);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest6.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    const std::string path("D:\\geometry\\intersectData\\BEZIER_INTtest6.json");
    read_chat_points(test_data, path);
    EXPECT_EQ(curve_count, test_data.m_int_chats.size());
    EXPECT_EQ(isolate_point_count, test_data.m_isolate_points.size());

    for (size_t index = 0; index < curve_count; ++index)
    {
        surfs_int_points_chat<double, 3>& int_chat = test_data.m_int_chats[index];
        size_t curve_index;
        bool flag = ts.is_point_in_intcurve(int_chat.m_inter_points[1].m_uv, int_chat.m_is_transversal, curve_index);
        EXPECT_TRUE(flag);
        size_t points_count = int_chat.m_inter_points.size();
        EXPECT_TRUE(points_count < 1.2 * intersection.m_int_chats[curve_index].m_inter_points.size());
        for (size_t i = 2; i < points_count; ++i)
        {
            size_t current_index;
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[i].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}



int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    
    ::testing::FLAGS_gtest_filter = "BEZIER_INT.test2";
    return RUN_ALL_TESTS();
}
