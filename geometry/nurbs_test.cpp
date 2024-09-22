
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
using namespace tnurbs;

//vs中assert中止调试查看堆栈在代码中添加下面的函数
//_set_error_mode(_OUT_TO_MSGBOX);
const std::string src_path = "../../intersectData/";

using nunbs_curve3d = nurbs_curve<double, 3, false, -1, -1>;
using nurbs_curve3d = nurbs_curve<double, 3, true, -1, -1>;
using nunbs_surface3d = nurbs_surface<double, 3, -1, -1, -1, -1, false>;
using nurbs_surface3d = nurbs_surface<double, 3, -1, -1, -1, -1, true>;

class CreateNurbsCurve : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::Vector<double, 3> center{0, 0, 0};
        Eigen::Vector<double, 3> u_dir{1, 0, 0};
        Eigen::Vector<double, 3> v_dir{0, 1, 0};
        double radius = 10;
        double start_angles = 0;
        double end_angles = M_PI + 0.78;
        nurbs_curve<double, 3, true, -1, -1> nurbs;
        pointss.resize(3, 5);
        ders.resize(3, 5);
        create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, nurbs);
        for (int i = 0; i < 5; ++i)
        {
            Eigen::Vector<Eigen::Vector3d, 2> point;   
            nurbs.derivative_on_curve<1>(0.2 * i, point);
            pointss.col(i) = point[0];
            ders.col(i) = point[1];
        }
        many_points.resize(3, 20);
        many_ders.resize(3, 20);
        for (int i = 0; i < 20; ++i)
        {
            Eigen::Vector<Eigen::Vector3d, 2> point;   
            nurbs.derivative_on_curve<1>(0.05 * i, point);
            many_points.col(i) = point[0];
            many_ders.col(i) = point[1];
        }
    }

    // nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    Eigen::Matrix<double, 3, Eigen::Dynamic> pointss;
    Eigen::Matrix<double, 3, Eigen::Dynamic> ders;
    Eigen::Matrix<double, 3, Eigen::Dynamic> many_points;
    Eigen::Matrix<double, 3, Eigen::Dynamic> many_ders;
};

class CreateNurbsCurve3 : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::Vector<double, 3> center{0, 0, 0};
        Eigen::Vector<double, 3> u_dir{1, 0, 0};
        Eigen::Vector<double, 3> v_dir{0, 1, 0};
        double radius = 10;
        double start_angles = 0;
        double end_angles = M_PI + 0.78;
        nurbs_curve<double, 3, true, -1, -1> nurbs;
        create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, nurbs);
        many_points.resize(3, 40);
        many_ders.resize(3, 40);
        for (int i = 0; i < 40; ++i)
        {
            Eigen::Vector<Eigen::Vector3d, 2> point;   
            nurbs.derivative_on_curve<1>(0.025 * i, point);
            many_points.col(i) = point[0];
            many_ders.col(i) = point[1];
        }
    }

    Eigen::Matrix<double, 3, Eigen::Dynamic> many_points;
    Eigen::Matrix<double, 3, Eigen::Dynamic> many_ders;
};

class CreateNurbsSurface : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::VectorX<double> u_knots_vector(7);
        u_knots_vector<< 0, 0, 0, 2.5, 5, 5, 5;
        Eigen::VectorX<double> v_knots_vector(7);
        v_knots_vector << 0, 0, 0, 1.5, 3, 3, 3;
        Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 4);
        points1 << 0, 0, 0, 0,
                2, 6, 2, 8,
                5, 4, 0, 2,
                1, 2, 1, 2;

        Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 4);
        points2 << 4, 12, 4, 8,
                6, 24, 10, 28,
                8, 12, 0, 0,
                2, 6,  2, 4;

        Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 4);
        points3 << 4, 8, 4, 12,
                2, 6, 4, 12,
                4, 4, 0, -3,
                1, 2, 1, 3;

        Eigen::Matrix<double, 4, Eigen::Dynamic> points4(4, 4);
        points4 << 4, 8, 4, 12,
                2, 6, 4, 12,
                4, 8, 4, 12,
                1, 2, 1, 3;

        Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(4);
        control_points(0) = points1;
        control_points(1) = points2;
        control_points(2) = points3;
        control_points(3) = points4;
        test_surface =  nurbs_surface<double, 3, -1, -1, -1, -1, true>(u_knots_vector, v_knots_vector, control_points);
        points.resize(5);
        ders.resize(5);
        for (int i = 0; i < 5; ++i)
        {
            points[i].resize(3, 5);
            ders[i].resize(3, 5);
            for (int j = 0; j < 5; ++j)
            {
                Eigen::Matrix<Eigen::Vector<double, 3>, 2, 2> ders_points;
                test_surface.derivative_on_surface<1>(i * 1.0, j * 0.6, ders_points);
                points[i].col(j) = ders_points(0, 0);
                ders[i].col(j) = ders_points(1, 0);
            }
        }
    }

    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface;
    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> points;
    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> ders;
};

class CreateNurbsCurve2 : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::Vector<double, 3> center{0, 0, 0};
        Eigen::Vector<double, 3> u_dir{1, 0, 0};
        Eigen::Vector<double, 3> v_dir{0, 1, 0};
        double radius = 10;
        double start_angles = 0;
        double end_angles = M_PI + 0.78;
        nurbs_curve<double, 3, true, -1, -1> nurbs;
        create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, nurbs);

        Eigen::Matrix<double, 3, Eigen::Dynamic> pointss(3, 5);
        for (int i = 0; i < 5; ++i)
        {
            Eigen::Vector3d point;   
            nurbs.point_on_curve(0.2 * i, point);
            pointss.col(i) = point;
        }
        global_curve_interpolate<double, 3, ENPARAMETERIEDTYPE::CHORD>(pointss, 3, m_nurbs);
        Eigen::VectorX<double> insert_knots(5);
        insert_knots << 0.2, 0.2 , 0.5, 0.7, 0.7;
        m_nurbs.refine_knots_vector(insert_knots);
        
    }
    
    nurbs_curve<double, 3, false, -1, -1> m_nurbs;
     

};

class CreateNurbsCurve4 : public testing::Test
{
protected:
    void SetUp() override
    {
        Eigen::Vector<double, 2> v1{0, 1};
        Eigen::Vector<double, 2> v2{1.0 / 2.0, 2};
        Eigen::Vector<double, 2> v3{1, 1};
        Eigen::Matrix<double, 2, Eigen::Dynamic> mat(2, 3);
        mat.col(0) = v1;
        mat.col(1) = v2;
        mat.col(2) = v3;
        Eigen::VectorX<double> knots_vector(6);
        knots_vector << 0, 0, 0, 1, 1, 1;
        m_nurbs1.set_control_points(mat);//  = nurbs_curve<double, 2, false, -1, -1>(knots_vector, mat);
        m_nurbs1.set_knots_vector(knots_vector);
        m_nurbs1.set_degree(2);

        Eigen::Vector<double, 3> wv1{0, 0, 1};
        Eigen::Vector<double, 3> wv2{0.3, 0.2, 0.5};
        Eigen::Vector<double, 3> wv3{1.4 / 2.0, 1, 2};
        Eigen::Vector<double, 3> wv4{2.6, 0.7, 3};
        Eigen::Vector<double, 3> wv5{1.5, 0.5, 0.5};
        Eigen::Vector<double, 3> wv6{3.4, 1.2, 1};
        Eigen::Matrix<double, 3, Eigen::Dynamic> mat2(3, 6);
        mat2.col(0) = wv1;
        mat2.col(1) = wv2;
        mat2.col(2) = wv3;
        mat2.col(3) = wv4;
        mat2.col(4) = wv5;
        mat2.col(5) = wv6;


        Eigen::VectorX<double> knots_vector2(10);
        knots_vector2 << 0, 0, 0, 0, 0.3, 0.7, 1, 1, 1, 1;


        m_nurbs2.set_knots_vector(knots_vector2);
        m_nurbs2.set_control_points(mat2);
        m_nurbs2.set_degree(3);
    }
    
    nurbs_curve<double, 2, false, -1, -1> m_nurbs1;
    nurbs_curve<double, 2, true, -1, -1> m_nurbs2;

};



TEST_F(CreateNurbsCurve2, RemoveKnots1)
{
    std::vector<double> params(9);
    for (int index = 1; index < 10; ++index)
        params[index - 1] = 0.1 * index;
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    std::vector<double> errors;
    m_nurbs.remove_knots_bound_curve(params, errors, new_nurbs, 2.0);
    std::vector<double> real_errors;
    for (int index = 1; index < 10; ++index)
    {
        Eigen::Vector3d p1, p2;
        m_nurbs.point_on_curve(0.1 * index, p1);
        new_nurbs.point_on_curve(0.1 * index, p2);
        double d = (p1 - p2).norm();
        real_errors.push_back(d);
        EXPECT_TRUE(d < errors[index - 1] + DEFAULT_ERROR);
    }

    ASSERT_TRUE(true == true);
}
TEST_F(CreateNurbsCurve2, RemoveKnots2)
{
    std::vector<double> params(5);
    for (int index = 0; index < 5; ++index)
        params[index] = 0.2 * index;
    Eigen::VectorX<double> knots(7);
    knots << 0, 0, 0, 0.522576, 1, 1, 1;
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 4);
    points << 0, 0, 0, 0,
              2, 2.99661, 2.45643, 3.11111,
              5 ,3.78986, 0.588487, 0.722222;
    nurbs_curve<double, 3, false, -1, -1> new_nurbs(knots, points);
    nurbs_curve<double, 3, false, -1, -1> new_nurbs2;
    std::vector<double> errors;
    std::vector<double> real_errors;
    new_nurbs.remove_knots_bound_curve(params, errors, new_nurbs2, 2.0);
    for (int index = 0; index < 5; ++index)
    {
        Eigen::Vector3d p1, p2;
        new_nurbs2.point_on_curve(0.2 * index, p1);
        new_nurbs.point_on_curve(0.2 * index, p2);
        double d = (p1 - p2).norm();
        real_errors.push_back(d);
        EXPECT_TRUE(d <= errors[index] + DEFAULT_ERROR);
    }

    ASSERT_TRUE(true == true);
}

TEST_F(CreateNurbsCurve, GlobalCurveApproximationErrBnd1)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs2;

    global_curve_approximation_err_bnd<double, 3, ENPARAMETERIEDTYPE::CHORD>(many_points, 5, new_nurbs2, 1);
}


TEST_F(CreateNurbsCurve, InterpolateWithEndsTangent)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    Eigen::Vector3d D0{-1, 9, 0};
    Eigen::Vector3d D1{-5, -5, 0};
    global_curve_interpolate_with_ends_tangent<double, 3, ENPARAMETERIEDTYPE::CHORD>(pointss, D0, D1, 3, new_nurbs);

    ASSERT_TRUE(true == true);
}

TEST_F(CreateNurbsCurve2, FitWithConic)
{
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 100);
    Eigen::Matrix<double, 3, Eigen::Dynamic> tangents(3, 100);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector<Eigen::Vector3d, 2> p;
    
        m_nurbs.derivative_on_curve<1>(i * 0.01, p);
        points.col(i) = p[0];
        tangents.col(i) = p[1];
    }

    nurbs_curve<double, 3, true, -1, -1> test_nurbs;
    fit_with_conic<double, 3>(points, tangents, test_nurbs, 0.5);
    curve_nearest<nurbs_curve<double, 3, true, -1, -1>> curve_neareast_tool;
    // curve_neareast_tool.set_dis_eps(0.1);
    curve_neareast_tool.init(&test_nurbs);
    double temp_u, temp_min_dis;
    Eigen::Vector3d temp_point;
    for (int i = 0; i < 100; ++i)
    {
        curve_neareast_tool.find_nearst_point_on_curve(points.col(i), temp_u, temp_point, temp_min_dis);
        // std::cout << i << std::endl;
        double d = (temp_point - points.col(i)).norm();
        EXPECT_NEAR(d, 0.0, 0.5);
    }
}

TEST_F(CreateNurbsCurve2, InterpolateWithEndsTangent3)
{
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 10);
    Eigen::Matrix<double, 3, Eigen::Dynamic> tangents(3, 10);
    for (int i = 0; i < 10; ++i)
    {
        Eigen::Vector<Eigen::Vector3d, 2> p;
    
        m_nurbs.derivative_on_curve<1>(i * 0.1, p);
        points.col(i) = p[0];
        tangents.col(i) = p[1];
    }
    nurbs_curve<double, 3, false, -1, -1> test_nurbs;
    fit_with_cubic<double, 3>(points, tangents, test_nurbs, 0.1);
    for (int i = 0; i < 10; ++i)
    {
        Eigen::Vector3d p;
        double u;
        test_nurbs.find_nearst_point_on_curve(points.col(i), u, p);
        double d = (p - points.col(i)).norm();
        EXPECT_NEAR(d, 0.0, 0.1);
    }
}

TEST_F(CreateNurbsCurve2, InterpolateWithEndsTangent4)
{
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 10);
    Eigen::Matrix<double, 3, Eigen::Dynamic> tangents(3, 10);
    double step = M_PI / 5;
    for (int index = 0; index < 10; ++index)
    {

        double x = 5.0 * std::cos(index * step);
        double y = 5.0 * std::sin(index * step);
        double z = 0.5 * (double)index;
        points(0, index) = x;
        points(1, index) = y;
        points(2, index) = z;
        tangents(0, index) = -y;
        tangents(1, index) = x;
        tangents(2, index) = 0.5;
    }
    nurbs_curve<double, 3, false, -1, -1> test_nurbs;
    fit_with_cubic<double, 3>(points, tangents, test_nurbs, 0.1);

    curve_nearest<nurbs_curve<double, 3, false, -1, -1>> curve_neareast_tool;
    curve_neareast_tool.init(&test_nurbs);
    double temp_u, temp_min_dis;
    Eigen::Vector3d temp_point;
    
    for (int i = 0; i < 10; ++i)
    {
        curve_neareast_tool.find_nearst_point_on_curve(points.col(i), temp_u, temp_point, temp_min_dis);
        double d = (temp_point - points.col(i)).norm();
        EXPECT_NEAR(d, 0.0, 0.1);
    }

}


TEST_F(CreateNurbsCurve, InterpolateWithEndsTangent2)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    Eigen::Vector3d D0{-1, 9, 0};
    Eigen::Vector3d D1{-5, -5, 0};
    global_3degree_curve_interpolate_with_ends_tangent<double, 3, ENPARAMETERIEDTYPE::CHORD>(pointss, D0, D1, new_nurbs);
    Eigen::Vector3d point;
    new_nurbs.point_on_curve(0.3, point);
    Eigen::Vector3d test_point(5.9121941939644582, 8.239425086235828, 0.0);
    double distance = (point - test_point).norm();
    EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);

    new_nurbs.point_on_curve(0.4, point);
    test_point = Eigen::Vector3d(3.1448187971936421, 9.6107198158545994, 0.0);
    distance = (point - test_point).norm();
    EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);

    new_nurbs.point_on_curve(0.5, point);
    test_point = Eigen::Vector3d(0.047821952388199605, 10.000425510364781, 0.0);
    distance = (point - test_point).norm();
    EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);
}

TEST_F(CreateNurbsCurve, InterpolateHermite)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    global_2or3degree_hermite_curve<double, 3, 3, ENPARAMETERIEDTYPE::CHORD>(pointss, ders, new_nurbs);
    Eigen::Vector3d point;

    for (int index = 0; index < 5; ++index)
    {
        Eigen::Vector<Eigen::Vector3d, 2> point;
        double u;
        Eigen::Vector3d p;
        ENUM_NURBS flag = new_nurbs.find_nearst_point_on_curve(pointss.col(index), u, p);
        ASSERT_EQ(ENUM_NURBS::NURBS_SUCCESS, flag);
        new_nurbs.derivative_on_curve<1>(u, point);

        double distance = (point[0] - pointss.col(index)).norm();
        EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);

        distance = (point[1] - ders.col(index)).norm();
        EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);
    }

}

TEST_F(CreateNurbsCurve, Local2degreeInterpolate)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    local_2degree_parabola_interpolate<double, 3, true>(pointss, new_nurbs);

    std::vector<Eigen::Vector3d> pointss1;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;   
        new_nurbs.point_on_curve(0.01 * i, point);
        pointss1.push_back(point);
    }
}

TEST_F(CreateNurbsCurve, LocalArcInterpolate)
{
    nurbs_curve<double, 3, true, -1, -1> new_nurbs;
    local_2degree_arc_interpolate<double, 3>(pointss, ders, new_nurbs);

}

TEST_F(CreateNurbsCurve, LocalArcInterpolate2)
{
    nurbs_curve<double, 3, true, -1, -1> new_nurbs;
    local_2degree_arc_interpolate<double, 3, true>(pointss, new_nurbs);
}

TEST_F(CreateNurbsCurve, Local3degreeInterpolate1)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    local_3degree_interpolate<double, 3>(pointss, ders, new_nurbs);

    const Eigen::VectorX<double> knots_vector = new_nurbs.get_knots_vector();
    int knots_count = knots_vector.size();
    std::vector<double> params;
    params.push_back(knots_vector[0]);
    for (int index = 1; index < knots_count; ++index)
    {
        if (knots_vector[index] != params.back())
        {
            params.push_back(knots_vector[index]);
        }
    }
    int params_count = params.size();
    std::vector<double> mid_params(params_count - 1);
    for (int index = 1; index < params_count; ++index)
    {
        mid_params[index - 1] = (params[index] + params[index - 1]) / 2.0;
    }
    Eigen::Vector<Eigen::Vector<double, 3>, 2> n_ders;
    new_nurbs.derivative_on_curve<1>(params[0], n_ders);
    double tangent_norm = n_ders[1].norm();
    for (int index = 1; index < params_count; ++index)
    {
        new_nurbs.derivative_on_curve<1>(params[index], n_ders);
        double t1 = n_ders[1].norm();
        EXPECT_NEAR(t1, tangent_norm, DEFAULT_ERROR);

        Eigen::Vector<Eigen::Vector<double, 3>, 2> n_ders2;
        new_nurbs.derivative_on_curve<1, ENUM_LIMITDIRECTION::LEFT>(params[index], n_ders2);
        double t = (n_ders2[1] - n_ders[1]).norm();
        EXPECT_NEAR(t, 0.0, DEFAULT_ERROR);

        new_nurbs.derivative_on_curve<1>(mid_params[index - 1], n_ders);
        double t2 = n_ders[1].norm();
        EXPECT_NEAR(t2, tangent_norm, DEFAULT_ERROR);
        
    }
}

TEST_F(CreateNurbsCurve, Local3degreeInterpolate2)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    local_3degree_interpolate<double, 3>(pointss, new_nurbs);

    const Eigen::VectorX<double> knots_vector = new_nurbs.get_knots_vector();
    int knots_count = knots_vector.size();
    std::vector<double> params;
    params.push_back(knots_vector[0]);
    for (int index = 1; index < knots_count; ++index)
    {
        if (knots_vector[index] != params.back())
        {
            params.push_back(knots_vector[index]);
        }
    }
    int params_count = params.size();
    std::vector<double> mid_params(params_count - 1);
    for (int index = 1; index < params_count; ++index)
    {
        mid_params[index - 1] = (params[index] + params[index - 1]) / 2.0;
    }
    Eigen::Vector<Eigen::Vector<double, 3>, 2> n_ders;
    new_nurbs.derivative_on_curve<1>(params[0], n_ders);
    double tangent_norm = n_ders[1].norm();
    for (int index = 1; index < params_count; ++index)
    {
        new_nurbs.derivative_on_curve<1>(params[index], n_ders);
        double t1 = n_ders[1].norm();
        EXPECT_NEAR(t1, tangent_norm, DEFAULT_ERROR);

        Eigen::Vector<Eigen::Vector<double, 3>, 2> n_ders2;
        new_nurbs.derivative_on_curve<1, ENUM_LIMITDIRECTION::LEFT>(params[index], n_ders2);
        double t = (n_ders2[1] - n_ders[1]).norm();
        EXPECT_NEAR(t, 0.0, DEFAULT_ERROR);

        new_nurbs.derivative_on_curve<1>(mid_params[index - 1], n_ders);
        double t2 = n_ders[1].norm();
        EXPECT_NEAR(t2, tangent_norm, DEFAULT_ERROR);
        
    }
}

TEST_F(CreateNurbsCurve, globalLeastSquaresCurveApproximation1)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    global_least_squares_curve_approximation<double, 3, ENPARAMETERIEDTYPE::CHORD>(pointss, 3, 4, new_nurbs);

    double u = 0.2;
    Eigen::Vector3d pos;
    new_nurbs.point_on_curve(u, pos);
    Eigen::Vector3d expect_pos(8.1336662545570189, 6.2001553478557714, 0);
    double distance = (pos - expect_pos).norm();
    EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);

    u += 0.2;
    new_nurbs.point_on_curve(u, pos);
    Eigen::Vector3d expect_pos2(3.0951955345688122, 9.304855647197515, 0);
    distance = (pos - expect_pos2).norm();
    EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);
}

TEST_F(CreateNurbsCurve, globalLeastSquaresCurveApproximation2)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;

    std::vector<double> w{-1, 100, 100, -1, 1, 1, -3, 4, 2, 1, 3, -2, 1, 1, 1, 1, 1, 2, 3, -3};
    std::vector<int> index_set;
    for (int index = 0; index < 20; ++index)
        index_set.push_back(index);
    std::vector<double> wd{-1, 2, 4, 1, 1, 1, 3, 4, 2, 1, 3, 2, 1, 1, 1, 1, 1, 2, 3, -3};
    ENUM_NURBS flag = global_wc_least_squares_curve_approximation<double, 3, ENPARAMETERIEDTYPE::CHORD>(many_points, w, many_ders, index_set, wd, 3, 9, new_nurbs);
    EXPECT_EQ(flag, ENUM_NURBS::NURBS_SUCCESS);
    Eigen::Vector2<Eigen::Vector3d> point;
    new_nurbs.derivative_on_curve<1>(0.0, point);
    double d1 = (point[0] - many_points.col(0)).norm();
    double d2 = (point[1] - many_ders.col(0)).norm();
    EXPECT_NEAR(d1, 0.0, DEFAULT_ERROR);
    EXPECT_NEAR(d2, 0.0, DEFAULT_ERROR);

    new_nurbs.derivative_on_curve<1>(0.3, point);
    Eigen::Vector3d t1(4.344122821982781, 8.9977669582779232, 0.0), t2(-33.048646954648341, 16.473034652625508, 0.0);
    d1 = (point[0] - t1).norm();
    d2 = (point[1] - t2).norm();
    EXPECT_NEAR(d1, 0.0, DEFAULT_ERROR);
    EXPECT_NEAR(d2, 0.0, DEFAULT_ERROR);

    new_nurbs.derivative_on_curve<1>(1.0, point);
    d1 = (point[0] - many_points.col(19)).norm();
    d2 = (point[1] - many_ders.col(19)).norm();
    EXPECT_NEAR(d1, 0.0, DEFAULT_ERROR);
    EXPECT_NEAR(d2, 0.0, DEFAULT_ERROR);

}

TEST_F(CreateNurbsSurface,SurfaceSimplify)
{
    nurbs_surface<double, 3, -1, -1, -1, -1, false> new_nurbs;

    global_surface_approximation<double, 3, ENPARAMETERIEDTYPE::CHORD>(points, 3, 2, 4, 4, new_nurbs);
    new_nurbs.surface_knots_insert(0.2, 1, ENUM_DIRECTION::U_DIRECTION);
    new_nurbs.surface_knots_insert(0.8, 1, ENUM_DIRECTION::U_DIRECTION);
    new_nurbs.surface_knots_insert(0.2, 2, ENUM_DIRECTION::V_DIRECTION);
    new_nurbs.surface_knots_insert(0.7, 2, ENUM_DIRECTION::V_DIRECTION);
    nurbs_surface<double, 3, -1, -1, -1, -1, false> new_nurbs2;
    std::vector<std::array<double, 2>> params;
    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> fit_params(5);
    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> fit_params2(5);
    for (int i = 0; i < 5; ++i)
    {
        fit_params[i].resize(3, 5);
        fit_params2[i].resize(3, 5);
    }
    for (int index = 0; index < 5; ++index)
    {
        for (int v_index = 0; v_index < 5; ++v_index)
        {
            std::array<double, 2> temp{0.2 * index, 0.2 * v_index};
            params.push_back(std::move(temp));
        }

    }
    std::vector<double> errors;
    new_nurbs.remove_knots_bound_surface(params, errors, new_nurbs2, 2.0);

    std::vector<double> real_error;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            Eigen::Vector3d p1, p2;
            new_nurbs2.point_on_surface(0.2 * i, 0.2 * j, p1);
            new_nurbs.point_on_surface(0.2 * i, 0.2 * j, p2);
            fit_params[i].col(j) = p1;
            fit_params2[i].col(j) = p2;
            double error = (p1 - p2).norm();
            real_error.push_back(error);
            EXPECT_TRUE(errors[i * 5 + j] + DEFAULT_ERROR > error);
        }
    }
}


TEST_F(CreateNurbsSurface, globalLeastSquaresSurfaceApproximation1)
{
    nurbs_surface<double, 3, -1, -1, -1, -1, false> new_nurbs;

    ENUM_NURBS flag = global_surface_approximation<double, 3, ENPARAMETERIEDTYPE::CHORD>(points, 3, 2, 4, 4, new_nurbs);

    EXPECT_EQ(flag, ENUM_NURBS::NURBS_SUCCESS);
    Eigen::Vector3d point;
    new_nurbs.point_on_surface(0.0, 0.0, point);
    double d = (point - points[0].col(0)).norm();
    EXPECT_NEAR(d, 0.0, DEFAULT_ERROR);

    new_nurbs.point_on_surface(1.0, 0.0, point);
    d = (point - points[0].col(4)).norm();
    EXPECT_NEAR(d, 0.0, DEFAULT_ERROR);

    new_nurbs.point_on_surface(0.0, 1.0, point);
    d = (point - points[4].col(0)).norm();
    EXPECT_NEAR(d, 0.0, DEFAULT_ERROR);
    
    new_nurbs.point_on_surface(1.0, 1.0, point);
    d = (point - points[4].col(4)).norm();
    EXPECT_NEAR(d, 0.0, DEFAULT_ERROR);

    Eigen::Vector3d test(0.64504075572135855, 3.484104533792812, 1.4768769670319279);
    
    new_nurbs.point_on_surface(0.2, 0.7, point);
    d = (point - test).norm();
    EXPECT_NEAR(d, 0.0, DEFAULT_ERROR);
}




TEST_F(CreateNurbsSurface, LoacalInterpolateSurface)
{
    nurbs_surface<double, 3, -1, -1, -1, -1, false> new_nurbs;

    local_bi3degree_interpolate<double, 3, ENPARAMETERIEDTYPE::CHORD>(points, new_nurbs);
    surface_nearest<nurbs_surface<double, 3, -1, -1, -1, -1, false>> surface_nearest_tool;
    surface_nearest_tool.init(&new_nurbs);
    double temp_min_dis;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            Eigen::Vector3d point;
            double u, v;
            Eigen::Vector3d p;
            ENUM_NURBS flag = surface_nearest_tool.find_nearst_point_on_surface(points[i].col(j), u, v, p, temp_min_dis);
            ASSERT_EQ(ENUM_NURBS::NURBS_SUCCESS, flag);

            double distance = (p - points[i].col(j)).norm();
            EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);
        }
    }

}



TEST_F(CreateNurbsSurface, InterpolateSurface)
{
    nurbs_surface<double, 3, -1, -1, -1, -1, false> new_nurbs;

    global_surface_interpolate<double, 3, ENPARAMETERIEDTYPE::CHORD>(points, 3, 4, new_nurbs);
    surface_nearest<nurbs_surface<double, 3, -1, -1, -1, -1, false>> surface_nearest_tool;
    
    surface_nearest_tool.init(&new_nurbs);
    double temp_min_dis;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            Eigen::Vector3d point;
            double u, v;
            Eigen::Vector3d p;
            ENUM_NURBS flag = surface_nearest_tool.find_nearst_point_on_surface(points[i].col(j), u, v, p, temp_min_dis);
            EXPECT_EQ(ENUM_NURBS::NURBS_SUCCESS, flag);

            double distance = (p - points[i].col(j)).norm();
            EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);
        }
    }
}

TEST_F(CreateNurbsCurve4, SwungSurface)
{

    nurbs_surface<double, 3, -1, -1, -1, -1, true> sur;
    swung_surface<double, true, false>(1.0, m_nurbs2, m_nurbs1, sur);
    nurbs_curve<double, 3, false, -1, -1> new_nurbs1;
    nurbs_curve<double, 3, true, -1, -1> new_nurbs2;
    m_nurbs1.dimension_elevate(2, new_nurbs1);
    m_nurbs2.dimension_elevate(1, new_nurbs2);
}

TEST_F(CreateNurbsCurve2, SkinSurface1)
{
    std::vector<nurbs_curve<double, 3, false, -1, -1> *> nurbs_curves;
    nurbs_curves.push_back(&m_nurbs);
    nurbs_curve<double, 3, false, -1, -1> nurbs2;
    m_nurbs.move(Eigen::Vector3d(10, 20, 10), nurbs2);
    nurbs_curve<double, 3, false, -1, -1> nurbs3;
    nurbs2.move(Eigen::Vector3d(0, 10, 10), nurbs3);
    nurbs_curves.push_back(&nurbs2);
    nurbs_curves.push_back(&nurbs3);


    nurbs_surface<double, 3, -1 ,-1, -1, -1, false> skin;
    std::vector<double> v_params{0, 0.5, 1};
    skin_surface<double, 3, false>(2, nurbs_curves, skin);
    Eigen::Vector3d pos;
    Eigen::Vector3d test_pos(17.011133746658562, 24.452491384379741, 8.3548178734335181);
    skin.point_on_surface(0.23, 0.56, pos);
    double d = (pos - test_pos).norm();
    EXPECT_NEAR(d, 0.0, DEFAULT_ERROR);

    skin_surface<double, 3, false>(2, v_params, nurbs_curves, skin);
    test_pos = Eigen::Vector3d(18.052542683375325, 28.339082447662982, 11.199999999999999);
    skin.point_on_surface(0.23, 0.56, pos);
    d = (pos - test_pos).norm();
    EXPECT_NEAR(d, 0.0, DEFAULT_ERROR);
}

TEST_F(CreateNurbsCurve2, SkinSurface2)
{
    std::vector<nurbs_curve<double, 3, false, -1, -1> *> nurbs_curves;
    nurbs_curves.push_back(&m_nurbs);
    nurbs_curve<double, 3, false, -1, -1> nurbs2;
    m_nurbs.move(Eigen::Vector3d(10, 20, 10), nurbs2);
    nurbs_curve<double, 3, false, -1, -1> nurbs3;
    nurbs2.move(Eigen::Vector3d(0, 10, 10), nurbs3);
    nurbs_curves.push_back(&nurbs2);
    nurbs_curves.push_back(&nurbs3);
    
    nurbs_curve<double, 2, false, -1, -1> m_nurbs1;
    Eigen::Vector<double, 2> v1{0, 10};
    Eigen::Vector<double, 2> v2{10.0 / 2.0, 20};
    Eigen::Vector<double, 2> v3{10, 25};
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat(2, 3);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    Eigen::VectorX<double> knots_vector(6);
    knots_vector << 0, 0, 0, 1, 1, 1;
    m_nurbs1.set_control_points(mat);//  = nurbs_curve<double, 2, false, -1, -1>(knots_vector, mat);
    m_nurbs1.set_knots_vector(knots_vector);
    m_nurbs1.set_degree(2);
    nurbs_curve<double, 3, false, -1, -1> new_nurbs1;
    m_nurbs1.dimension_elevate(0, new_nurbs1);

    nurbs_surface<double, 3, -1 ,-1, -1, -1, false> skin;
    std::vector<double> v_params{0, 0.5, 1};
    skin_surface<double, 3, false, 2>(nurbs_curves, new_nurbs1, v_params, skin);
    Eigen::Matrix<Eigen::Vector3d, 2, 2> tangents;
    skin.derivative_on_surface<1>(0, 0.63397459621556151, tangents);
    Eigen::Vector2<Eigen::Vector3d> tangents2;
    new_nurbs1.derivative_on_curve<1>(0.5, tangents2);
    double cp = tangents2[1].cross(tangents(0, 1)).norm();
    EXPECT_NEAR(cp, 0.0, DEFAULT_ERROR);
}

TEST_F(CreateNurbsCurve2, sweep)
{
    Eigen::Vector<double, 3> center{0, 0, 0};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 0, 1};
    double radius = 10;
    double start_angles = 0;
    double end_angles = M_PI * 2.0;
    nurbs_curve<double, 3, true, -1, -1> spine;
    create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, spine);

    Eigen::Vector<double, 3> x_dir{0, 1, 0};
    radius = 2;
    nurbs_curve<double, 3, true, -1, -1> profile;
    create_nurbs_circle(center, x_dir, v_dir, radius, start_angles, end_angles, profile);

    nurbs_surface<double, 3, -1, -1, -1, -1, true> surf;
    sweep_surface<double, true, true>(profile, spine, 10, surf);

    Eigen::Vector3d pos;
    surf.point_on_surface(0, 0.63397459621556151, pos);
    Eigen::Vector3d test_pos(-5.3106937415249327, 0, -5.9830203061437084);
    double cp = (pos - test_pos).norm();
    EXPECT_NEAR(cp, 0.0, DEFAULT_ERROR);
}

TEST_F(CreateNurbsCurve2, sweep2)
{
    Eigen::Vector<double, 3> center{0, 0, 0};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 0, 1};
    double radius = 10;
    double start_angles = 0;
    double end_angles = M_PI * 2.0;
    nurbs_curve<double, 3, true, -1, -1> spine;
    create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, spine);

    Eigen::Vector<double, 3> x_dir{0, 1, 0};
    radius = 2;
    nurbs_curve<double, 3, true, -1, -1> profile;
    create_nurbs_circle(center, x_dir, v_dir, radius, start_angles, end_angles, profile);

    nurbs_surface<double, 3, -1, -1, -1, -1, true> surf;
    sweep_surface<double, true, true>(profile, spine, 6, 10, surf);
    Eigen::Vector3d pos;
    surf.point_on_surface(0, 0.63397459621556151, pos);
    Eigen::Vector3d test_pos(-5.3202257954846646, 0, -5.9357274224120795);
    double cp = (pos - test_pos).norm();
    EXPECT_NEAR(cp, 0.0, DEFAULT_ERROR);
}

TEST_F(CreateNurbsCurve2, GordonSurface2)
{    
    std::vector<nurbs_curve<double, 3, false, -1, -1> *> nurbs_curves;
    nurbs_curves.push_back(&m_nurbs);
    nurbs_curve<double, 3, false, -1, -1> nurbs2;
    m_nurbs.move(Eigen::Vector3d(10, 20, 10), nurbs2);
    nurbs_curve<double, 3, false, -1, -1> nurbs3;
    nurbs2.move(Eigen::Vector3d(0, 10, 10), nurbs3);
    nurbs_curves.push_back(&nurbs2);
    nurbs_curves.push_back(&nurbs3);
    double step = 1.0 / 3.0;
    std::vector<nurbs_curve<double, 3, false, -1, -1> *> nurbs_curves2(4, nullptr);
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 3);
    std::vector<double> u_params{0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0};
    std::vector<double> v_params{0.0, 0.5, 1.0};
    for (int i = 0; i < 4; ++i)
    {
        double param = i * step;
        if (i == 3)
            param = 1.0;
        Eigen::Vector3d p;
        for (int j = 0; j < 3; ++j)
        {
            nurbs_curves[j]->point_on_curve(param, p);
            points.col(j) = p;
        }
        nurbs_curve<double, 3, false, -1, -1> *temp = new nurbs_curve<double, 3, false, -1, -1>();
        global_curve_interpolate(points, 2, v_params, *temp);   
        nurbs_curves2[i] = temp;
    }
    nurbs_surface<double, 3, -1, -1, -1 ,-1, false> gsurf;
    gordon_surface(nurbs_curves, nurbs_curves2, u_params, v_params, 1, 2, gsurf);
    Eigen::Vector3d pos;
    gsurf.point_on_surface(0.3, 0.63397459621556151, pos);
    Eigen::Vector3d test_pos(16.858787768994848, 31.717702537988608, 12.679491924311229);
    double cp = (pos - test_pos).norm();
    EXPECT_NEAR(cp, 0.0, DEFAULT_ERROR);
}


TEST_F(CreateNurbsCurve2, Coons1)
{    
    
    Eigen::Vector<double, 3> center{0, 0, 0};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 1, 0};
    double radius = 10;
    double start_angles = 0;
    double end_angles = M_PI + 0.78;
    nurbs_curve<double, 3, true, -1, -1> nurbs;
    create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, nurbs);

    std::vector<nurbs_curve<double, 3, true, -1, -1> *> nurbs_curves(3, nullptr);
    nurbs_curves[0] = &nurbs;
    nurbs_curve<double, 3, true, -1, -1> nurbs2;
    nurbs.move(Eigen::Vector3d(10, 20, 10), nurbs2);
    nurbs_curves[1] = &nurbs2;
    nurbs_curve<double, 3, true, -1, -1> nurbs3;
    nurbs2.move(Eigen::Vector3d(0, 10, 0), nurbs3);
    nurbs_curves[2] = &nurbs3;
    std::array<nurbs_curve<double, 3, true, -1, -1> *, 2> nurbs_curves2;
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 3);
    std::vector<double> u_params{0.0, 1.0};
    std::vector<double> v_params{0.0, 1.0};
    std::vector<double> params{0, 0.5, 1.0};
    for (int i = 0; i < 2; ++i)
    {
        Eigen::Vector3d p;
        for (int j = 0; j < 3; ++j)
        {
            nurbs_curves[j]->point_on_curve(v_params[i], p);
            points.col(j) = p;
        }
        
        nurbs_curve<double, 3, false, -1, -1> temp;
        global_curve_interpolate(points, 2, params, temp);
        nurbs_curve<double, 3, true, -1, -1> *temp1 = new nurbs_curve<double, 3, true, -1, -1>();
        temp.to_rational_nurbs(*temp1);
        nurbs_curves2[i] = temp1;
    }
    nurbs_surface<double, 3, -1, -1, -1 ,-1, true> csurf;
    std::array<nurbs_curve<double, 3, true, -1, -1>*, 2> new_nurbs{nurbs_curves[0], nurbs_curves[2]};
    coons_surface(new_nurbs, nurbs_curves2, csurf);
    delete nurbs_curves2[0];
    delete nurbs_curves2[1];
    Eigen::Vector3d pos;
    csurf.point_on_surface(0.3, 0.63397459621556151, pos);
    Eigen::Vector3d test_pos(14.939318509409294, 33.097216638049076, 11.159567219564558);
    double cp = (pos - test_pos).norm();
    EXPECT_NEAR(cp, 0.0, DEFAULT_ERROR);

    Eigen::MatrixX<nurbs_surface<double, 3, -1, -1, -1, -1, true> *> segment;
}

TEST_F(CreateNurbsCurve2, FindNearstPoint1)
{    
    _set_error_mode(_OUT_TO_MSGBOX);
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1, 3, 3, 3;
    Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 4);
    points1 << 0, 0, 0, 0,
               2, 6, 2, 8,
               5, 4, 0, 2,
               1, 2, 1, 2;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 4);
    points2 << 4, 12, 4, 8,
               6, 24, 10, 28,
               8, 12, 0, 0,
               2, 6,  2, 4;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 4);
    points3 << 4, 8, 4, 12,
               2, 6, 4, 12,
               4, 4, 0, -3,
               1, 2, 1, 3;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points4(4, 4);
    points4 << 4, 8, 4, 12,
               2, 6, 4, 12,
               4, 8, 4, 12,
               1, 2, 1, 3;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface(u_knots_vector, v_knots_vector, control_points);


    std::vector<Eigen::Vector3d> points22;
    Eigen::Vector<Eigen::Vector3d, 2> normals;
    normals[0] = Eigen::Vector3d(0, 0, 1);
    normals[1] = Eigen::Vector3d(0, 0, -1);
    for (int i = 30; i < 35; ++i)
    {
        for (int j = 50; j < 55; ++j)
        {
            Eigen::Vector3d point;
            test_surface.point_on_surface(0.05 * i, 0.03 * j, point);
        }

    }

    surface_nearest<nurbs_surface<double, 3, -1, -1, -1, -1, true>> surface_nearest_tool;
    surface_nearest_tool.init(&test_surface);
    double temp_u, temp_v, temp_min_dis;
    Eigen::Vector3d temp_min_point;
    
    for (int i = 30; i < 35; ++i)
    {
        for (int j = 50; j < 55; ++j)
        {
            Eigen::Vector3d point;
            test_surface.point_on_surface(0.05 * i, 0.03 * j, point);
            Eigen::Vector3d p = point + 1 * normals[i % 2], nearst_point;
            double u_p, v_p;
            surface_nearest_tool.find_nearst_point_on_surface(p, temp_u, temp_v, temp_min_point, temp_min_dis);
            test_surface.find_nearst_point_on_surface(p, u_p, v_p, nearst_point);
            double d1 = (nearst_point - p).norm();
            double d2 = (temp_min_point - p).norm();
            double dis = d1 - d2;
            // std::cout << dis << std::endl;
            EXPECT_NEAR(dis, 0.0, 1e-4);
            // ASSERT_TRUE(std::abs(dis) < 1e-4);
            
        }
    } 
}

TEST_F(CreateNurbsCurve2, DiscSurface)
{    
    
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1, 3, 3, 3;
    Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 4);
    points1 << 0, 0, 0, 0,
               2, 6, 2, 8,
               5, 4, 0, 2,
               1, 2, 1, 2;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 4);
    points2 << 4, 12, 4, 8,
               6, 24, 10, 28,
               8, 12, 0, 0,
               2, 6,  2, 4;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 4);
    points3 << 4, 8, 4, 12,
               2, 6, 4, 12,
               4, 4, 0, -3,
               1, 2, 1, 3;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points4(4, 4);
    points4 << 4, 8, 4, 12,
               2, 6, 4, 12,
               4, 8, 4, 12,
               1, 2, 1, 3;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface(u_knots_vector, v_knots_vector, control_points);


    surface_mesh_helper<nurbs_surface<double, 3, -1, -1, -1, -1, true>> mh;
    disc_surface(&test_surface, mh, TDEFAULT_ERROR<double>::value, 0.5, 0.1, 1.0);
    mesh<2> surface_mesh;
    mesh_help_to_mesh(mh, surface_mesh);
}


TEST_F(CreateNurbsCurve2, FindNearstPoint2)
{    
    
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, Eigen::Dynamic> mat(3, 4);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::VectorX<double> knots_vector(7);
    knots_vector << 0, 0, 0, 0.4, 1, 1, 1;
    nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);
    Eigen::Vector<Eigen::Vector3d, 2> normals;
    normals[0] = Eigen::Vector3d(1, -1, 0);
    normals[1] = Eigen::Vector3d(-1, 1, 0);
    
    curve_nearest<nurbs_curve<double, 3, false, -1, -1>> curve_neareast_tool;
    curve_neareast_tool.set_dis_eps(0.1);
    curve_neareast_tool.init(&curve1);
    double temp_u, temp_min_dis;
    Eigen::Vector3d temp_point;
    
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector3d p = point + 1 * normals[i % 2], nearst_point;
        curve_neareast_tool.find_nearst_point_on_curve(p, temp_u, temp_point, temp_min_dis);
        double u_p;
        curve1.find_nearst_point_on_curve(p, u_p, nearst_point);
        double d1 = (temp_point - nearst_point).norm();
        EXPECT_NEAR(d1, 0.0, 1e-4);
    }

}


TEST_F(CreateNurbsCurve2, DiscCurve)
{    
    std::vector<nurbs_curve<double, 3, false, -1, -1>*> curves;
    Eigen::VectorX<double> insert_knots(1);
    insert_knots << 0.2;
    m_nurbs.refine_knots_vector(insert_knots);
    m_nurbs.decompose_to_nurbs(curves);
    curve_mesh_helper<nurbs_curve<double, 3, false, -1, -1>> mh;
    disc_curve(&m_nurbs, mh, TDEFAULT_ERROR<double>::value, 20.0, 0.1, 0.1);
    ASSERT_TRUE(true);
}

TEST(CONVEXHELL, test1)
{
    std::vector<Eigen::Vector2<double>> test{ {0, 0}, { 1, 4}, {3, 3}, { 3, 1}, {5, 5}, { 5, 2 }, {7, 0}, {9, 6} };
    std::vector<Eigen::Vector2<double>> result = graham_scan(test);
    for (auto elem : result)
    {
        std::cout << "[" << elem[0] << " " << elem[1] << "], ";
    }
    std::cout << std::endl;
} 

TEST(BEZIER_INT, test1)
{
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat1(2, 2);
    mat1 << -1, 1,
            -1, 1;
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat2(2, 2);
    mat2 << -0.75, 1,
            0, 0;

    bezier_curve<double, 2, false, -1> b1(mat1);
    bezier_curve<double, 2, false, -1> b2(mat2);
    std::vector<Eigen::Vector<double, 2>> int_params = beziers_int(b1, b2);
    for (auto param : int_params)
    {
        std::cout << "[" << param[0] << " " << param[1] << "], ";
    }
    std::cout << std::endl;
}  


TEST(BEZIER_INT, test2)
{
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

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
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

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    control_points2(3) = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false> *test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    // save_obj2(*test_surface, "test_surface.obj");
    // save_obj2(*test_surface2, "test_surface2.obj");
    _set_error_mode(_OUT_TO_MSGBOX);
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface, test_surface2);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest2.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    read_chat_points(test_data, src_path + "BEZIER_INTtest2.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
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

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
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

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    control_points2(3) = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points2);
    // save_obj2(*test_surface, "test_surface.obj");
    // save_obj2(*test_surface2, "test_surface2.obj");
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
    read_chat_points(test_data, src_path + "BEZIER_INTtest3.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
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

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
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

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    control_points2(3) = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    _set_error_mode(_OUT_TO_MSGBOX);
    
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface2, test_surface);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    //save_chat_points(intersection, "BEZIER_INTtest4.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    read_chat_points(test_data, src_path + "BEZIER_INTtest4.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}

TEST(BEZIER_NORMAL, test_normal)
{
    _set_error_mode(_OUT_TO_MSGBOX);
    Eigen::VectorX<double> v_knots_vector2(8);
    v_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;

    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector << 0, 0, 0, 1, 1, 1;

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

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);

    nurbs_surface<double, 3, -1, -1, -1, -1, false> test_surface_normal;
    test_surface->eval_normal_surface(test_surface_normal);
    Eigen::MatrixX<Eigen::Vector<double, 3>> ders;
    for (int i = 0; i <= 10; ++i)
    {
        for (int j = 0; j <= 10; ++j)
        {
            double u = 0.1 * i;
            double v = 0.1 * j;
            test_surface->derivative_on_surface(1, u, v, ders);
            Eigen::Vector3d normal = ders(1, 0).cross(ders(0, 1));
            Eigen::Vector3d point;
            test_surface_normal.point_on_surface(u, v, point);
            double len1 = normal.norm();
            double len2 = point.norm();
            normal.normalize();
            point.normalize();
            double dot_value = point.dot(normal);
            EXPECT_NEAR(dot_value, 1.0, 1e-4);
            EXPECT_NEAR(len1, len2, 1e-4);
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
    points1 <<  0, 10, 20,
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
    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);

    Eigen::Matrix<double, 3, Eigen::Dynamic> points11(3, 4);
    points11 << 0, 2, 4, 6,
                0, 10, 10, 0,
                -200, -200, -200, -200;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points12(3, 4);
    points12 << 0, 2, 4, 6,
                13, 17,17, 13,
                5, 10, 15, 20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points13(3, 4);
    points13 << 0, 2, 4, 6,
                17, 13, 13, 17,
                5, 10, 15, 20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points14(3, 4);
    points14 << 0, 2, 4, 6,
                15, 10, 10, 15,
                -200, -200, -200, -200;

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    control_points2(3) = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    int i = 0;
    _set_error_mode(_OUT_TO_MSGBOX);
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface, test_surface2);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest5.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    read_chat_points(test_data, src_path + "BEZIER_INTtest5.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
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
    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
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

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    control_points2(3) = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    _set_error_mode(_OUT_TO_MSGBOX);
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface2, test_surface);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest6.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    read_chat_points(test_data, src_path + "BEZIER_INTtest6.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}


TEST(BEZIER_INT, test7)
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
        0, -5, -10,
        1.5625, 1.5625, 1.5625;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points2(3, 3);
    points2 << 0, 10, 20,
        10, -5, 0,
        1.5625, 1.5625, 1.5625;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points3(3, 3);
    points3 << 0, 10, 20,
        20, 15, 10,
        1.5625, 1.5625, 1.5625;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points4(3, 3);
    points4 << 0, 10, 20,
        30, 25, 20,
        1.5625, 1.5625, 1.5625;
    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);

    Eigen::Matrix<double, 3, Eigen::Dynamic> points11(3, 4);
    points11 << 0, 4, 8, 12,
        5, 10, 10, 5,
        -20, -20, -20, -20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points12(3, 4);
    points12 << 0, 0, 0, 12,
        8, 12, 12, 8,
        20, 5, 5, 20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points13(3, 4);
    points13 << 0, 12, 12, 12,
        12, 8, 8, 12,
        20, 5, 5, 20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points14(3, 4);
    points14 << 0, 4, 8, 12,
        15, 10, 10, 15,
        -20, -20, -20, -20;

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    control_points2(3) = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    // save_obj2(*test_surface, "test_surface.obj");
    // save_obj2(*test_surface2, "test_surface2.obj");
    _set_error_mode(_OUT_TO_MSGBOX);
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface, test_surface2);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest7.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    read_chat_points(test_data, src_path + "BEZIER_INTtest7.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }

}



TEST(BEZIER_INT, test8)
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
        9, 9, 9;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points2(3, 3);
    points2 << 0, 10, 20,
        10, 10, 10,
        9, 9, 9;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points3(3, 3);
    points3 << 0, 10, 20,
        20, 20, 20,
        9, 9, 9;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points4(3, 3);
    points4 << 0, 10, 20,
        30, 30, 30,
        9, 9, 9;
    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);

    Eigen::Matrix<double, 3, Eigen::Dynamic> points11(3, 4);
    points11 << 0, 4, 8, 12,
        5, 10, 10, 5,
        -20, -20, -20, -20;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points12(3, 4);
    points12 << 0, 4, 8, 12,
        5, 10, 10, 5,
        20, 10, 30, 10;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points13(3, 4);
    points13 << 0, 4, 8, 12,
        15, 10, 10, 15,
        20, 10, 30, 10;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points14(3, 4);
    points14 << 0, 4, 8, 12,
        15, 10, 10, 15,
        -20, -20, -20, -20;

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    control_points2(3) = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    _set_error_mode(_OUT_TO_MSGBOX);
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface, test_surface2);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest8.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    read_chat_points(test_data, src_path + "BEZIER_INTtest8.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}

TEST(BEZIER_INT, test9)
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
    points1 <<  -10, 0, 10,
                -30, -30, -30,
                3, -3, 3;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points2(3, 3);
    points2 <<  -10, 0, 10,
                -10, -10, -10,
                3, -3, 3;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points3(3, 3);
    points3 << -10, 0, 10,
                10, 10, 10,
                3, -3, 3;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points4(3, 3);
    points4 << -10, 0, 10,
                30, 30, 15,
                3, -3, 3;

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector2, control_points);
    test_surface->reverse_uv();

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

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(4);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    control_points2(3) = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector2, v_knots_vector2, control_points2);
    _set_error_mode(_OUT_TO_MSGBOX);
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface, test_surface2);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest9.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    read_chat_points(test_data, src_path + "BEZIER_INTtest9.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}

TEST(BEZIER_INT, test10)
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
    points1 <<  1 / 7.0, 0, 3 / 5.0,
               3.0 / 5, 1.0 / 5, 3.0 / 4,
               1, 0, 7.0 / 10;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points2(3, 3);
    points2 << 3.0 / 8, 4.0 / 9, 2.0 / 3,
        2.0 / 3, 3.0 / 4, 1 / 3.0,
        6.0 / 7, 3.0 / 8, 5.0 / 7;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points3(3, 3);
    points3 <<  1 / 5.0,  6 / 7.0 ,4 / 7.0,
        4 / 3.0,  7 / 8.0 ,3 / 4.0,
        7 / 8.0,  7 / 9.0 ,5 / 8.0;

    points1.transposeInPlace();
    points2.transposeInPlace();
    points3.transposeInPlace();

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points(3);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector, control_points);

    Eigen::Matrix<double, 3, Eigen::Dynamic> points11(3, 3);
    points11 <<  2 / 7.0, 1 / 7.0, 2 / 5.0,
               3.0 / 5, 1.0 / 10, 2.0 / 3,
               1 , 0 , 4.0 / 5;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points12(3, 3);
    points12 << 3 / 8.0, 4.0 / 9, 2.0 / 3,
               1.0 / 3, 1.0 / 2, 5.0 / 5,
               5.0 / 7, 3.0 / 8, 2.0 / 7;

    Eigen::Matrix<double, 3, Eigen::Dynamic> points13(3, 3);
    points13 << 1 / 5.0, 7 / 6.0, 3 / 7.0,
               3.0 / 4, 7.0 / 8, 5.0 / 8,
               7 / 8.0, 4 / 7.0, 5.0 / 10;

    points11.transposeInPlace();
    points12.transposeInPlace();
    points13.transposeInPlace();

    Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> control_points2(3);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    nurbs_surface<double, 3, -1, -1, -1, -1, false>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, false>(u_knots_vector, v_knots_vector, control_points2);
    _set_error_mode(_OUT_TO_MSGBOX);
    nurbs_surfaces_intersect<nunbs_surface3d, nunbs_surface3d> ts;
    ts.init(test_surface, test_surface2);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest10.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    read_chat_points(test_data, src_path + "BEZIER_INTtest10.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}

TEST(BEZIER_INT, test11)
{
    Eigen::VectorX<double> u_knots_vector2(8);
    u_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector2(8);
    v_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;

    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector << 0, 0, 0, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 1, 1, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 3);
    points2 << 3.0 / 8, 4.0 / 9, 2.0 / 3,
        2.0 / 3, 3.0 / 4, 1 / 3.0,
        6.0 / 7, 3.0 / 8, 5.0 / 7,
        3.123, 2.12312, 1.23;


    nurbs_curve<double, 3, true, -1, -1>* test_curve = new nurbs_curve<double, 3, true, -1, -1>(u_knots_vector, points2);

    Eigen::Matrix<double, Eigen::Dynamic, 4> points11(3, 4);
    points11 << 2 / 7.0, 1 / 7.0, 2 / 5.0, 2.1,
        3.0 / 5, 1.0 / 10, 2.0 / 3, 2.4,
        1, 0, 4.0 / 5, 1.3;
    Eigen::Matrix<double, 4, Eigen::Dynamic> points21 = points11.transpose();

    Eigen::Matrix<double, Eigen::Dynamic, 4> points12(3, 4);
    points12 << 3 / 8.0, 4.0 / 9, 2.0 / 3,1.3, 
               1.0 / 3, 1.0 / 2, 5.0 / 5,2.21,
               5.0 / 7, 3.0 / 8, 2.0 / 7, 2.43;
               
    Eigen::Matrix<double, 4, Eigen::Dynamic> points22 = points12.transpose();

    Eigen::Matrix<double, Eigen::Dynamic, 4> points13(3, 4);
    points13 << 1 / 5.0, 7 / 6.0, 3 / 7.0,2.13, 
               3.0 / 4, 7.0 / 8, 5.0 / 8,3.24,
               7 / 8.0, 4 / 7.0, 5.0 / 10, 3.24;
               
    Eigen::Matrix<double, 4, Eigen::Dynamic> points23 = points13.transpose();
    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points2(3);
    control_points2(0) = points21;
    control_points2(1) = points22;
    control_points2(2) = points23;
    nurbs_surface<double, 3, -1, -1, -1, -1, true>* test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, true>(u_knots_vector, v_knots_vector, control_points2);
    // save_obj2(*test_curve, "test_curve.obj");
    // save_obj2(*test_surface2, "test_surface2.obj");
    _set_error_mode(_OUT_TO_MSGBOX);
    std::vector<Eigen::Vector<double, 3>> int_point = bezier_curve_int_bezier_surface<double, 3, true>(*test_curve, *test_surface2);
}

TEST(BEZIER_INT, test12)
{
    Eigen::VectorX<double> u_knots_vector2(8);
    u_knots_vector2<< 0, 0, 0, 0, 1, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector2(8);
    v_knots_vector2 << 0, 0, 0, 0, 1, 1, 1, 1;

    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector<< 0, 0, 0, 1, 1, 1;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 1, 1, 1;


    Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 3);
    points1 << -10, 0, 10,
               -30, -30, -30,
               3, -3, 3,
               1.32, 1.231, 1.321;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 3);
    points2 << -10, 0, 10,
               -10, -10, -10,
               3, -3, 3,
               1.32, 1.312, 1.3213;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 3);
    points3 << -10, 0, 10,
               10, 10, 10,
               3, -3, 3,
               0.9423, 1.3213, 0.9423;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points4(4, 3);
    points4 << -10, 0, 10,
               30,  30, 15,
               3,  -3,  3,
               0.9423, 1.32134, 1.4234;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(4);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    control_points(3) = points4;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> *test_surface = new nurbs_surface<double, 3, -1, -1, -1, -1, true>(u_knots_vector, v_knots_vector2, control_points);
    
    Eigen::Matrix<double, 4, Eigen::Dynamic> points11(4, 4);
    points11 << 20,  10, -10, -20,
                -10, -10, -10, -10,
                3,  3, 3, 3,
               0.7423, 1.3213, 0.9423, 1.32;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points12(4, 4);
    points12 << 20,  10, -10, -20,
                -10, -10, -10, -10,
                -1,  -1, -1, -1,
               0.8423, 1.32134, 2.4234, 1.32;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points13(4, 4);
    points13 << 20,  10, -10, -20,
                10, 10, 10, 10,
                -1,  -1, -1, -1,
               2.32, 1.231, 3.321, 1.32;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points14(4, 4);
    points14 << 20,  10, -10, -20,
                10, 10,  10,  10,
                3,  3, 3, 3,
               1.32, 2.312, 1.3213, 1.32;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points2(4);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    control_points2(3) = points14;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> *test_surface2 = new nurbs_surface<double, 3, -1, -1, -1, -1, true>(u_knots_vector2, v_knots_vector2, control_points2);
    _set_error_mode(_OUT_TO_MSGBOX);
    nurbs_surfaces_intersect<nurbs_surface3d, nurbs_surface3d> ts;
    ts.init(test_surface, test_surface2);
    ts.surafces_intersection2();
    surf_surf_int<double, 3> intersection = ts.m_result;
    // save_chat_points(intersection, "BEZIER_INTtest12.json");
    size_t curve_count = intersection.m_int_chats.size();
    size_t isolate_point_count = intersection.m_isolate_points.size();
    surf_surf_int<double, 3> test_data;
    read_chat_points(test_data, src_path + "BEZIER_INTtest12.json");
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
            flag = ts.is_point_in_intcurve(int_chat.m_inter_points[index].m_uv, int_chat.m_is_transversal, current_index);
            EXPECT_TRUE(flag);
            EXPECT_EQ(current_index, curve_index);
        }
    }
}  

TEST(BEZIER_INT, test13)
{
    Eigen::VectorX<double> u_knots_vector2(10);
    u_knots_vector2<< 0, 0, 0, 0, 1, 2, 3.123, 3.123, 3.123, 3.123;

    Eigen::VectorX<double> v_knots_vector2(10);
    v_knots_vector2 << 0, 0, 0, 0, 1.2, 2.0, 2.3, 2.3, 2.3, 2.3;

    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 1, 1, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 6);
    points2 << -10, 0, 10, 0, -10, -12,
               -10, -10, 0, 0, 0, 0,
               3, -3, 3, 4, 5, 3,
               1.32, 1.312, 1.3213, 1.0, 1.5, 1.0;

    nurbs_curve<double, 3, true, -1, -1> test_curve(u_knots_vector2, points2);
    
    Eigen::Matrix<double, 4, Eigen::Dynamic> points11(4, 6);
    points11 << 20,  10, -10, -20, -25, -30,
                -10, -10, -10, -10, -10, -10,
                3,  3, 3, 3, 3, 3,
               0.7423, 1.3213, 0.9423, 1.32, 1.0, 1.0;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points12(4, 6);
    points12 << 20,  10, -10, -20, -25, -30,
                -10, -10, -10, -10, -10, -10,
                -1,  -1, -1, -1, -1, -1,
               0.8423, 1.32134, 2.4234, 1.32, 1.0, 1.0;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points13(4, 6);
    points13 << 20,  10, -10, -20, -22, -23,
                10, 10, 10, 10, 10, 10,
                -1,  -1, -1, -1, -1, -1,
               2.32, 1.231, 3.321, 1.32, 1.0, 1.0;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points2(3);
    control_points2(0) = points11;
    control_points2(1) = points12;
    control_points2(2) = points13;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface2(v_knots_vector2, v_knots_vector, control_points2);
    _set_error_mode(_OUT_TO_MSGBOX);

    curve_surface_int<nurbs_curve<double, 3, true, -1, -1>, nurbs_surface<double, 3, -1, -1, -1, -1, true>> curve_surface_intersection;
    curve_surface_intersection.init(&test_curve, &test_surface2);
    curve_surface_intersection.run_intersect();
    std::vector<Eigen::Vector3d> int_points;
    EXPECT_EQ(2, curve_surface_intersection.m_int_points.size());
    EXPECT_NEAR((curve_surface_intersection.m_int_points[0].m_int_point - Eigen::Vector3d(3.1889989613324179, -3.8871574189335987, 0.12885501236201899)).norm(), 0.0, 1e-8);
    EXPECT_NEAR((curve_surface_intersection.m_int_points[1].m_int_point - Eigen::Vector3d(-6.5295318280773662, -7.5569111797466855, 1.6730177355028792)).norm(), 0.0, 1e-8);
}



int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    
    // ::testing::FLAGS_gtest_filter = "BEZIER_INT.test2";
    return RUN_ALL_TESTS();
}
