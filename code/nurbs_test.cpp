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
// #include <memory>
using namespace tnurbs;

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

TEST_F(CreateNurbsCurve2, InterpolateWithEndsTangent2)
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
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d p;
        double u;
        test_nurbs.find_nearst_point_on_curve(points.col(i), u, p);
        double d = (p - points.col(i)).norm();
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
    for (int i = 0; i < 10; ++i)
    {
        Eigen::Vector3d p;
        double u;
        test_nurbs.find_nearst_point_on_curve(points.col(i), u, p);
        double d = (p - points.col(i)).norm();
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
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            Eigen::Vector3d point;
            double u, v;
            Eigen::Vector3d p;
            ENUM_NURBS flag = new_nurbs.find_nearst_point_on_surface(points[i].col(j), u, v, p);
            ASSERT_EQ(ENUM_NURBS::NURBS_SUCCESS, flag);
            new_nurbs.point_on_surface(u, v, point);

            double distance = (point - points[i].col(j)).norm();
            EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);
        }
    }

}



TEST_F(CreateNurbsSurface, InterpolateSurface)
{
    nurbs_surface<double, 3, -1, -1, -1, -1, false> new_nurbs;

    global_surface_interpolate<double, 3, ENPARAMETERIEDTYPE::CHORD>(points, 3, 4, new_nurbs);

    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            Eigen::Vector3d point;
            double u, v;
            Eigen::Vector3d p;
            ENUM_NURBS flag = new_nurbs.find_nearst_point_on_surface(points[i].col(j), u, v, p);
            ASSERT_EQ(ENUM_NURBS::NURBS_SUCCESS, flag);
            new_nurbs.point_on_surface(u, v, point);

            double distance = (point - points[i].col(j)).norm();
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


    std::vector<Eigen::Vector3d> pointss;
    std::vector<Eigen::Vector3d> points22;
    Eigen::Vector<Eigen::Vector3d, 2> normals;
    normals[0] = Eigen::Vector3d(0, 0, 1);
    normals[1] = Eigen::Vector3d(0, 0, -1);
    for (int i = 30; i < 35; ++i)
    {
        for (int j = 50; j < 55; ++j)
        {
            std::cout << i << " " << j << std::endl;
            Eigen::Vector3d point;
            test_surface.point_on_surface(0.05 * i, 0.03 * j, point);
            Eigen::Vector3d p = point + 1 * normals[i % 2], nearst_point;
            pointss.push_back(p);
        }

    }
    std::vector<Eigen::Vector3d> nearst_points;
    std::vector<double> us, vs;
    clock_t start_time = clock();
    find_nearst_point_on_surface(test_surface, pointss, us, vs, nearst_points);
    
    clock_t start_time2 = clock();
    for (int i = 30; i < 35; ++i)
    {
        for (int j = 50; j < 55; ++j)
        {
            Eigen::Vector3d point;
            test_surface.point_on_surface(0.05 * i, 0.03 * j, point);
            Eigen::Vector3d p = point + 1 * normals[i % 2], nearst_point;
            double u_p, v_p;
            test_surface.find_nearst_point_on_surface(p, u_p, v_p, nearst_point);
            double d1 = (nearst_point - p).norm();
            double d2 = (nearst_points[(i - 30) * 5 + (j - 50)] - p).norm();
            double dis = d1 - d2;
            std::cout << dis << std::endl;
            ASSERT_TRUE(std::abs(dis) < 1e-4);
            
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
    std::vector<Eigen::Vector3d> pointss;
    std::vector<Eigen::Vector3d> points2;
    Eigen::Vector<Eigen::Vector3d, 2> normals;
    normals[0] = Eigen::Vector3d(1, -1, 0);
    normals[1] = Eigen::Vector3d(-1, 1, 0);
    std::vector<double> us;
    std::vector<Eigen::Vector3d> nearst_points;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector3d p = point + 1 * normals[i % 2], nearst_point;
        pointss.push_back(p);
    }
    find_nearst_point_on_curve(curve1, pointss, us, nearst_points);
    
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector3d p = point + 1 * normals[i % 2], nearst_point;
        double u_p;
        curve1.find_nearst_point_on_curve(p, u_p, nearst_point);
        double d1 = (nearst_point - p).norm();
        double d2 = (nearst_points[i] - p).norm();
        std::cout << (d1 - d2) << std::endl;
        std::cout << i << std::endl;
        ASSERT_TRUE(std::abs(d1 - d2) < 1e-4);
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

  


int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    
    ::testing::FLAGS_gtest_filter = "CreateNurbsCurve2.FindNearstPoint*";
    return RUN_ALL_TESTS();
}
