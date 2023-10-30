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


int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    // ::testing::FLAGS_gtest_filter = "CreateNurbsCurve2.InterpolateWithEndsTangent4";
    return RUN_ALL_TESTS();
}
