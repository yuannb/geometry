#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "convert_nubrs_with_polynomial.h" 
#include "create_nurbs_arc.h"
#include "contruct_primitive_nurbs_surface.h"
#include "fit_nurbs.h"
#include "gtest/gtest.h"
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
    }

    // nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    Eigen::Matrix<double, 3, Eigen::Dynamic> pointss;
    Eigen::Matrix<double, 3, Eigen::Dynamic> ders;
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


TEST_F(CreateNurbsCurve, InterpolateWithEndsTangent)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    Eigen::Vector3d D0{-1, 9, 0};
    Eigen::Vector3d D1{-5, -5, 0};
    global_curve_interpolate_with_ends_tangent<double, 3, ENPARAMETERIEDTYPE::CHORD>(pointss, D0, D1, 3, new_nurbs);

    ASSERT_TRUE(true == true);
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

TEST_F(CreateNurbsCurve, local_2degree_interpolate)
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
    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pointss1)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();

    std::string dir("view2.obj");
    std::ofstream outfile(dir);
    for (int index = 0; index < 5; ++index)
    {
        Eigen::Vector3d point = pointss.col(index);
        outfile << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;    
    }
    outfile.close();

    // Eigen::Vector3d point;
    // for (int index = 0; index < 5; ++index)
    // {
    //     Eigen::Vector<Eigen::Vector3d, 2> point;
    //     double u;
    //     Eigen::Vector3d p;
    //     ENUM_NURBS flag = new_nurbs.find_nearst_point_on_curve(pointss.col(index), u, p);
    //     ASSERT_EQ(ENUM_NURBS::NURBS_SUCCESS, flag);
    //     new_nurbs.derivative_on_curve<1>(u, point);

    //     double distance = (point[0] - pointss.col(index)).norm();
    //     EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);

    //     distance = (point[1] - ders.col(index)).norm();
    //     EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);
    // }

}

TEST_F(CreateNurbsCurve, local_arc_interpolate)
{
    nurbs_curve<double, 3, true, -1, -1> new_nurbs;
    local_2degree_arc_interpolate<double, 3>(pointss, ders, new_nurbs);

    std::vector<Eigen::Vector3d> pointss1;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;   
        new_nurbs.point_on_curve(0.01 * i, point);
        pointss1.push_back(point);
    }
    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pointss1)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();

    std::string dir("view2.obj");
    std::ofstream outfile(dir);
    for (int index = 0; index < 5; ++index)
    {
        Eigen::Vector3d point = pointss.col(index);
        outfile << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;    
    }
    outfile.close();

    // Eigen::Vector3d point;
    // for (int index = 0; index < 5; ++index)
    // {
    //     Eigen::Vector<Eigen::Vector3d, 2> point;
    //     double u;
    //     Eigen::Vector3d p;
    //     ENUM_NURBS flag = new_nurbs.find_nearst_point_on_curve(pointss.col(index), u, p);
    //     ASSERT_EQ(ENUM_NURBS::NURBS_SUCCESS, flag);
    //     new_nurbs.derivative_on_curve<1>(u, point);

    //     double distance = (point[0] - pointss.col(index)).norm();
    //     EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);

    //     distance = (point[1] - ders.col(index)).norm();
    //     EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);
    // }

}

TEST_F(CreateNurbsCurve, local_arc_interpolate2)
{
    nurbs_curve<double, 3, true, -1, -1> new_nurbs;
    local_2degree_arc_interpolate<double, 3, true>(pointss, new_nurbs);

    std::vector<Eigen::Vector3d> pointss1;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;   
        new_nurbs.point_on_curve(0.01 * i, point);
        pointss1.push_back(point);
    }
    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pointss1)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();

    std::string dir("view2.obj");
    std::ofstream outfile(dir);
    for (int index = 0; index < 5; ++index)
    {
        Eigen::Vector3d point = pointss.col(index);
        outfile << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;    
    }
    outfile.close();

    // Eigen::Vector3d point;
    // for (int index = 0; index < 5; ++index)
    // {
    //     Eigen::Vector<Eigen::Vector3d, 2> point;
    //     double u;
    //     Eigen::Vector3d p;
    //     ENUM_NURBS flag = new_nurbs.find_nearst_point_on_curve(pointss.col(index), u, p);
    //     ASSERT_EQ(ENUM_NURBS::NURBS_SUCCESS, flag);
    //     new_nurbs.derivative_on_curve<1>(u, point);

    //     double distance = (point[0] - pointss.col(index)).norm();
    //     EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);

    //     distance = (point[1] - ders.col(index)).norm();
    //     EXPECT_NEAR(distance, 0.0, DEFAULT_ERROR);
    // }

}


TEST_F(CreateNurbsSurface, InterpolateSurface)
{
    nurbs_surface<double, 3, -1, -1, -1, -1, false> new_nurbs;

    global_surface_interpolate<double, 3, ENPARAMETERIEDTYPE::CHORD>(points, 3, 4, new_nurbs);

    std::vector<Eigen::Vector3d> pss;

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

