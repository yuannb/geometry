#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "convert_nubrs_with_polynomial.h" 
#include "create_nurbs_arc.h"
#include "contruct_primitive_nurbs_surface.h"
#include "fit_nurbs.h"
#include "gtest/gtest.h"
using namespace tnurbs;

class CreateNurbs : public testing::Test
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


TEST_F(CreateNurbs, InterpolateWithEndsTangent)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    Eigen::Vector3d D0{-1, 9, 0};
    Eigen::Vector3d D1{-5, -5, 0};
    global_curve_interpolate_with_ends_tangent<double, 3, ENPARAMETERIEDTYPE::CHORD>(pointss, D0, D1, 3, new_nurbs);

    ASSERT_TRUE(true == true);
}


TEST_F(CreateNurbs, InterpolateWithEndsTangent2)
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

TEST_F(CreateNurbs, InterpolateHermite)
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

