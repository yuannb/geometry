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
        create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, nurbs);
        for (int i = 0; i < 5; ++i)
        {
            Eigen::Vector3d point;   
            nurbs.point_on_curve(0.2 * i, point);
            pointss.col(i) = point;
        }
    }

    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    Eigen::Matrix<double, 3, Eigen::Dynamic> pointss;
};


TEST_F(CreateNurbs, InterpolateWithEndsTangent)
{
    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    Eigen::Vector3d D0{-1, 9, 0};
    Eigen::Vector3d D1{-5, -5, 0};
    global_curve_interpolate_with_ends_tangent<double, 3, ENPARAMETERIEDTYPE::CHORD>(pointss, D0, D1, 3, new_nurbs);

    std::vector<Eigen::Vector3d> new_points;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        new_nurbs.point_on_curve(0.01 * i + 0.0099999, point);
        new_points.push_back(point);
    }

    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : new_points)
    {
        outfile << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    std::cout << "google_test_1" << std::endl;
    ASSERT_TRUE(true == false);
}

