#include "build.h"
#include "query.h"
#include <iostream>
#include <fstream>
#include "planeequ.h"
// #include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "split.h"
#include <sstream>
#include "delaunay.h"
#include "nurbs_tool.h"
#include "nurbs_curve.h"
#include "bezier_curve.h"
#include "bezier_surface.h"
#include "nurbs_surface.h"
#include "time.h"
#include "polynomial_curve.h"
#include "convert_nubrs_with_polynomial.h" 
#include "create_nurbs_arc.h"
#include "contruct_primitive_nurbs_surface.h"
#include "fit_nurbs.h"
using namespace tnurbs;
// void test_DeCasteljaul_t()
// {
//     Eigen::Vector<double, 3> v1{0, 0, 0};
//     Eigen::Vector<double, 3> v2{1, 0, 0};
//     Eigen::Vector<double, 3> v3{1, 1, 0};
//     std::vector<Eigen::Vector<double, 3>> point {v1, v2, v3};
//     bezier_curve<double, 3, false, -1> curve1(point);

//     std::vector<Eigen::Vector3d> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         Eigen::Vector3d point;
//         curve1.point_on_curve(0.01 * i, point);
//         pointss.push_back(point);
//     }
//     //     // // write doc
//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (auto point : pointss)
//     {
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }

// void test_DeCasteljaul_t_t()
// {
//     Eigen::Vector<double, 3> v1{0, 0, 0};
//     Eigen::Vector<double, 3> v2{1, 0, 0};
//     Eigen::Vector<double, 3> v3{1, 1, 0};
//     Eigen::Vector<Eigen::Vector<double, 3>, 3> point {v1, v2, v3};
//     bezier_curve<double, 3, false, 3> curve1(point);

//     std::vector<Eigen::Vector3d> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         Eigen::Vector3d point;
//         curve1.point_on_curve(0.01 * i, point);
//         pointss.push_back(point);
//     }
//     //     // // write doc
//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (auto point : pointss)
//     {
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }

// void test_AllBernstein_t()
// {
//     Eigen::Vector<double, 3> v1{0, 0, 0};
//     Eigen::Vector<double, 3> v2{1, 0, 0};
//     Eigen::Vector<double, 3> v3{1, 1, 0};
//     std::vector<Eigen::Vector<double, 3>> point {v1, v2, v3};
//     bezier_curve<double, 3, false, -1> curve1(point);

//     std::vector<Eigen::Vector3d> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         Eigen::Vector3d point;
//         curve1.point_on_curve_b(0.01 * i, point);
//         pointss.push_back(point);
//     }
//     //     // // write doc
//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (auto point : pointss)
//     {
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }

// void test_AllBernstein_t_t()
// {
//     Eigen::Vector<double, 3> v1{0, 0, 0};
//     Eigen::Vector<double, 3> v2{1, 0, 0};
//     Eigen::Vector<double, 3> v3{1, 1, 0};
//     Eigen::Vector<Eigen::Vector<double, 3>, 3> point {v1, v2, v3};
//     bezier_curve<double, 3, false, 3> curve1(point);

//     std::vector<Eigen::Vector3d> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         Eigen::Vector3d point;
//         curve1.point_on_curve_b(0.01 * i, point);
//         pointss.push_back(point);
//     }
//     //     // // write doc
//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (auto point : pointss)
//     {
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }


// void test_DeCasteljaul2_t()
// {
//     Eigen::Vector<double, 3> v1{0, 0, 0};
//     Eigen::Vector<double, 3> v2{1, 0, 0};
//     Eigen::Vector<double, 3> v3{1, 1, 0};

//     Eigen::Vector<double, 3> v4{0, 0, 10};
//     Eigen::Vector<double, 3> v5{1, 0, 10};
//     Eigen::Vector<double, 3> v6{1, 1, 10};
//     std::vector<Eigen::Vector<double, 3>> point {v1, v2, v3, v4, v5, v6};
//     bezier_surface<double, 3, -1, -1, false> curve1(2, 3, point);

//     std::vector<Eigen::Vector3d> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         for (int j = 0; j < 100; ++j)
//         {
//             Eigen::Vector3d point;
//             curve1.point_on_surface((double)0.01 * i,(double) 0.01 * j, point);
//             pointss.push_back(point);
//         }

//     }
//     //     // // write doc
//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (auto point : pointss)
//     {
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }

// void test_DeCasteljaul2_t_t()
// {
//     Eigen::Vector<double, 3> v1{0, 0, 0};
//     Eigen::Vector<double, 3> v2{1, 0, 0};
//     Eigen::Vector<double, 3> v3{1, 1, 0};

//     Eigen::Vector<double, 3> v4{0, 0, 10};
//     Eigen::Vector<double, 3> v5{1, 0, 10};
//     Eigen::Vector<double, 3> v6{1, 1, 10};
//     std::vector<Eigen::Vector<double, 3>> point {v1, v2, v3, v4, v5, v6};
//     bezier_surface<double, 3, 2, 3, false> curve1(point);

//     std::vector<Eigen::Vector3d> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         for (int j = 0; j < 100; ++j)
//         {
//             Eigen::Vector3d point;
//             curve1.point_on_surface((double)0.01 * i,(double) 0.01 * j, point);
//             pointss.push_back(point);
//         }

//     }
//     //     // // write doc
//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (auto point : pointss)
//     {
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }

// void test_basis_functions_t_t()
// {
//     Eigen::Vector<double, 3> v1{0, 0, 0};
//     Eigen::Vector<double, 3> v2{1, 0, 0};
//     Eigen::Vector<double, 3> v3{1, 1, 0};
//     Eigen::Matrix<double, 3, 3> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;

//     Eigen::Vector<double, 6> knots_vector{0, 0, 0, 1, 1, 1};
//     nurbs_curve<double, 3, false, 3, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector3d> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         Eigen::Vector3d point;
//         curve1.point_on_curve((double)0.01 * i, point);
//         pointss.push_back(point);
//     }
//     //     // // write doc
//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (auto point : pointss)
//     {
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }


// void test_basis_functions_t()
// {
//     Eigen::Vector<double, 3> v1{0, 0, 0};
//     Eigen::Vector<double, 3> v2{1, 0, 0};
//     Eigen::Vector<double, 3> v3{1, 1, 0};
//     Eigen::Matrix<double, 3, 3> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;

//     Eigen::Vector<double, 6> knots_vector{0, 0, 0, 1, 1, 1};
//     nurbs_curve<double, 3, false, -1, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector3d> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         Eigen::Vector3d point;
//         curve1.point_on_curve((double)0.01 * i, point);
//         pointss.push_back(point);
//     }
//     //     // // write doc
//     std::string dir("view1.obj");
//     std::ofstream outfile(dir);

//     for (auto point : pointss)
//     {
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }


// void test_basis_functions()
// {
//     Eigen::Vector<double, 3> v1{0, 0, 0};
//     Eigen::Vector<double, 3> v2{1, 0, 0};
//     Eigen::Vector<double, 3> v3{1, 1, 0};
//     Eigen::Matrix<double, 3, 3> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;

//     Eigen::Vector<double, 6> knots_vector{0, 0, 0, 1, 1, 1};
//     nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector3d> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         Eigen::Vector3d point;
//         curve1.point_on_curve((double)0.01 * i, point);
//         pointss.push_back(point);
//     }
//     //     // // write doc
//     std::string dir("view2.obj");
//     std::ofstream outfile(dir);

//     for (auto point : pointss)
//     {
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }

// void ders_basis_funs_t()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, 8, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::Vector<Eigen::Vector<double, 3>, Eigen::Dynamic> reuslt(2);
//         curve1.derivative_on_curve(i * 0.005, 1, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void ders_basis_funs_t_t()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, 8, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::Vector<Eigen::Vector<double, 3>, 2> reuslt;
//         curve1.derivative_on_curve<1>(i * 0.005, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void ders_basis_funs_t_1()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, -1, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::Vector<Eigen::Vector<double, 3>, Eigen::Dynamic> reuslt(2);
//         curve1.derivative_on_curve(i * 0.005, 1, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void ders_basis_funs_t_t_1()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, -1, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::Vector<Eigen::Vector<double, 3>, Eigen::Dynamic> reuslt(2);
//         curve1.derivative_on_curve<1>(i * 0.005, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void ders_basis_funs_t_2()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::Vector<Eigen::Vector<double, 3>, Eigen::Dynamic> reuslt(2);
//         curve1.derivative_on_curve(i * 0.005, 1, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void ders_basis_funs_t_t_2()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::Vector<Eigen::Vector<double, 3>, Eigen::Dynamic> reuslt(2);
//         curve1.derivative_on_curve<1>(i * 0.005, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void curve_derivs_alg2_t_t()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, 8, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::Vector<Eigen::Vector<double, 3>, 2> reuslt;
//         curve1.curve_derivs_alg2<1>(i * 0.005, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void curve_derivs_alg2_t()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, 8, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::VectorX<Eigen::Vector<double, 3>> reuslt(2);
//         curve1.curve_derivs_alg2(1, i * 0.005, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void curve_derivs_alg2_t_2()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, -1, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::Vector<Eigen::Vector<double, 3>, 2> reuslt;
//         curve1.curve_derivs_alg2<1>(i * 0.005, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void curve_derivs_alg2_2()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, -1, 2> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::VectorX<Eigen::Vector<double, 3>> reuslt(2);
//         curve1.curve_derivs_alg2(1, i * 0.005, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));
//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void curve_derivs_alg2_t_3()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::Vector<Eigen::Vector<double, 3>, 2> reuslt;
//         curve1.curve_derivs_alg2<1>(i * 0.005, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));

//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

// void curve_derivs_alg2_3()
// {
//     Eigen::Vector<double, 3> v1{0.5, -0.5, 1};
//     Eigen::Vector<double, 3> v2{1, 0, 2};
//     Eigen::Vector<double, 3> v3{1, 1, 3};
//     Eigen::Vector<double, 3> v4{-1, 0, 3};
//     Eigen::Vector<double, 3> v5{0, -1, 5};
//     Eigen::Vector<double, 3> v6{2, 3, 6};
//     Eigen::Vector<double, 3> v7{4, 3, 10};
//     Eigen::Vector<double, 3> v8{5, 6, 13};
//     Eigen::Matrix<double, 3, 8> mat;
//     mat.col(0) = v1;
//     mat.col(1) = v2;
//     mat.col(2) = v3;
//     mat.col(3) = v4;
//     mat.col(4) = v5;
//     mat.col(5) = v6;
//     mat.col(6) = v7;
//     mat.col(7) = v8;

//     Eigen::Vector<double, 11> knots_vector{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
//     nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);

//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 1000; ++i)
//     {
//         Eigen::VectorX<Eigen::Vector<double, 3>> reuslt(2);
//         curve1.curve_derivs_alg2(1, i * 0.005, reuslt);
//         pointss.push_back(reuslt(0));
//         pointss.push_back(reuslt(1) + reuslt(0));
//     }

//     //     // // write doc
//     std::string dir1("view1.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size() / 2;
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index * 2];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     std::string dir("view.obj");
//     std::ofstream outfile(dir);

//     for (int index = 0; index < 100; ++index)
//     {
//         auto point = pointss[index * 20];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;

//         point = pointss[index * 20 + 1];
//         outfile << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }

//     for (int index = 0; index < 99; ++index)
//     {
//         outfile << "l " <<  2 * index + 1 << " " <<
//           2 * index +  2 << std::endl;
//     }
// }

void nurbs_surface_point()
{
    Eigen::VectorX<double> u_knots_vector(10);
    u_knots_vector<< 0, 0, 0, 1, 2, 3, 4, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(8);
    v_knots_vector << 0, 0, 0, 1, 2, 3, 3, 3;
    Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 7);
    points1 <<0, 0, 0, 4, 4, 0, 0,
              0, 0, 2, 6, 2, 0, 0,
              0, 0, 4, 8, 4, 0, 0,
              1, 1, 1, 2, 1, 1, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 7);
    points2 << 0, 0, 0, 12, 8, 0, 0,
               0, 0, 6, 24, 6, 0, 0,
               0, 0, 4, 12, 4, 0, 0,
               1, 1, 2, 6,  2, 1, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 7);
    points3 <<0, 0, 0, 4, 4, 0, 0,
              0, 0, 2, 6, 2, 0, 0,
              0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 2, 1, 1, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points4(4, 7);
    points4 <<0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 1, 1, 1, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points5(4, 7);
    points5 <<0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 1, 1, 1, 1;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(5);
    control_points(0) = points4;
    control_points(4) = points5;
    control_points(1) = points1;
    control_points(2) = points2;
    control_points(3) = points3;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    Eigen::Vector3d vec;
    test_surface.point_on_surface(5.0 / 2.0, 1, vec);
    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss.push_back(reuslt);
        }
    }

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

// void nurbs_surface_point_1()
// {
//     Eigen::VectorX<double> u_knots_vector(6);
//     u_knots_vector<< 0, 0, 0, 5, 5, 5;
//     Eigen::VectorX<double> v_knots_vector(6);
//     v_knots_vector << 0, 0, 0, 3, 3, 3;
//     Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 3);
//     points1 << 0, 0, 0,
//                2, 6, 2,
//                5, 4, 0,
//                1, 2, 1;

//     Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 3);
//     points2 << 4, 12, 4,
//                6, 24, 6,
//                8, 12, 0,
//                2, 6,  2;

//     Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 3);
//     points3 << 4, 8, 4,
//                2, 6, 2,
//                4, 4, 0,
//                1, 2, 1;

//     Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(3);
//     control_points(0) = points1;
//     control_points(1) = points2;
//     control_points(2) = points3;
//     nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface(u_knots_vector, v_knots_vector, control_points);
//     std::vector<Eigen::Vector<double, 3>> pointss;
//     for (int i = 0; i < 100; ++i)
//     {
//         for (int j = 0; j < 100; ++j)
//         {
//             Eigen::Vector<double, 3> reuslt;
//             test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
//             pointss.push_back(reuslt);
//         }
//     }

//     //     // // write doc
//     std::string dir1("view.obj");
//     std::ofstream outfile1(dir1);
//     int pointsCount = pointss.size();
//     for (int index = 0; index < pointsCount; ++index)
//     {
//         auto point = pointss[index];
//         outfile1 << "v " << point[0] << " " <<
//            point[1] << " " << point[2] << std::endl;
//     }
// }


void nurbs_surface_ders_point()
{
    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector<< 0, 0, 0, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 3, 3, 3;
    Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 3);
    points1 << 0, 0, 0,
               2, 6, 2,
               5, 4, 0,
               1, 2, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 3);
    points2 << 4, 12, 4,
               6, 24, 6,
               8, 12, 0,
               2, 6,  2;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 3);
    points3 << 4, 8, 4,
               2, 6, 2,
               4, 4, 0,
               1, 2, 1;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(3);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    std::vector<Eigen::Vector<double, 3>> pointss;
    clock_t start_time = clock();
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Matrix<Eigen::Vector<double, 3>, 2, 2> tanget_vector;
            test_surface.derivative_on_surface<1>(i * 0.05, j * 0.03, tanget_vector);
            pointss.push_back(tanget_vector(0, 0));
            tanget_vector(1, 0) += tanget_vector(0, 0);
            pointss.push_back(tanget_vector(1, 0));
        }
    }
    clock_t end_time=clock();
    std::cout << "The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();

    for (int index = 0; index < pointsCount / 20; ++index)
    {
        auto point = pointss[20 * index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
        point = pointss[20 * index + 1];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir("view1.obj");
     std::ofstream outfile(dir);
    for (int index = 0; index < pointsCount / 2; ++index)
    {
        auto point = pointss[2 * index];
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    for (int index = 0; index < pointsCount / 20 - 1; ++index)
    {
        outfile1 << "l " << 2 * index + 1 << " " <<
           2 * index + 2 << std::endl;
    }
}

void nurbs_surface_ders_point_1()
{
    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector<< 0, 0, 0, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 3, 3, 3;
    Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 3);
    points1 << 0, 0, 0,
               2, 6, 2,
               5, 4, 0,
               1, 2, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 3);
    points2 << 4, 12, 4,
               6, 24, 6,
               8, 12, 0,
               2, 6,  2;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 3);
    points3 << 4, 8, 4,
               2, 6, 2,
               4, 4, 0,
               1, 2, 1;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(3);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    std::vector<Eigen::Vector<double, 3>> pointss;
    clock_t start_time = clock();
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::MatrixX<Eigen::Vector<double, 3>> tanget_vector;
            test_surface.derivative_on_surface(1, i * 0.05, j * 0.03, tanget_vector);
            pointss.push_back(tanget_vector(0, 0));
            tanget_vector(1, 0) += tanget_vector(0, 0);
            pointss.push_back(tanget_vector(1, 0));
        }
    }
    clock_t end_time=clock();
    std::cout << "The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();

    for (int index = 0; index < pointsCount / 20; ++index)
    {
        auto point = pointss[20 * index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
        point = pointss[20 * index + 1];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir("view1.obj");
     std::ofstream outfile(dir);
    for (int index = 0; index < pointsCount / 2; ++index)
    {
        auto point = pointss[2 * index];
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    for (int index = 0; index < pointsCount / 20 - 1; ++index)
    {
        outfile1 << "l " << 2 * index + 1 << " " <<
           2 * index + 2 << std::endl;
    }
}

void nurbs_surface_ders_point_2()
{
    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector<< 0, 0, 0, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 3, 3, 3;
    Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 3);
    points1 << 0, 0, 0,
               2, 6, 2,
               5, 4, 0,
               1, 2, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 3);
    points2 << 4, 12, 4,
               6, 24, 6,
               8, 12, 0,
               2, 6,  2;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 3);
    points3 << 4, 8, 4,
               2, 6, 2,
               4, 4, 0,
               1, 2, 1;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(3);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    nurbs_surface<double, 3, -1, -1, 2, 2, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    std::vector<Eigen::Vector<double, 3>> pointss;
    clock_t start_time = clock();
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Matrix<Eigen::Vector<double, 3>, 2, 2> tanget_vector;
            test_surface.derivative_on_surface<1>(i * 0.05, j * 0.03, tanget_vector);
            pointss.push_back(tanget_vector(0, 0));
            tanget_vector(1, 0) += tanget_vector(0, 0);
            pointss.push_back(tanget_vector(1, 0));
        }
    }
    clock_t end_time=clock();
    std::cout << "The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount / 20; ++index)
    {
        auto point = pointss[20 * index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
        point = pointss[20 * index + 1];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir("view1.obj");
     std::ofstream outfile(dir);
    for (int index = 0; index < pointsCount / 2; ++index)
    {
        auto point = pointss[2 * index];
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    for (int index = 0; index < pointsCount / 20 - 1; ++index)
    {
        outfile1 << "l " << 2 * index + 1 << " " <<
           2 * index + 2 << std::endl;
    }
}

void nurbs_surface_ders_point_3()
{
    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector<< 0, 0, 0, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 3, 3, 3;
    Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 3);
    points1 << 0, 0, 0,
               2, 6, 2,
               5, 4, 0,
               1, 2, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 3);
    points2 << 4, 12, 4,
               6, 24, 6,
               8, 12, 0,
               2, 6,  2;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 3);
    points3 << 4, 8, 4,
               2, 6, 2,
               4, 4, 0,
               1, 2, 1;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(3);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    nurbs_surface<double, 3, -1, -1, 2, 2, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    std::vector<Eigen::Vector<double, 3>> pointss;
    clock_t start_time = clock();
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::MatrixX<Eigen::Vector<double, 3>> tanget_vector;
            test_surface.derivative_on_surface(1, i * 0.05, j * 0.03, tanget_vector);
            pointss.push_back(tanget_vector(0, 0));
            tanget_vector(1, 0) += tanget_vector(0, 0);
            pointss.push_back(tanget_vector(1, 0));
        }
    }
    clock_t end_time=clock();
    std::cout << "The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();

    for (int index = 0; index < pointsCount / 20; ++index)
    {
        auto point = pointss[20 * index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
        point = pointss[20 * index + 1];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir("view1.obj");
     std::ofstream outfile(dir);
    for (int index = 0; index < pointsCount / 2; ++index)
    {
        auto point = pointss[2 * index];
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    for (int index = 0; index < pointsCount / 20 - 1; ++index)
    {
        outfile1 << "l " << 2 * index + 1 << " " <<
           2 * index + 2 << std::endl;
    }
}

void test_insert_knots()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Matrix<double, 3, 3> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;

    Eigen::Vector<double, 6> knots_vector{0, 0, 0, 1, 1, 1};
    nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);

    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }

    curve1.insert_knots(0.5, 2);
    std::vector<Eigen::Vector3d> pointss2;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        pointss2.push_back(point);
    }
    //     // // write doc
    std::string dir("view.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir1("view1.obj");
    std::ofstream outfile1(dir1);

    for (auto point : pointss2)
    {
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_insert_knots_1()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Matrix<double, 3, 3> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;

    Eigen::Vector<double, 6> knots_vector{0, 0, 0, 1, 1, 1};
    nurbs_curve<double, 3, false, 3, 2> curve1(knots_vector, mat);

    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }
    nurbs_curve<double, 3, false, 5, 2> new_nurbs;
    curve1.insert_knots<2>(0.5, new_nurbs);
    std::vector<Eigen::Vector3d> pointss2;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        new_nurbs.point_on_curve((double)0.01 * i, point);
        pointss2.push_back(point);
    }
    //     // // write doc
    std::string dir("view.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir1("view1.obj");
    std::ofstream outfile1(dir1);

    for (auto point : pointss2)
    {
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_pnt_by_corner_cut()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Matrix<double, 3, 3> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;

    Eigen::Vector<double, 6> knots_vector{0, 0, 0, 1, 1, 1};
    nurbs_curve<double, 3, false, 3, 2> curve1(knots_vector, mat);

    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.pnt_by_corner_cut((double)0.01 * i, point);
        pointss.push_back(point);
    }

    //     // // write doc
    std::string dir("view.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_pnt_by_corner_cut_1()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Matrix<double, 3, 3> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;

    Eigen::Vector<double, 6> knots_vector{0, 0, 0, 1, 1, 1};
    nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);

    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.pnt_by_corner_cut((double)0.01 * i, point);
        pointss.push_back(point);
    }

    //     // // write doc
    std::string dir("view.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_nurbs_surface_insert_knots_1()
{
    Eigen::VectorX<double> u_knots_vector(6);
    u_knots_vector<< 0, 0, 0, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(6);
    v_knots_vector << 0, 0, 0, 3, 3, 3;
    Eigen::Matrix<double, 4, Eigen::Dynamic> points1(4, 3);
    points1 << 0, 0, 0,
               2, 6, 2,
               5, 4, 0,
               1, 2, 1;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points2(4, 3);
    points2 << 4, 12, 4,
               6, 24, 6,
               8, 12, 0,
               2, 6,  2;

    Eigen::Matrix<double, 4, Eigen::Dynamic> points3(4, 3);
    points3 << 4, 8, 4,
               2, 6, 2,
               4, 4, 0,
               1, 2, 1;

    Eigen::VectorX<Eigen::Matrix<double, 4, Eigen::Dynamic>> control_points(3);
    control_points(0) = points1;
    control_points(1) = points2;
    control_points(2) = points3;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    std::vector<Eigen::Vector<double, 3>> pointss;
    test_surface.surface_knots_insert(2.5, 2, ENUM_DIRECTION::V_DIRECTION);
    clock_t start_time = clock();
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            std::cout << i << " " << j << std::endl;
            Eigen::Vector<double, 3> tanget_vector;
            test_surface.point_on_surface(i * 0.05, j * 0.03, tanget_vector);
            pointss.push_back(tanget_vector);
        }
    }
    clock_t end_time=clock();
    std::cout << "The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();

    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_refine_knots_vector_1()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, 4> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::Vector<double, 7> knots_vector{0, 0, 0, 0.4, 1, 1, 1};
    nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    curve1.refine_knots_vector(insert_knots);
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        std::cout << i << std::endl;
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}


void test_refine_knots_vector_2()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, 4> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::Vector<double, 7> knots_vector{0, 0, 0, 0.4, 1, 1, 1};
    nurbs_curve<double, 3, false, 4, 2> curve1(knots_vector, mat);
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    nurbs_curve<double, 3, false, -1, 2> curve2;
    curve1.refine_knots_vector(insert_knots, curve2);
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve2.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_refine_knots_vector_3()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, 4> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::Vector<double, 7> knots_vector{0, 0, 0, 0.4, 1, 1, 1};
    nurbs_curve<double, 3, false, -1, 2> curve1(knots_vector, mat);
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    curve1.refine_knots_vector(insert_knots);
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void nurbs_surface_refine_point()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    nurbs_surface<double, 3, -1, -1, 2, 2, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    Eigen::Vector<double, 3> insert_knots{1.3, 2.3, 3.7};
    Eigen::Vector<double, 3> insert_knots_v{0.8, 1.333, 2};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss.push_back(reuslt);
        }
    }

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void nurbs_surface_refine_point_1()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    Eigen::Vector<double, 3> insert_knots{1.3, 2.3, 3.7};
    Eigen::Vector<double, 3> insert_knots_v{0.8, 1.333, 2};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss.push_back(reuslt);
        }
    }

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_decompose_vector_2()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, 4> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::Vector<double, 7> knots_vector{0, 0, 0, 0.4, 1, 1, 1};
    nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    curve1.refine_knots_vector(insert_knots);
    std::vector<bezier_curve<double, 3, false, -1> *> beziers_curves;
    curve1.decompose_to_bezier(beziers_curves);
    std::vector<std::vector<Eigen::Vector3d>>  pointss(4);
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 100; ++i)
        {
            Eigen::Vector3d point;
            beziers_curves[j]->point_on_curve((double)0.01 * i, point);
            pointss[j].push_back(point);
        }
    }

    for (int i = 0; i < 4; ++i)
    {
            //     // // write doc
        std::string dir("view" + std::to_string(i) + ".obj");
        std::ofstream outfile(dir);
        for (auto point : pointss[i])
        {
            outfile << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile.close();
    }

}

void test_decompose_vector_3()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, 4> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::Vector<double, 7> knots_vector{0, 0, 0, 0.4, 1, 1, 1};
    nurbs_curve<double, 3, false, -1, 2> curve1(knots_vector, mat);
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    curve1.refine_knots_vector(insert_knots);
    std::vector<bezier_curve<double, 3, false, 3> *> beziers_curves;
    curve1.decompose_to_bezier(beziers_curves);
    std::vector<std::vector<Eigen::Vector3d>>  pointss(4);
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 100; ++i)
        {
            Eigen::Vector3d point;
            beziers_curves[j]->point_on_curve((double)0.01 * i, point);
            pointss[j].push_back(point);
        }
    }

    for (int i = 0; i < 4; ++i)
    {
            //     // // write doc
        std::string dir("view" + std::to_string(i) + ".obj");
        std::ofstream outfile(dir);
        for (auto point : pointss[i])
        {
            outfile << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile.close();
    }

}

void test_decompose_vector_4()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, 4> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::Vector<double, 7> knots_vector{0, 0, 0, 0.4, 1, 1, 1};
    nurbs_curve<double, 3, false, 4, 2> curve1(knots_vector, mat);
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    // curve1.refine_knots_vector(insert_knots);
    std::vector<bezier_curve<double, 3, false, 3> *> beziers_curves;
    curve1.decompose_to_bezier(beziers_curves);
    std::vector<std::vector<Eigen::Vector3d>>  pointss(4);
    for (int j = 0; j < 2; ++j)
    {
        for (int i = 0; i < 100; ++i)
        {
            Eigen::Vector3d point;
            beziers_curves[j]->point_on_curve((double)0.01 * i, point);
            pointss[j].push_back(point);
        }
    }

    for (int i = 0; i < 4; ++i)
    {
            //     // // write doc
        std::string dir("view" + std::to_string(i) + ".obj");
        std::ofstream outfile(dir);
        for (auto point : pointss[i])
        {
            outfile << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile.close();
    }

}

void nurbs_surface_decompose_1()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    Eigen::Vector<double, 3> insert_knots{1.3, 2.3, 3.7};
    Eigen::Vector<double, 3> insert_knots_v{0.8, 1.333, 2};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    Eigen::MatrixX<bezier_surface<double, 3 ,-1, -1, true> *> bezier_surfaces;
    test_surface.decompose_to_bezier(bezier_surfaces);
    Eigen::Matrix<std::vector<Eigen::Vector3d>, 4, 4>  pointss;
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 4; ++i)
        {
            for (int k = 0; k <= 100; ++k)
            {
                for (int l = 0; l <= 100; ++l)
                {
                    Eigen::Vector3d point;
                    bezier_surfaces(i, j)->point_on_surface(0.01 * k, 0.01 * l, point);
                    pointss(i, j).push_back(point);
                }

            }
        }
    }

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            //     // // write doc
            std::string dir("view" + std::to_string(i * 4 + j) + ".obj");
            std::ofstream outfile(dir);
            for (auto point : pointss(i, j))
            {
                outfile << "v " << point[0] << " " <<
                point[1] << " " << point[2] << std::endl;
            }
            outfile.close();
        }

    }
}

void nurbs_surface_decompose_2()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    nurbs_surface<double, 3, -1, -1, 2, 2, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    Eigen::Vector<double, 3> insert_knots{1.3, 2.3, 3.7};
    Eigen::Vector<double, 3> insert_knots_v{0.8, 1.333, 2};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    Eigen::MatrixX<bezier_surface<double, 3, 3, 3, true> *> bezier_surfaces;
    test_surface.decompose_to_bezier(bezier_surfaces);
    Eigen::Matrix<std::vector<Eigen::Vector3d>, 4, 4>  pointss;
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 4; ++i)
        {
            for (int k = 0; k <= 100; ++k)
            {
                for (int l = 0; l <= 100; ++l)
                {
                    Eigen::Vector3d point;
                    bezier_surfaces(i, j)->point_on_surface(0.01 * k, 0.01 * l, point);
                    pointss(i, j).push_back(point);
                }
            }
        }
    }

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            //     // // write doc
            std::string dir("view" + std::to_string(i * 4 + j) + ".obj");
            std::ofstream outfile(dir);
            for (auto point : pointss(i, j))
            {
                outfile << "v " << point[0] << " " <<
                point[1] << " " << point[2] << std::endl;
            }
            outfile.close();
        }

    }
}

void test_remove_knots_vector_1()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, 4> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::Vector<double, 7> knots_vector{0, 0, 0, 0.4, 1, 1, 1};
    nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    curve1.refine_knots_vector(insert_knots);
    int time1, time2;
    curve1.remove_knots(0.22222, 3, time1);
    curve1.remove_knots(0.4, 3, time2);
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_remove_knots_vector_2()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, 4> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::Vector<double, 7> knots_vector{0, 0, 0, 0.4, 1, 1, 1};
    nurbs_curve<double, 3, false, 4, 2> curve1(knots_vector, mat);
    int time2;
    nurbs_curve<double, 3, false, -1, 2> curve2;
    curve1.remove_knots(0.4, 3, time2, curve2);
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve2.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void nurbs_surface_reomve_point_1()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    Eigen::Vector<double, 3> insert_knots{1.3, 2.3, 3.7};
    Eigen::Vector<double, 3> insert_knots_v{0.8, 1.333, 2};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    int num1, num2;
    test_surface.remove_knots(1.3, 2, ENUM_DIRECTION::U_DIRECTION, num1);
    test_surface.remove_knots(2.3, 2, ENUM_DIRECTION::U_DIRECTION, num2);
    test_surface.remove_knots(1.333, 2, ENUM_DIRECTION::V_DIRECTION, num1);
    test_surface.remove_knots(0.8, 2, ENUM_DIRECTION::V_DIRECTION, num2);
    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss.push_back(reuslt);
        }
    }

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void nurbs_surface_reomve_point_2()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    nurbs_surface<double, 3, -1, -1, 2, 2, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    Eigen::Vector<double, 3> insert_knots{1.3, 2.3, 3.7};
    Eigen::Vector<double, 3> insert_knots_v{0.8, 1.333, 2};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    int num1, num2;
    test_surface.remove_knots(1.3, 2, ENUM_DIRECTION::U_DIRECTION, num1);
    test_surface.remove_knots(2.3, 2, ENUM_DIRECTION::U_DIRECTION, num2);
    test_surface.remove_knots(1.333, 2, ENUM_DIRECTION::V_DIRECTION, num1);
    test_surface.remove_knots(0.8, 2, ENUM_DIRECTION::V_DIRECTION, num2);
    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss.push_back(reuslt);
        }
    }

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_degree_elevate_curve_1()
{
    Eigen::Vector<double, 3> v1{0, 0, 0};
    Eigen::Vector<double, 3> v2{1, 0, 0};
    Eigen::Vector<double, 3> v3{1, 1, 0};
    Eigen::Vector<double, 3> v4{2, 2, 0};
    Eigen::Matrix<double, 3, 4> mat;
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;

    Eigen::Vector<double, 7> knots_vector{0, 0, 0, 0.4, 1, 1, 1};
    nurbs_curve<double, 3, false, -1, -1> curve1(knots_vector, mat);
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    curve1.refine_knots_vector(insert_knots);
    curve1.degree_elevate(2);
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_degree_reduce_curve_1()
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
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    curve1.refine_knots_vector(insert_knots);
    curve1.degree_elevate(2);
    curve1.degree_reduce();
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_degree_reduce_curve_2()
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
    nurbs_curve<double, 3, false, -1, 2> curve1(knots_vector, mat);
    Eigen::VectorX<double> insert_knots(3);
    insert_knots << 0.22222, 0.4, 0.77777;
    curve1.refine_knots_vector(insert_knots);
    nurbs_curve<double, 3, false, -1, 4> curve2;
    curve1.degree_elevate<2>(curve2);
    nurbs_curve<double, 3, false, -1, 3> curve3;
    nurbs_curve<double, 3, false, -1, 2> curve4;
    curve2.degree_reduce(curve3);
    curve3.degree_reduce(curve4);
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve4.point_on_curve((double)0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}


void nurbs_surface_elevate_degree_1()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    Eigen::Vector<double, 3> insert_knots{1.3, 2.3, 3.7};
    Eigen::Vector<double, 3> insert_knots_v{0.8, 1.333, 2};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    test_surface.degree_elevate(2, ENUM_DIRECTION::U_DIRECTION);
    test_surface.degree_elevate(3, ENUM_DIRECTION::V_DIRECTION);
    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss.push_back(reuslt);
        }
    }

    //     // // write doc
    std::string dir1("view3.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void nurbs_surface_elevate_degree_2()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    nurbs_surface<double, 3, -1, -1, 2, 2, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    Eigen::Vector<double, 3> insert_knots{1.3, 2.3, 3.7};
    Eigen::Vector<double, 3> insert_knots_v{0.8, 1.333, 2};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    nurbs_surface<double, 3, -1, -1, -1, -1, true> new_nurbs_surface;
    test_surface.degree_elevate(2, ENUM_DIRECTION::U_DIRECTION, new_nurbs_surface);
    new_nurbs_surface.degree_elevate(3, ENUM_DIRECTION::V_DIRECTION);
    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            new_nurbs_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss.push_back(reuslt);
        }
    }

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void nurbs_surface_reduce_degree_1()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    Eigen::Vector<double, 3> insert_knots{1.3, 2.3, 3.7};
    Eigen::Vector<double, 3> insert_knots_v{0.8, 1.333, 2};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    test_surface.degree_elevate(2, ENUM_DIRECTION::U_DIRECTION);
    test_surface.degree_elevate(3, ENUM_DIRECTION::V_DIRECTION);
    test_surface.degree_reduce(ENUM_DIRECTION::U_DIRECTION);
    test_surface.degree_reduce(ENUM_DIRECTION::U_DIRECTION);
    test_surface.degree_reduce(ENUM_DIRECTION::V_DIRECTION);
    test_surface.degree_reduce(ENUM_DIRECTION::V_DIRECTION);
    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss.push_back(reuslt);
        }
    }

    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_find_nearst_point_curve_1()
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
    nurbs_curve<double, 3, false, 4, 2> curve1(knots_vector, mat);
    std::vector<Eigen::Vector3d> pointss;
    std::vector<Eigen::Vector3d> points2;
    Eigen::Vector<Eigen::Vector3d, 2> normals;
    normals[0] = Eigen::Vector3d(1, -1, 0);
    normals[1] = Eigen::Vector3d(-1, 1, 0);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector3d p = point + 1 * normals[i % 2], nearst_point;
        double u_p;
        curve1.find_nearst_point_on_curve(p, u_p, nearst_point);
        pointss.push_back(point);
        points2.push_back(p);
        points2.push_back(nearst_point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
    int size_count = points2.size() / 2;
    for (int index = 0; index < size_count; ++index)
    {
        outfile2 << "l " << index * 2 + 1 << " " << index * 2 + 2 << std::endl;
    }
}


void test_find_nearst_point_surface_1()
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
    for (int i = 30; i < 50; ++i)
    {
        for (int j = 50; j < 60; ++j)
        {
            std::cout << i << " " << j << std::endl;
            Eigen::Vector3d point;
            test_surface.point_on_surface(0.05 * i, 0.03 * j, point);
            Eigen::Vector3d p = point + 1 * normals[i % 2], nearst_point;
            double u_p, v_p;
            test_surface.find_nearst_point_on_surface(p, u_p, v_p, nearst_point);
            pointss.push_back(point);
            points22.push_back(p);
            points22.push_back(nearst_point);
        }

    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points22)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
    int size_count = points22.size() / 2;
    for (int index = 0; index < size_count; ++index)
    {
        outfile2 << "l " << index * 2 + 1 << " " << index * 2 + 2 << std::endl;
    }
}

void test_parallel_projection_1()
{
    Eigen::Vector<double, 4> v1{0, 2, 4, 1};
    Eigen::Vector<double, 4> v2{4, 6, 10, 2};
    Eigen::Vector<double, 4> v3{4, 2, 15, 1};
    // Eigen::Vector<double, 4> v4{2, 2, 0};
    Eigen::Matrix<double, 4, Eigen::Dynamic> mat(4, 3);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    // mat.col(3) = v4;

    Eigen::VectorX<double> knots_vector(6);
    knots_vector << 0, 0, 0, 3, 3, 3;
    nurbs_curve<double, 3, true, -1, -1> curve1(knots_vector, mat);
    std::vector<Eigen::Vector3d> pointss;
    std::vector<Eigen::Vector3d> points2;
    Eigen::Vector<Eigen::Vector3d, 2> normals;
    normals[0] = Eigen::Vector3d(1, -1, 0);
    normals[1] = Eigen::Vector3d(-1, 1, 0);
    nurbs_curve<double, 3, true, -1, -1> new_curve;
    curve1.parallel_projection(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 1, 0), new_curve);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.03 * i, point);
        Eigen::Vector3d project_point;
        new_curve.point_on_curve(0.03 * i, project_point);
        pointss.push_back(point);
        points2.push_back(project_point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}


void test_perspective_projection_1()
{
    Eigen::Vector<double, 4> v1{0, 2, 4, 1};
    Eigen::Vector<double, 4> v2{4, 6, 8, 2};
    Eigen::Vector<double, 4> v3{4, 2, 4, 1};
    Eigen::Matrix<double, 4, Eigen::Dynamic> mat(4, 3);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;

    Eigen::VectorX<double> knots_vector(6);
    knots_vector << 0, 0, 0, 3, 3, 3;
    nurbs_curve<double, 3, true, -1, -1> curve1(knots_vector, mat);
    std::vector<Eigen::Vector3d> pointss;
    std::vector<Eigen::Vector3d> points2;
    nurbs_curve<double, 3, true, -1, -1> new_curve;
    curve1.perspective_projection(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(0, 0, 10), new_curve);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.03 * i, point);
        Eigen::Vector3d project_point;
        new_curve.point_on_curve(0.03 * i, project_point);
        pointss.push_back(point);
        points2.push_back(project_point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_reparameter_1()
{
    Eigen::Vector<double, 2> v1{0, 0};
    Eigen::Vector<double, 2> v2{1.0 / 2.0, 1};
    Eigen::Vector<double, 2> v3{1, 0};
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat(2, 3);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;

    Eigen::Vector<double, 1> v12{0};
    Eigen::Vector<double, 1> v23{7.0 / 24.0};
    Eigen::Vector<double, 1> v33{1};
    Eigen::Matrix<double, 1, Eigen::Dynamic> mat1(1, 3);
    mat1.col(0) = v12;
    mat1.col(1) = v23;
    mat1.col(2) = v33;

    Eigen::VectorX<double> knots_vector(6);
    knots_vector << 0, 0, 0, 1, 1, 1;

    // Eigen::VectorX<double> knots_vector_1(6);
    // knots_vector_1 << 1, 1, 1, 4, 4, 4;
    nurbs_curve<double, 2, false, -1, -1> curve1(knots_vector, mat);
    std::vector<Eigen::Vector2d> pointss;
    std::vector<Eigen::Vector2d> points2;
    nurbs_curve<double, 1, false, -1, -1> parameter_function(knots_vector, mat1);
    nurbs_curve<double, 2, false, -1, -1> new_curve;
    curve1.curve_reparameter(parameter_function, new_curve);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector2d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector2d project_point;
        new_curve.point_on_curve(0.01 * i, project_point);
        pointss.push_back(point);
        points2.push_back(project_point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }
}

void test_reparameter_2()
{
    Eigen::Vector<double, 2> v1{0, 0};
    Eigen::Vector<double, 2> v2{0.3, 0.2};
    Eigen::Vector<double, 2> v3{1.0 / 2.0, 1};
    Eigen::Vector<double, 2> v4{0.7, 0.5};
    Eigen::Vector<double, 2> v5{1, 0};
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat(2, 5);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;
    mat.col(4) = v5;

    Eigen::Vector<double, 1> v12{0};
    Eigen::Vector<double, 1> v23{0.4};
    Eigen::Vector<double, 1> v33{1};
    Eigen::Matrix<double, 1, Eigen::Dynamic> mat1(1, 3);
    mat1.col(0) = v12;
    mat1.col(1) = v23;
    mat1.col(2) = v33;

    Eigen::VectorX<double> knots_vector(9);
    knots_vector << 0, 0, 0, 0, 0.5, 1, 1, 1, 1;

    Eigen::VectorX<double> new_knots_vector(6);
    new_knots_vector << 0, 0, 0, 1, 1, 1;

    // Eigen::VectorX<double> knots_vector_1(6);
    // knots_vector_1 << 1, 1, 1, 4, 4, 4;
    nurbs_curve<double, 2, false, -1, -1> curve1(knots_vector, mat);
    std::vector<Eigen::Vector2d> pointss;
    std::vector<Eigen::Vector2d> points2;
    nurbs_curve<double, 1, false, -1, -1> parameter_function(new_knots_vector, mat1);
    nurbs_curve<double, 2, false, -1, -1> new_curve;
    curve1.curve_reparameter(parameter_function, new_curve);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector2d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector2d project_point;
        new_curve.point_on_curve(0.01 * i, project_point);
        pointss.push_back(point);
        points2.push_back(project_point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }
}

void test_reparameter_3()
{
    Eigen::Vector<double, 2> v1{0, 0};
    Eigen::Vector<double, 2> v2{0.3, 0.2};
    Eigen::Vector<double, 2> v3{1.0 / 2.0, 1};
    Eigen::Vector<double, 2> v4{0.7, 0.5};
    Eigen::Vector<double, 2> v5{1, 0};
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat(2, 5);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;
    mat.col(4) = v5;

    Eigen::Vector<double, 1> v12{0};
    Eigen::Vector<double, 1> v23{0.4};
    Eigen::Vector<double, 1> v33{0.9};
    Eigen::Vector<double, 1> v34{1.0};
    Eigen::Matrix<double, 1, Eigen::Dynamic> mat1(1, 4);
    mat1.col(0) = v12;
    mat1.col(1) = v23;
    mat1.col(2) = v33;
    mat1.col(3) = v34;

    Eigen::VectorX<double> knots_vector(9);
    knots_vector << 0, 0, 0, 0, 0.5, 1, 1, 1, 1;

    Eigen::VectorX<double> new_knots_vector(7);
    new_knots_vector << 0, 0, 0, 0.5, 1, 1, 1;

    // Eigen::VectorX<double> knots_vector_1(6);
    // knots_vector_1 << 1, 1, 1, 4, 4, 4;
    nurbs_curve<double, 2, false, -1, -1> curve1(knots_vector, mat);
    std::vector<Eigen::Vector2d> pointss;
    std::vector<Eigen::Vector2d> points2;
    nurbs_curve<double, 1, false, -1, -1> parameter_function(new_knots_vector, mat1);
    nurbs_curve<double, 2, false, -1, -1> new_curve;
    curve1.curve_reparameter(parameter_function, new_curve);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector2d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector2d project_point;
        new_curve.point_on_curve(0.01 * i, project_point);
        pointss.push_back(point);
        points2.push_back(project_point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }
}


void test_reparameter_4()
{
    Eigen::Vector<double, 2> v1{0, 0};
    Eigen::Vector<double, 2> v2{1.0 / 2.0, 1};
    Eigen::Vector<double, 2> v3{1, 0};
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat(2, 3);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;

    Eigen::Vector<double, 2> v12{0, 1};
    Eigen::Vector<double, 2> v23{7.0 / 24.0, 1};
    Eigen::Vector<double, 2> v33{1, 1};
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat1(2, 3);
    mat1.col(0) = v12;
    mat1.col(1) = v23;
    mat1.col(2) = v33;

    Eigen::VectorX<double> knots_vector(6);
    knots_vector << 0, 0, 0, 1, 1, 1;

    nurbs_curve<double, 2, false, -1, -1> curve1(knots_vector, mat);
    std::vector<Eigen::Vector2d> pointss;
    std::vector<Eigen::Vector2d> points2;
    nurbs_curve<double, 1, true, -1, -1> parameter_function(knots_vector, mat1);
    nurbs_curve<double, 2, true, -1, -1> new_curve;
    curve1.curve_reparameter(parameter_function, new_curve);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector2d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector2d project_point;
        new_curve.point_on_curve(0.01 * i, project_point);
        pointss.push_back(point);
        points2.push_back(project_point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }
}

void test_reparameter_5()
{
    Eigen::Vector<double, 4> v1{0, 0, 0, 1};
    Eigen::Vector<double, 4> v2{0.3, 0.2, 0.2, 2};
    Eigen::Vector<double, 4> v3{1.0 / 2.0, 1, 0.5, 3};
    Eigen::Vector<double, 4> v4{0.7, 0.5,0.7, 2};
    Eigen::Vector<double, 4> v5{1, 0, 1, 1};
    Eigen::Matrix<double, 4, Eigen::Dynamic> mat(4, 5);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;
    mat.col(4) = v5;

    Eigen::Vector<double, 2> v12{0, 1};
    Eigen::Vector<double, 2> v23{0.4, 2};
    Eigen::Vector<double, 2> v33{0.8, 2};
    Eigen::Vector<double, 2> v34{1, 1};
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat1(2, 3);
    mat1.col(0) = v12;
    mat1.col(1) = v23;
    // mat1.col(2) = v33;
    mat1.col(2) = v34;

    Eigen::VectorX<double> knots_vector(9);
    knots_vector << 0, 0, 0, 0, 0.5, 1, 1, 1, 1;

    Eigen::VectorX<double> new_knots_vector(6);
    new_knots_vector << 0, 0, 0, 1, 1, 1;
    nurbs_curve<double, 3, true, -1, -1> curve1(knots_vector, mat);
    std::vector<Eigen::Vector3d> pointss;
    std::vector<Eigen::Vector3d> points2;
    nurbs_curve<double, 1, true, -1, -1> parameter_function(new_knots_vector, mat1);
    nurbs_curve<double, 3, true, -1, -1> new_curve;
    curve1.curve_reparameter(parameter_function, new_curve);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector3d project_point;
        new_curve.point_on_curve(0.01 * i, project_point);
        pointss.push_back(point);
        points2.push_back(project_point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_reparameter_6()
{
    Eigen::Vector<double, 3> v1{0, 0, 1};
    Eigen::Vector<double, 3> v2{0.3, 0.2, 0.5};
    Eigen::Vector<double, 3> v3{1.0 / 2.0, 1, 2};
    Eigen::Vector<double, 3> v4{0.7, 0.5, 3};
    Eigen::Vector<double, 3> v5{1, 0, 0.5};
    Eigen::Vector<double, 3> v6{1.2, -0.3, 1};
    Eigen::Matrix<double, 3, Eigen::Dynamic> mat(3, 6);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;
    mat.col(4) = v5;
    mat.col(5) = v6;


    Eigen::VectorX<double> knots_vector(10);
    knots_vector << 0, 0, 0, 0, 0.3, 0.7, 1, 1, 1, 1;


    nurbs_curve<double, 2, true, -1, -1> curve1(knots_vector, mat);
    std::vector<Eigen::Vector2d> pointss;
    std::vector<Eigen::Vector2d> points2;
    Eigen::Matrix<double, 2, 2> parameter_function;
    parameter_function << 2, -1, 3, 2;
    nurbs_curve<double, 2, true, -1, -1> new_curve;
    curve1.curve_reparameter_with_linear_function(parameter_function, new_curve);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector2d point;
        curve1.point_on_curve((double)0.01 * i, point);
        Eigen::Vector2d project_point;
        new_curve.point_on_curve(-0.5 + 0.007 * i, project_point);
        pointss.push_back(point);
        points2.push_back(project_point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }
}

void nurbs_surface_reparameter_1()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    
    Eigen::Vector<double, 2> v12{0, 1};
    // Eigen::Vector<double, 2> v23{1.2, 2};
    Eigen::Vector<double, 2> v34{10, 2};
    Eigen::Matrix<double, 2, Eigen::Dynamic> mat1(2, 2);
    mat1.col(0) = v12;
    // mat1.col(1) = v23;
    mat1.col(1) = v34;

    Eigen::VectorX<double> new_knots_vector(4);
    new_knots_vector << 0,  0, 1, 1;
    nurbs_curve<double, 1, true, -1, -1> parameter_function(new_knots_vector, mat1);

    nurbs_surface<double, 3, -1, -1, -1, -1, true> new_nurbs_surface;
    test_surface.surface_reparameter(parameter_function, ENUM_DIRECTION::U_DIRECTION, new_nurbs_surface);

    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            new_nurbs_surface.point_on_surface(i * 0.01, j * 0.03, reuslt);
            pointss.push_back(reuslt);
        }
    }

    std::vector<Eigen::Vector<double, 3>> pointss2;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss2.push_back(reuslt);
        }
    }


    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir2("view1.obj");
    std::ofstream outfile2(dir2);
    pointsCount = pointss2.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss2[index];
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void nurbs_surface_reparameter_2()
{
    Eigen::VectorX<double> u_knots_vector(7);
    u_knots_vector<< 0, 0, 0, 2.3, 5, 5, 5;
    Eigen::VectorX<double> v_knots_vector(7);
    v_knots_vector << 0, 0, 0, 1.333, 3, 3, 3;
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
    
    nurbs_surface<double, 3, -1, -1, -1, -1, true> new_nurbs_surface;
    Eigen::Matrix<double, 2, 2> parameter_function;
    parameter_function << 2, 1, 3, 2;
    test_surface.surface_reparameter_with_linear_function(parameter_function, ENUM_DIRECTION::V_DIRECTION, new_nurbs_surface);


    std::vector<Eigen::Vector<double, 3>> pointss;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            new_nurbs_surface.point_on_surface(i * 0.05, 0.5 + j * 0.00147, reuslt);
            pointss.push_back(reuslt);
        }
    }

    std::vector<Eigen::Vector<double, 3>> pointss2;
    for (int i = 0; i < 100; ++i)
    {
        for (int j = 0; j < 100; ++j)
        {
            Eigen::Vector<double, 3> reuslt;
            test_surface.point_on_surface(i * 0.05, j * 0.03, reuslt);
            pointss2.push_back(reuslt);
        }
    }


    //     // // write doc
    std::string dir1("view.obj");
    std::ofstream outfile1(dir1);
    int pointsCount = pointss.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss[index];
        outfile1 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }

    std::string dir2("view1.obj");
    std::ofstream outfile2(dir2);
    pointsCount = pointss2.size();
    for (int index = 0; index < pointsCount; ++index)
    {
        auto point = pointss2[index];
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << point[2] << std::endl;
    }
}

void test_curve_reverse_1()
{
    Eigen::Vector<double, 3> v1{0, 0, 1};
    Eigen::Vector<double, 3> v2{0.3, 0.2, 0.5};
    Eigen::Vector<double, 3> v3{1.0 / 2.0, 1, 2};
    Eigen::Vector<double, 3> v4{0.7, 0.5, 3};
    Eigen::Vector<double, 3> v5{1, 0, 0.5};
    Eigen::Vector<double, 3> v6{1.2, -0.3, 1};
    Eigen::Matrix<double, 3, Eigen::Dynamic> mat(3, 6);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;
    mat.col(4) = v5;
    mat.col(5) = v6;


    Eigen::VectorX<double> knots_vector(10);
    knots_vector << 0, 0, 0, 0, 0.3, 0.7, 1, 1, 1, 1;


    nurbs_curve<double, 2, true, -1, -1> *curve1 = new nurbs_curve<double, 2, true, -1, -1>(knots_vector, mat);
    std::vector<Eigen::Vector2d> pointss;
    std::vector<Eigen::Vector2d> points2;
    nurbs_curve<double, 2, true, -1, -1> new_curve(*curve1);
    new_curve.curve_reverse();
    for (int i = 0; i < 70; ++i)
    {
        Eigen::Vector2d point;
        curve<nurbs_curve<double, 2, true, -1, -1>> * curve2 = (static_cast<curve<nurbs_curve<double, 2, true, -1, -1>> *>(curve1));
        
        curve2->point_on_curve(0.01 * i, point);
        Eigen::Vector2d project_point;
        new_curve.point_on_curve(0.01 * i, project_point);
        pointss.push_back(point);
        points2.push_back(project_point);
    }
    delete curve1;
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : points2)
    {
        outfile2 << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }
}


void test_nurbs_to_polynomial_1()
{
    Eigen::Vector<double, 3> v1{0, 0, 1};
    Eigen::Vector<double, 3> v2{0.3, 0.2, 0.5};
    Eigen::Vector<double, 3> v3{1.0 / 2.0, 1, 2};
    Eigen::Vector<double, 3> v4{0.7, 0.5, 3};
    Eigen::Vector<double, 3> v5{1, 0, 0.5};
    Eigen::Vector<double, 3> v6{1.2, -0.3, 1};
    Eigen::Matrix<double, 3, Eigen::Dynamic> mat(3, 6);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;
    mat.col(4) = v5;
    mat.col(5) = v6;


    Eigen::VectorX<double> knots_vector(10);
    knots_vector << 0, 0, 0, 0, 0.3, 0.7, 1, 1, 1, 1;


    nurbs_curve<double, 2, true, -1, -1> *curve1 = new nurbs_curve<double, 2, true, -1, -1>(knots_vector, mat);
    std::vector<Eigen::Vector2d> pointss;
    std::vector<std::vector<Eigen::Vector2d>> points2(3);
    std::vector<polynomial_curve<double, 2, true, -1, -1>*> poly_curves(3);
    convert_nubrs_to_polynomial(*curve1, poly_curves);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector2d point, p1, p2, p3;
        (static_cast<curve<nurbs_curve<double, 2, true, -1, -1>> *>(curve1))->point_on_curve(0.01 * i, point);
        pointss.push_back(point);
        poly_curves[0]->point_on_curve(0.003 * i, p1);
        poly_curves[1]->point_on_curve(0.3 +  0.004 * i, p2);
        poly_curves[2]->point_on_curve(0.7 + 0.003 * i, p3);
        points2[0].push_back(p1);
        points2[1].push_back(p2);
        points2[2].push_back(p3);
    }
    delete curve1;
    delete poly_curves[0];
    delete poly_curves[1];
    delete poly_curves[2];

    //     // // write doc
    std::string dir("view.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }

    for (int index = 0; index < 3; ++index)
    {
        std::string dir2("view" + std::to_string(index) + ".obj");
        std::ofstream outfile2(dir2);
        for (auto point : points2[index])
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << 0 << std::endl;
        }
    }

}

void nurbs_surface_to_polynomial_1()
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
    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    Eigen::Vector<double, 2> insert_knots{1.25, 3.75};
    Eigen::Vector<double, 2> insert_knots_v{0.75, 2.25};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    Eigen::MatrixX<polynomial_surface<double, 3, -1, -1, true> *> bezier_surfaces;
    // Eigen::MatrixX<nurbs_surface<double, 3, -1, -1,-1, -1, true> *> bezier_surfaces;
    convert_nubrs_to_polynomial(test_surface, bezier_surfaces);
    // test_surface.decompose_to_bezier(bezier_surfaces);
    Eigen::Matrix<std::vector<Eigen::Vector3d>, 4, 4>  pointss;
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 4; ++i)
        {
            for (int k = 0; k <= 100; ++k)
            {
                for (int l = 0; l <= 100; ++l)
                {
                    Eigen::Vector3d point;
                    bezier_surfaces(i, j)->point_on_surface(1.25 * i + 0.0125 * k, 0.75 * j + 0.0075 * l, point);
                    pointss(i, j).push_back(point);
                }
            }
        }
    }

    std::vector<Eigen::Vector3d> pss;

    for (int k = 0; k <= 100; ++k)
    {
        for (int l = 0; l <= 100; ++l)
        {
            Eigen::Vector3d point;
            test_surface.point_on_surface(0.05 * k, 0.03 * l, point);
            pss.push_back(point);
        }
    }

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            //     // // write doc
            std::string dir("view" + std::to_string(i * 4 + j) + ".obj");
            std::ofstream outfile(dir);
            for (auto point : pointss(i, j))
            {
                outfile << "v " << point[0] << " " <<
                point[1] << " " << point[2] << std::endl;
            }
            outfile.close();
        }

    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pss)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();
}

void test_polynomial_to_nurbs_1()
{
    Eigen::Vector<double, 3> v1{0, 0, 1};
    Eigen::Vector<double, 3> v2{0.3, 0.2, 0.5};
    Eigen::Vector<double, 3> v3{1.0 / 2.0, 1, 2};
    Eigen::Vector<double, 3> v4{0.7, 0.5, 3};
    Eigen::Vector<double, 3> v5{1, 0, 0.5};
    Eigen::Vector<double, 3> v6{1.2, -0.3, 1};
    Eigen::Matrix<double, 3, Eigen::Dynamic> mat(3, 6);
    mat.col(0) = v1;
    mat.col(1) = v2;
    mat.col(2) = v3;
    mat.col(3) = v4;
    mat.col(4) = v5;
    mat.col(5) = v6;


    Eigen::VectorX<double> knots_vector(10);
    knots_vector << 0, 0, 0, 0, 0.3, 0.7, 1, 1, 1, 1;


    nurbs_curve<double, 2, true, -1, -1> *curve1 = new nurbs_curve<double, 2, true, -1, -1>(knots_vector, mat);
    std::vector<Eigen::Vector2d> pointss;
    std::vector<std::vector<Eigen::Vector2d>> points2(3);
    std::vector<polynomial_curve<double, 2, true, -1, -1>*> poly_curves(3);
    nurbs_curve<double, 2, true, -1, -1> curve11;
    nurbs_curve<double, 2, true, -1, -1> curve12;
    nurbs_curve<double, 2, true, -1, -1> curve13;
    convert_nubrs_to_polynomial(*curve1, poly_curves);
    convert_polynomial_to_nubrs(*(poly_curves[0]), Interval(0.0, 0.3), curve11);
    convert_polynomial_to_nubrs(*(poly_curves[1]), Interval(0.3, 0.7), curve12);
    convert_polynomial_to_nubrs(*(poly_curves[2]), Interval(0.7, 1.0), curve13);
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector2d point, p1, p2, p3;
        (static_cast<curve<nurbs_curve<double, 2, true, -1, -1>> *>(curve1))->point_on_curve(0.01 * i, point);
        pointss.push_back(point);
        curve11.point_on_curve(0.003 * i, p1);
        curve12.point_on_curve(0.3 +  0.004 * i, p2);
        curve13.point_on_curve(0.7 + 0.003 * i, p3);
        points2[0].push_back(p1);
        points2[1].push_back(p2);
        points2[2].push_back(p3);
    }
    delete curve1;
    delete poly_curves[0];
    delete poly_curves[1];
    delete poly_curves[2];

    //     // // write doc
    std::string dir("view.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
           point[1] << " " << 0 << std::endl;
    }

    for (int index = 0; index < 3; ++index)
    {
        std::string dir2("view" + std::to_string(index) + ".obj");
        std::ofstream outfile2(dir2);
        for (auto point : points2[index])
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << 0 << std::endl;
        }
    }

}

void polynomial_to_nurbs_surface_1()
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
    nurbs_surface<double, 3, -1, -1, -1, -1, true> test_surface(u_knots_vector, v_knots_vector, control_points);
    Eigen::Vector<double, 2> insert_knots{1.25, 3.75};
    Eigen::Vector<double, 2> insert_knots_v{0.75, 2.25};
    test_surface.refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
    test_surface.refine_knots_vector(insert_knots_v, ENUM_DIRECTION::V_DIRECTION);
    Eigen::MatrixX<polynomial_surface<double, 3, -1, -1, true> *> bezier_surfaces;
    convert_nubrs_to_polynomial(test_surface, bezier_surfaces);
    Eigen::Matrix<std::vector<Eigen::Vector3d>, 4, 4>  pointss;
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 4; ++i)
        {
            nurbs_surface<double, 3, -1, -1, -1, -1, true> nurbs_surf;
            convert_polynomial_to_nubrs(*(bezier_surfaces(i, j)), Interval(1.25 * i, 1.25 * (i + 1)), Interval(0.75 * j, 0.75 * (j + 1)), nurbs_surf);
            for (int k = 0; k <= 100; ++k)
            {
                for (int l = 0; l <= 100; ++l)
                {
                    Eigen::Vector3d point;
                    nurbs_surf.point_on_surface(1.25 * i + 0.0125 * k, 0.75 * j + 0.0075 * l, point);
                    pointss(i, j).push_back(point);
                }
            }
        }
    }

    std::vector<Eigen::Vector3d> pss;

    for (int k = 0; k <= 100; ++k)
    {
        for (int l = 0; l <= 100; ++l)
        {
            Eigen::Vector3d point;
            test_surface.point_on_surface(0.05 * k, 0.03 * l, point);
            pss.push_back(point);
        }
    }

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            //     // // write doc
            std::string dir("view" + std::to_string(i * 4 + j) + ".obj");
            std::ofstream outfile(dir);
            for (auto point : pointss(i, j))
            {
                outfile << "v " << point[0] << " " <<
                point[1] << " " << point[2] << std::endl;
            }
            outfile.close();
        }

    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pss)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();
}

void creat_nurbs_arc_1()
{
    Eigen::Vector<double, 3> center{0, 0, 0};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 1, 0};
    double radius = 10;
    double start_angles = 0;
    double end_angles = PI + 0.78;
    nurbs_curve<double, 3, true, -1, -1> nurbs;
    create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, nurbs);
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;   
        nurbs.point_on_curve(0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
}

void creat_nurbs_arc_2()
{
    Eigen::Vector<double, 3> center{0, 0, 0};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 1, 0};
    nurbs_curve<double, 3, true, -1, -1> nurbs;
    Eigen::Vector<double, 3> P0 = center + 10 * u_dir;
    Eigen::Vector<double, 3> P2 = center - 10 * v_dir;
    Eigen::Vector<double, 3> S = (5.0 * std::sqrt(2.0)) * (u_dir + v_dir);
    create_open_conic<double, 3>(P0, v_dir, P2, u_dir, S, nurbs);
    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;   
        nurbs.point_on_curve(0.01 * i, point);
        pointss.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
}

void creat_nurbs_arc_3()
{
    Eigen::Vector<double, 3> center{0, 0, 0};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 1, 0};
    double radius = 10;
    double start_angles = 0;
    double end_angles = PI + 0.78;
    nurbs_curve<double, 3, true, -1, -1> nurbs;
    create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, nurbs);
    bezier_curve<double, 3, true, -1> bezier;
    create_bezier_circle(radius, bezier);
    std::vector<Eigen::Vector3d> pointss;
    std::vector<Eigen::Vector3d> pointss2;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;   
        nurbs.point_on_curve(0.01 * i, point);
        pointss.push_back(point);
        Eigen::Vector3d point1;
        bezier.point_on_curve(0.01 * i, point1);
        pointss2.push_back(point1);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : pointss)
    {
        outfile << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (auto point : pointss2)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
}

void test_bilinear_surface_1()
{
    Eigen::Vector<double, 3> P00{10, 0, 10};
    Eigen::Vector<double, 3> P01{0, 0, 0};
    Eigen::Vector<double, 3> P10{10, 10, 0};
    Eigen::Vector<double, 3> P11{0, 10, 10};

    nurbs_surface<double, 3, -1, -1, -1, -1, false> test_surface;
    create_bilinear_nurbs_surface(P00, P01, P10, P11, test_surface);

    std::vector<Eigen::Vector3d> pss;
    for (int k = 0; k <= 100; ++k)
    {
        for (int l = 0; l <= 100; ++l)
        {
            Eigen::Vector3d point;
            test_surface.point_on_surface(0.01 * k, 0.01 * l, point);
            pss.push_back(point);
        }
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pss)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();
}

void creat_cylinder_1()
{
    Eigen::Vector<double, 3> n_dir{0, 0, 1};
    double radius = 10;

    nurbs_curve<double, 3, true, -1, -1> bezier;
    create_bezier_circle(radius, bezier);
    std::vector<Eigen::Vector3d> pointss;
    std::vector<Eigen::Vector3d> pointss2;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> surf;
    create_cylinder_surface(bezier, 10.0, n_dir, surf);
    
    std::vector<Eigen::Vector3d> pss;
    for (int k = 0; k <= 100; ++k)
    {
        for (int l = 0; l <= 100; ++l)
        {
            Eigen::Vector3d point;
            surf.point_on_surface(0.01 * k, 0.01 * l, point);
            pss.push_back(point);
        }
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pss)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();
}

void create_ruled_surface_1()
{
    Eigen::Vector<double, 3> center{0, 0, 10};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 1, 0};
    double radius = 10;

    nurbs_curve<double, 3, true, -1, -1> nurbs;
    Eigen::Vector<double, 3> P0 = center + 10 * u_dir;
    Eigen::Vector<double, 3> P2 = center - 10 * v_dir;
    Eigen::Vector<double, 3> S = (5.0 * std::sqrt(2.0)) * (u_dir + v_dir);
    create_open_conic<double, 3>(P0, v_dir, P2, u_dir, S, nurbs);


    nurbs_curve<double, 3, true, -1, -1> bezier;
    create_bezier_circle(radius, bezier);
    std::vector<Eigen::Vector3d> pointss;
    std::vector<Eigen::Vector3d> pointss2;
    nurbs_surface<double, 3, -1, -1, -1, -1, true> surf;
    create_ruled_surface(nurbs, bezier, surf);
    std::vector<Eigen::Vector3d> pss;
    for (int k = 0; k <= 100; ++k)
    {
        for (int l = 0; l <= 100; ++l)
        {
            Eigen::Vector3d point;
            surf.point_on_surface(0.01 * k, 0.01 * l, point);
            pss.push_back(point);
        }
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pss)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();
}

void create_revolved_surface_1()
{
    Eigen::Vector<double, 3> center{0, 0, 10};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 0, 1};
    double radius = 10;
    nurbs_curve<double, 3, true, -1, -1> nurbs;
    create_nurbs_circle(center, u_dir, v_dir, radius, 0.0, M_PI, nurbs);

    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;   
        nurbs.point_on_curve(0.01 * i, point);
        pointss.push_back(point);
    }

    nurbs_surface<double, 3, -1, -1, -1, -1, true> surf;
    create_revolved_surface<double>(nurbs, center, u_dir, M_PI * 2.0, surf);
    std::vector<Eigen::Vector3d> pss;
    for (int k = 0; k <= 100; ++k)
    {
        for (int l = 0; l <= 100; ++l)
        {
            Eigen::Vector3d point;
            surf.point_on_surface(0.01 * k, 0.01 * l, point);
            pss.push_back(point);
        }
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pss)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();
}

void create_scale_nurbs_surface_1()
{
    Eigen::Vector<double, 3> center{0, 0, 10};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 0, 1};
    Eigen::Vector<double, 3> z_dir{0, -1, 0};
    Eigen::Vector3d scales {2, 3, 1};
    double radius = 10;
    nurbs_curve<double, 3, true, -1, -1> nurbs;
    create_nurbs_circle(center, u_dir, v_dir, radius, 0.0, M_PI, nurbs);

    std::vector<Eigen::Vector3d> pointss;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;   
        nurbs.point_on_curve(0.01 * i, point);
        pointss.push_back(point);
    }

    nurbs_surface<double, 3, -1, -1, -1, -1, true> surf;
    create_revolved_surface<double>(nurbs, center, u_dir, M_PI * 2.0, surf);
    Eigen::Matrix<double, 3, 3> mat;
    mat << u_dir, v_dir, z_dir;
    surf.scale_surface(center, mat, scales);
    std::vector<Eigen::Vector3d> pss;
    for (int k = 0; k <= 100; ++k)
    {
        for (int l = 0; l <= 100; ++l)
        {
            Eigen::Vector3d point;
            surf.point_on_surface(0.01 * k, 0.01 * l, point);
            pss.push_back(point);
        }
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);
    for (auto point : pss)
    {
        outfile2 << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }
    outfile2.close();
}


void create_nurbs_fit_1()
{
    Eigen::Vector<double, 3> center{0, 0, 0};
    Eigen::Vector<double, 3> u_dir{1, 0, 0};
    Eigen::Vector<double, 3> v_dir{0, 1, 0};
    double radius = 10;
    double start_angles = 0;
    double end_angles = PI + 0.78;
    nurbs_curve<double, 3, true, -1, -1> nurbs;
    create_nurbs_circle(center, u_dir, v_dir, radius, start_angles, end_angles, nurbs);
    Eigen::Matrix<double, 3, Eigen::Dynamic> pointss(3, 5);
    for (int i = 0; i < 5; ++i)
    {
        Eigen::Vector3d point;   
        nurbs.point_on_curve(0.2 * i, point);
        pointss.col(i) = point;
    }


    nurbs_curve<double, 3, false, -1, -1> new_nurbs;
    global_curve_interpolate<double, 3, ENPARAMETERIEDTYPE::CHORD>(pointss, 3, new_nurbs);
    
    std::vector<Eigen::Vector3d> new_points;
    for (int i = 0; i < 100; ++i)
    {
        Eigen::Vector3d point;   
        new_nurbs.point_on_curve(0.01 * i, point);
        new_points.push_back(point);
    }
    //     // // write doc
    std::string dir("view2.obj");
    std::ofstream outfile(dir);

    for (auto point : new_points)
    {
        outfile << "v " << point[0] << " " <<
        point[1] << " " << point[2] << std::endl;
    }

    std::string dir2("view.obj");
    std::ofstream outfile2(dir2);

    for (int index = 0; index < 5; ++index)
    {
        outfile2 << "v " << pointss(0, index) << " " <<
        pointss(1, index) << " " << pointss(2, index) << std::endl;
    }

}






int main()
{

    create_nurbs_fit_1();
    return 0;
}



