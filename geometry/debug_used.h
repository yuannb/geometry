#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "nurbs_curve.h"
#include "nurbs_surface.h" 
#include "discret.h"
#include "vodes.h"
namespace tnurbs
{
    void save_box(const std::vector<Box<double, 2>>& boxes, const char* path);
    void printEigenVector(const Eigen::VectorX<double> &vec);

    void printEigenMatrix(const Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> &mat);

    void save_obj(const Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> &mat, const char *path);
    void save_obj(Eigen::Matrix<double, 3, Eigen::Dynamic> &mat, const char *path);

    void save_obj(std::vector<Eigen::Vector3d>& points, const char* path);

    void save_obj(Eigen::Matrix<double, 2, Eigen::Dynamic>& mat, const char* path);

    void printEigenMatrix(const Eigen::Matrix<double, 2, Eigen::Dynamic> &mat);

    void printEigenMatrix(const Eigen::Matrix<double, 3, Eigen::Dynamic> &mat);

    void printEigenMatrix(const Eigen::Matrix<double, 4, Eigen::Dynamic> &mat);

    void printEigenMatrix(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mat);
    void printEigenMatrix(const Eigen::SparseMatrix<double> &mat);

    ENUM_NURBS save_obj(const nurbs_curve<double, 2, false, -1, -1> &nurbs_cur, const char *path);

    ENUM_NURBS save_obj(const nurbs_curve<double, 2, true, -1, -1> &nurbs_cur, const char *path);

    ENUM_NURBS save_obj(const nurbs_curve<double, 3, false, -1, -1> &nurbs_cur, const char *path);

    ENUM_NURBS save_obj(const nurbs_curve<double, 3, true, -1, -1> &nurbs_cur, const char *path);


    ENUM_NURBS save_obj(const nurbs_surface<double, 2, -1, -1, -1, -1, false> &surf, const char *path);

    ENUM_NURBS save_obj(const nurbs_surface<double, 2, -1, -1, -1, -1, true> &surf, const char *path);

    ENUM_NURBS save_obj(const nurbs_surface<double, 3, -1, -1, -1, -1, false> &surf, const char *path);

    ENUM_NURBS save_obj(const nurbs_surface<double, 3, -1, -1, -1, -1, true> &surf, const char *path);

    ENUM_NURBS save_obj(const mesh<2> *surf, const char *path);

    ENUM_NURBS save_obj2(const nurbs_surface<double, 3, -1, -1, -1, -1, false>& surf, const char* path);

    ENUM_NURBS save_obj2(const nurbs_surface<double, 3, -1, -1, -1, -1, true>& surf, const char* path);

    ENUM_NURBS save_obj2(const nurbs_curve<double, 3, false, -1, -1>& curv, const char* path);

    ENUM_NURBS save_obj2(const nurbs_curve<double, 3, true, -1, -1>& curv, const char* path);

    ENUM_NURBS save_chat_points(const std::vector<surfs_int_points_chat<double, 3>>& points_chat, const char* path);
    
    ENUM_NURBS read_chat_points(std::vector<surfs_int_points_chat<double, 3>>& points_chat, const char* path);

    ENUM_NURBS save_chat_points(const surf_surf_int<double, 3>& points_chat, const char* path);
    
    ENUM_NURBS read_chat_points(std::vector<surfs_int_points_chat<double, 3>>& points_chat, const char* path);

    ENUM_NURBS read_chat_points(surf_surf_int<double, 3>& points_chat, const std::string& path);


} // namespace tnurbs

