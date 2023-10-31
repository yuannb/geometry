#pragma once

#include "fit_nurbs.h"

namespace tnurbs{

/// @brief 生成摆转曲面
/// @tparam T double, float
/// @tparam is_ratioal_a 轮廓线是否是有理nurbs(xz平面的曲线)
/// @tparam is_rational_b 轨道线是否是有理nurbs(xy平面的曲线)
/// @param alpha 缩放因子
/// @param profile_curve 轮廓线
/// @param trajectory_curve 轨道线
/// @return 错误码
template<typename T, bool is_ratioal_a, bool is_rational_b>
ENUM_NURBS swung_surface(T alpha, const nurbs_curve<T, 2, is_ratioal_a, -1, -1> &profile_curve, const nurbs_curve<T, 2, is_rational_b, -1, -1> &trajectory_curve,
    nurbs_surface<T, 3, -1, -1, -1, -1, is_ratioal_a || is_rational_b> &swung_nurbs)
{
    if constexpr (is_ratioal_a == false && is_rational_b == false)
    {
        Eigen::Matrix<T, 2, Eigen::Dynamic> profile_control_points = profile_curve.get_control_points();
        Eigen::Matrix<T, 2, Eigen::Dynamic> trajectory_control_points = trajectory_curve.get_control_points();
        int profile_control_points_count = profile_control_points.cols();
        int trajectory_control_points_count = trajectory_control_points.cols();
        Eigen::VectorX<Eigen::Matrix<T, 3, Eigen::Dynamic>> control_points(trajectory_control_points_count);
        for (int v_index = 0; v_index < trajectory_control_points_count; ++v_index)
        {
            control_points.resize(3, profile_control_points_count);
            for (int u_index = 0; u_index < profile_control_points_count; ++u_index)
            {
                control_points[v_index](0, u_index) = alpha * profile_control_points(0, u_index) * trajectory_control_points(0, v_index);
                control_points[v_index](1, u_index) = alpha * profile_control_points(0, u_index) * trajectory_control_points(1, v_index);
                control_points[v_index](2, u_index) = profile_control_points(1, u_index);
            }
        }
        Eigen::VectorX<T> u_knots = profile_curve.get_knots_vector();
        Eigen::VectorX<T> v_knots = trajectory_curve.get_knots_vector();
        swung_nurbs.set_control_points(control_points);
        swung_nurbs.set_uv_degree(profile_curve.get_degree(), trajectory_curve.get_degree());
        swung_nurbs.set_uv_knots(u_knots, v_knots);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    else
    {
        Eigen::Matrix<T, 3, Eigen::Dynamic> profile_control_points;// = profile_curve.get_control_points();
        Eigen::Matrix<T, 3, Eigen::Dynamic> trajectory_control_points;// = trajectory_curve.get_control_points();
        if constexpr (is_ratioal_a == true)
            profile_control_points = profile_curve.get_control_points();
        else
        {
            Eigen::Matrix<T, 2, Eigen::Dynamic> temp = profile_curve.get_control_points();
            int temp_cols = temp.cols();
            profile_control_points.resize(3, temp_cols);
            profile_control_points.block(0, 0, 2, temp_cols) = temp;
            profile_control_points.block(2, 0, 1, temp_cols).setConstant(1.0);
        }
        if constexpr (is_rational_b == true)
            trajectory_control_points = trajectory_curve.get_control_points();
        else
        {
            Eigen::Matrix<T, 2, Eigen::Dynamic> temp = trajectory_curve.get_control_points();
            int temp_cols = temp.cols();
            trajectory_control_points.resize(3, temp_cols);
            trajectory_control_points.block(0, 0, 2, temp_cols) = temp;
            trajectory_control_points.block(2, 0, 1, temp_cols).setConstant(1.0);
        }
        int profile_control_points_count = profile_control_points.cols();
        int trajectory_control_points_count = trajectory_control_points.cols();
        Eigen::VectorX<Eigen::Matrix<T, 4, Eigen::Dynamic>> control_points(trajectory_control_points_count);
        for (int v_index = 0; v_index < trajectory_control_points_count; ++v_index)
        {
            control_points[v_index].resize(4, profile_control_points_count);
            for (int u_index = 0; u_index < profile_control_points_count; ++u_index)
            {
                control_points[v_index](0, u_index) = alpha * profile_control_points(0, u_index) * trajectory_control_points(0, v_index);
                control_points[v_index](1, u_index) = alpha * profile_control_points(0, u_index) * trajectory_control_points(1, v_index);
                control_points[v_index](2, u_index) = profile_control_points(1, u_index) * trajectory_control_points(2, v_index);
                control_points[v_index](3, u_index) = profile_control_points(2, u_index) * trajectory_control_points(2, v_index);
            }
        }
        Eigen::VectorX<T> u_knots = profile_curve.get_knots_vector();
        Eigen::VectorX<T> v_knots = trajectory_curve.get_knots_vector();
        swung_nurbs.set_control_points(control_points);
        swung_nurbs.set_uv_degree(profile_curve.get_degree(), trajectory_curve.get_degree());
        swung_nurbs.set_uv_knots(u_knots, v_knots);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    return ENUM_NURBS::NURBS_ERROR;
}



}