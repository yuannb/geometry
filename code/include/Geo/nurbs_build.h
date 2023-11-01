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
/// @param swung_nurbs 摆转曲面
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

/// @brief 生成蒙皮曲面
/// @tparam T double, float
/// @tparam dim 曲线所在的欧式空间的维数
/// @tparam is_rational 曲线是否是有理(对于无理曲线和有理曲线混合作蒙皮的情况, 需要将无理的变成有理的)
/// @param section_curves 截面线
/// @param v_params v向节点矢量
/// @param skin_nurbs 蒙皮曲面
/// @return 错误码
template<typename T, int dim, bool is_rational>
ENUM_NURBS skin_surface(int v_degree, const std::vector<T> &v_params, const std::vector<nurbs_curve<T, dim, is_rational, -1, -1> *> &section_curves, nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &skin_nurbs)
{
    //获得nurbs的最高degree

    int max_degree = 0;
    for (auto &curve : section_curves)
    {
        int current_degree = curve->get_degree();
        if (max_degree < current_degree)
            max_degree = current_degree;
    }
    //拷贝一份升阶后的曲线
    std::vector<std::unique_ptr<nurbs_curve<T, dim, is_rational, -1, -1>>> section_curves_elevate;
    int curves_count = section_curves.size();
    for (int index = 0; index < curves_count; ++index)
    {
        const nurbs_curve<T, dim, is_rational, -1, -1> *curve = section_curves[index];
        std::unique_ptr<nurbs_curve<T, dim, is_rational, -1, -1>> elevate_curve = std::make_unique<nurbs_curve<T, dim, is_rational, -1, -1>>(*curve);
        elevate_curve->degree_elevate(max_degree - curve->get_degree());
        section_curves_elevate.push_back(std::move(elevate_curve));
    }

    //将参数域以第一条曲线为准对齐
    std::array<T, 2> ends_knots;
    section_curves_elevate[0]->get_ends_knots(ends_knots);
    for (int index = 1; index < curves_count; ++index)
    {
        std::array<T, 2> current_ends_knots;
        section_curves_elevate[index]->get_ends_knots(current_ends_knots);
        T alpha = (current_ends_knots[1] - current_ends_knots[0]) / (ends_knots[1] - ends_knots[0]);
        T beta = current_ends_knots[0] - alpha * ends_knots[0];
        std::unique_ptr<nurbs_curve<T, dim, is_rational, -1, -1>> reparameter_curve = std::make_unique<nurbs_curve<T, dim, is_rational, -1, -1>>();
        section_curves_elevate[index]->curve_reparameter_with_linear_function(alpha, beta, *reparameter_curve);
        std::swap(section_curves_elevate[index], reparameter_curve);
    }

    Eigen::VectorX<T> merge_knots_vector = section_curves_elevate[0]->get_knots_vector();
    for (int index = 1; index < curves_count; ++index)
    {
        std::vector<T> knots_vector1_add;
        std::vector<T> knots_vector2_add;
        Eigen::VectorX<T> current_knots = section_curves_elevate[index]->get_knots_vector();

        std::vector<T> current_merge_knots;
        merge_two_knots_vector<T>(max_degree, merge_knots_vector, current_knots, current_merge_knots, knots_vector1_add, knots_vector2_add);
        merge_knots_vector = Eigen::Map<Eigen::VectorX<T>>(current_merge_knots.data(), current_merge_knots.size());
    }

    //将每一条曲线的节点加细
    for (int index = 0; index < curves_count; ++index)
    {
        std::vector<T> knots_vector1_add;
        std::vector<T> knots_vector2_add;
        Eigen::VectorX<T> current_knots = section_curves_elevate[index]->get_knots_vector();

        std::vector<T> current_merge_knots;
        merge_two_knots_vector<T>(max_degree, merge_knots_vector, current_knots, current_merge_knots, knots_vector1_add, knots_vector2_add);
        Eigen::VectorX<T> refine_knots = Eigen::Map<Eigen::VectorX<T>>(knots_vector2_add.data(), knots_vector2_add.size());
        section_curves_elevate[index]->refine_knots_vector(refine_knots);
    }

    //将section_curves_elevate的每一条曲线的每一列控制点插值
    int u_points_count = merge_knots_vector.size() - max_degree - 1;
    int v_points_count = v_params.size();
    constexpr int points_size = is_rational ? dim + 1 : dim;
    Eigen::VectorX<Eigen::Matrix<T, points_size, Eigen::Dynamic>> control_points(v_points_count);
    for (int index = 0; index < v_points_count; ++index)
    {
        control_points[index].resize(points_size, u_points_count);
    }

    Eigen::VectorX<T> v_knots;
    for (int u_index = 0; u_index < u_points_count; ++u_index)
    {
        Eigen::Matrix<T, points_size, Eigen::Dynamic> col_control_points(points_size, curves_count);
        for (int v_index = 0; v_index < curves_count; ++v_index)
        {
            Eigen::Vector<T, points_size> p;
            section_curves_elevate[v_index]->get_control_point(u_index, p);
            col_control_points.col(v_index) = p;
        }
        nurbs_curve<T, points_size, false, -1, -1> nurbs;
        global_curve_interpolate<T, points_size>(col_control_points, v_degree, v_params, nurbs);
        Eigen::Matrix<T, points_size, Eigen::Dynamic> temp_control_points = nurbs.get_control_points();
        if (u_index == 0)
            v_knots = nurbs.get_knots_vector();
        for (int index = 0; index < v_points_count; ++index)
            control_points[index].col(u_index) = temp_control_points.col(index);
    }
    skin_nurbs.set_control_points(control_points);
    skin_nurbs.set_uv_degree(max_degree, v_degree);
    skin_nurbs.set_uv_knots(merge_knots_vector, v_knots);
    return ENUM_NURBS::NURBS_SUCCESS;
}


}