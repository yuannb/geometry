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
    /// @param v_params v向节点参数
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
            //下面的代码可以优化, 矩阵的LU分解只需要调用一次即可, 现在调用了u_points_count次
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
        if constexpr (is_rational == true)
        {
            //检查权重是否小于等于0
            for (int v_index = 0; v_index < v_points_count; ++v_index)
            {
                for (int u_index = 0; u_index < u_points_count; ++u_index)
                {
                    if (control_points[v_index](dim, u_index) <= MIN_WEIGHT<T>::value)
                        return ENUM_NURBS::NURBS_ERROR;
                }
            }
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 生成蒙皮曲面
    /// @tparam T double, float
    /// @tparam dim 曲线所在的欧式空间的维数
    /// @tparam is_rational 曲线是否是有理(对于无理曲线和有理曲线混合作蒙皮的情况, 需要将无理的变成有理的)
    /// @param section_curves 截面线
    /// @param skin_nurbs 蒙皮曲面
    /// @return 错误码
    template<typename T, int dim, bool is_rational>
    ENUM_NURBS skin_surface(int v_degree, const std::vector<nurbs_curve<T, dim, is_rational, -1, -1> *> &section_curves, nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &skin_nurbs)
    {
        int curves_count = section_curves.size();
        if (curves_count < 2)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
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

        //计算v向参数
        int u_points_count = merge_knots_vector.size() - max_degree - 1;
        int v_points_count = curves_count;
        std::vector<T> v_params(curves_count, 0.0);
        v_params[curves_count - 1] = 1.0;
        
        Eigen::ArrayX<T> d(u_points_count);
        d.setConstant(0.0);
        for (int v_index = 1; v_index < v_points_count; ++v_index)
        {
            Eigen::Matrix<T, dim, Eigen::Dynamic> mat = section_curves_elevate[v_index]->get_nonhomo_control_points() - section_curves_elevate[v_index - 1]->get_nonhomo_control_points();
            for (int u_index = 0; u_index < u_points_count; ++u_index)
                d[u_index] += mat.col(u_index).norm();
        }
        for (int v_index = 1; v_index < v_points_count - 1; ++v_index)
        {
            Eigen::Matrix<T, dim, Eigen::Dynamic> mat = section_curves_elevate[v_index]->get_nonhomo_control_points() - section_curves_elevate[v_index - 1]->get_nonhomo_control_points();
            Eigen::ArrayX<T> temp(u_points_count);
            for (int u_index = 0; u_index < u_points_count; ++u_index)
                temp[u_index] = mat.col(u_index).norm();
            temp /=  (d * (u_points_count));
            v_params[v_index] = v_params[v_index - 1];
            for (int index = 0; index < u_points_count; ++index)
                v_params[v_index] += temp[index];
        }

        //将section_curves_elevate的每一条曲线的每一列控制点插值
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
            //下面的代码可以优化, 矩阵的LU分解只需要调用一次即可, 现在调用了u_points_count次
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
        if constexpr (is_rational == true)
        {
            //检查权重是否小于等于0
            for (int v_index = 0; v_index < v_points_count; ++v_index)
            {
                for (int u_index = 0; u_index < u_points_count; ++u_index)
                {
                    if (control_points[v_index](dim, u_index) <= MIN_WEIGHT<T>::value)
                        return ENUM_NURBS::NURBS_ERROR;
                }
            }
        }
        return ENUM_NURBS::NURBS_SUCCESS;

    }



    /// @brief 带脊线的蒙皮曲面
    /// @tparam T double float
    /// @tparam dim 曲线所在的欧式空间的维数
    /// @tparam is_rational 
    /// @tparam spine_is_rational 
    /// @tparam v_degree 蒙皮曲线v方向的阶数
    /// @param section_curves 截面曲线(必须是非有理)
    /// @param spine_curve 脊线
    /// @param spline_params 脊线参数
    /// @param section_curves_params 截面线参数
    /// @param skin_nurbs 蒙皮曲面
    /// @return 错误码
    template<typename T, int dim, bool spine_is_rational, int v_degree>
    ENUM_NURBS skin_surface(const std::vector<nurbs_curve<T, dim, false, -1, -1> *> &section_curves, const nurbs_curve<T, dim, spine_is_rational, -1, -1> &spine_curve,
        const std::vector<T> &spline_params,  nurbs_surface<T, dim, -1, -1, -1, -1, false> &skin_nurbs)
    {
        //先大致判断参数合法性(暂时略过)
        static_assert(v_degree == 2 || v_degree == 3, "skin_surface : v_degree have to 2 or 3");
        int curves_count = section_curves.size();
        if (curves_count < 2)
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        //获得nurbs的最高degree
        int max_degree = 0;
        for (auto &curve : section_curves)
        {
            int current_degree = curve->get_degree();
            if (max_degree < current_degree)
                max_degree = current_degree;
        }
        //拷贝一份升阶后的曲线
        std::vector<std::unique_ptr<nurbs_curve<T, dim, false, -1, -1>>> section_curves_elevate;

        for (int index = 0; index < curves_count; ++index)
        {
            const nurbs_curve<T, dim, false, -1, -1> *curve = section_curves[index];
            std::unique_ptr<nurbs_curve<T, dim, false, -1, -1>> elevate_curve = std::make_unique<nurbs_curve<T, dim, false, -1, -1>>(*curve);
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
            std::unique_ptr<nurbs_curve<T, dim, false, -1, -1>> reparameter_curve = std::make_unique<nurbs_curve<T, dim, false, -1, -1>>();
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

        //计算v向参数
        int u_points_count = merge_knots_vector.size() - max_degree - 1;
        int v_points_count = curves_count;
        std::vector<T> v_params(curves_count, 0.0);
        v_params[curves_count - 1] = 1.0;
        
        Eigen::ArrayX<T> d(u_points_count);
        d.setConstant(0.0);
        for (int v_index = 1; v_index < v_points_count; ++v_index)
        {
            Eigen::Matrix<T, dim, Eigen::Dynamic> mat = section_curves_elevate[v_index]->get_nonhomo_control_points() - section_curves_elevate[v_index - 1]->get_nonhomo_control_points();
            for (int u_index = 0; u_index < u_points_count; ++u_index)
                d[u_index] += mat.col(u_index).norm();
        }
        for (int v_index = 1; v_index < v_points_count - 1; ++v_index)
        {
            Eigen::Matrix<T, dim, Eigen::Dynamic> mat = section_curves_elevate[v_index]->get_nonhomo_control_points() - section_curves_elevate[v_index - 1]->get_nonhomo_control_points();
            Eigen::ArrayX<T> temp(u_points_count);
            for (int u_index = 0; u_index < u_points_count; ++u_index)
                temp[u_index] = mat.col(u_index).norm();
            temp /=  (d * (u_points_count));
            v_params[v_index] = v_params[v_index - 1];
            for (int index = 0; index < u_points_count; ++index)
                v_params[v_index] += temp[index];
        }

        //计算v向切向
        Eigen::VectorX<Eigen::Vector<T, dim>> Dks(curves_count);
        Eigen::Vector2<Eigen::Vector<T, dim>> ders;
        spine_curve.derivative_on_curve<1>(spline_params[0], ders);
        Dks[0] = std::move(ders[1]);
        std::vector<T> lengths(curves_count - 1);
        Eigen::Vector<T, dim> pre_point = std::move(ders[0]);
        for (int index = 1; index < curves_count; ++index)
        {
            Eigen::Vector2<Eigen::Vector<T, dim>> temp_ders;
            spine_curve.derivative_on_curve<1>(spline_params[index], temp_ders);
            Dks[index] = std::move(temp_ders[1]);
            lengths[index - 1] = (temp_ders[0] - pre_point).norm();
            pre_point = std::move(temp_ders[0]);
        }

        //将section_curves_elevate的每一条曲线的每一列控制点插值
        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> control_points(v_points_count * 2);
        for (int index = 0; index < v_points_count * 2; ++index)
        {
            control_points[index].resize(dim, u_points_count);
        }

        Eigen::VectorX<T> v_knots;
        for (int u_index = 0; u_index < u_points_count; ++u_index)
        {
            Eigen::Matrix<T, dim, Eigen::Dynamic> col_control_points(dim, curves_count);
            Eigen::Matrix<T, dim, Eigen::Dynamic> col_control_tangent(dim, curves_count);
            for (int v_index = 0; v_index < curves_count; ++v_index)
            {

                Eigen::Vector<T, dim> p;
                section_curves_elevate[v_index]->get_control_point(u_index, p);
                if (v_index == 0)
                {
                    Eigen::Vector<T, dim> next_p;
                    section_curves_elevate[1]->get_control_point(u_index, next_p);
                    col_control_tangent.col(0) = ((p - next_p).norm() / lengths[0]) * Dks[0];
                }
                else if (v_index == curves_count - 1)
                {
                    Eigen::Vector<T, dim> pre_p;
                    section_curves_elevate[v_index - 1]->get_control_point(u_index, pre_p);
                    col_control_tangent.col(v_index) = ((p - pre_p).norm() / lengths.back()) * Dks[v_index];
                }
                else
                {
                    Eigen::Vector<T, dim> pre_p;
                    Eigen::Vector<T, dim> next_p;
                    section_curves_elevate[v_index - 1]->get_control_point(u_index, pre_p);
                    section_curves_elevate[v_index + 1]->get_control_point(u_index, next_p);
                    T temp = (p - pre_p).norm() + (p - next_p).norm();
                    col_control_tangent.col(v_index) = (temp / (lengths[v_index] + lengths[v_index - 1])) * Dks[v_index];
                }
                col_control_points.col(v_index) = std::move(p);
            }
            nurbs_curve<T, dim, false, -1, -1> nurbs;
            //下面的代码可以优化, 矩阵的LU分解只需要调用一次即可, 现在调用了u_points_count次
            global_2or3degree_hermite_curve<T, dim, v_degree>(col_control_points, col_control_tangent, v_params, nurbs);
            Eigen::Matrix<T, dim, Eigen::Dynamic> temp_control_points = nurbs.get_control_points();
            if (u_index == 0)
                v_knots = nurbs.get_knots_vector();
            for (int index = 0; index < v_points_count * 2; ++index)
                control_points[index].col(u_index) = temp_control_points.col(index);
        }
        skin_nurbs.set_control_points(control_points);
        skin_nurbs.set_uv_degree(max_degree, v_degree);
        skin_nurbs.set_uv_knots(merge_knots_vector, v_knots);
        return ENUM_NURBS::NURBS_SUCCESS;

    }


    /// @brief 通过法向投影法生成一系列标架(对于高维来说, 还需要更多的量来确定标架, 因此此方法仅用于三维情况, 故此处的维数固定为3即可)
    /// @tparam T double, float ...
    /// @tparam is_rational 是否是有理的
    /// @param params nurbs曲线上的参数点
    /// @param nurbs nurbs曲线
    /// @param frames nurbs曲线在参数点标架
    /// @return 错误码
    template<typename T, bool is_rational>
    ENUM_NURBS eval_frames_by_project(const std::vector<T> &params, const nurbs_curve<T, 3, is_rational, -1, -1> &nurbs, std::vector<frame<T, 3>> &frames)
    {
        //TODO: 首先大致检查参数合法性
        frames.clear();
        int params_count = params.size();
        frames.reserve(params_count);
        for (int index = 0; index < params_count; ++index)
        {
            // Eigen::VectorX<Eigen::Vector3<T>> ders;
            Eigen::Vector<Eigen::Vector<T, 3>, 2> ders;
            // nurbs.derivative_on_curve(params[index], 1, ders);
            nurbs.template derivative_on_curve<1>(params[index], ders);
            frame<T, 3> current_frame;
            current_frame.origin.col(0) = std::move(ders[0]);
            ders[1].normalize();
            current_frame.basis.col(0) = std::move(ders[1]);
            frames.push_back(current_frame);
        }

        //先将法向从前向后计算一遍
        Eigen::Vector3<T> bi_normal;
        Eigen::Vector3<T> tangent = frames[0].basis.col(0);
        if (std::abs(tangent[0]) > std::abs(tangent[1]) && std::abs(tangent[0] > std::abs(tangent[2])))
        {
            bi_normal = Eigen::Vector3<T>(-tangent[1], tangent[0], 0.0);
        }
        else
        {
            bi_normal = Eigen::Vector3<T>(0.0, -tangent[2], tangent[1]);
        }
        bi_normal.normalize();
        frames[0].basis.col(2) = std::move(bi_normal);

        for (int index = 1; index < params_count; ++index)
        {
            const Eigen::Vector3<T> &pre_bi_normal = frames[index - 1].basis.col(2);
            const Eigen::Vector3<T> &current_tangent = frames[index].basis.col(0);
            bi_normal =  pre_bi_normal - pre_bi_normal.dot(current_tangent) * current_tangent;
            bi_normal.normalize();
            frames[index].basis.col(2) = std::move(bi_normal);
        }

        T dis = (frames[0].origin - frames.back().origin).norm();
        if (dis < TDEFAULT_ERROR<T>::value)//闭曲线
        {
            frames.back().origin = frames[0].origin;
            frames.back().basis = frames[0].basis;

            //将法向从后往前在算一遍 
            Eigen::Vector3<T> next_bi_normal = frames[0].basis.col(2);
            for (int index = params_count - 2; index >= 0; --index)
            {
                const Eigen::Vector3<T> &current_tangent = frames[index].basis.col(0);
                bi_normal =  next_bi_normal - next_bi_normal.dot(current_tangent) * current_tangent;
                bi_normal.normalize();
                next_bi_normal = bi_normal;
                bi_normal += frames[index].basis.col(2);
                bi_normal.normalize();
                frames[index].basis.col(2) = std::move(bi_normal);
            }
        }
        //计算normal
        for (int index = 0; index < params_count; ++index)
        {
            frames[index].basis.col(1) = frames[index].basis.col(2).cross(frames[index].basis.col(0));
        }
        return ENUM_NURBS::NURBS_SUCCESS;

    }


    /// @brief 扫掠曲面(插值轨道曲线)
    /// @tparam T double float ...
    /// @tparam spine_is_rational 脊线是否是有理的
    /// @tparam scale_is_rational 放缩函数(nurbs)是否是有理的
    /// @param profile_curve 截面曲线(x轴为法向)
    /// @param spine_curve 轨道曲线
    /// @param K 截面曲线实例的最小值(K越大, 所得到的曲面约接近理论上的扫掠曲面)
    /// @param nurbs 扫掠曲面
    /// @param scale_function 放缩函数(为nullptr表示各个方向的放缩值为1, scale的中心为原点)
    /// @return 错误码
    template<typename T, bool profile_curve_is_rational, bool spine_is_rational, bool scale_is_rational = false>
    ENUM_NURBS sweep_surface(const nurbs_curve<T, 3, profile_curve_is_rational, -1, -1> &profile_curve, const nurbs_curve<T, 3, spine_is_rational, -1, -1> &spine_curve, int K, 
        nurbs_surface<T, 3, -1, -1, -1, -1, true> &nurbs, const nurbs_curve<T, 3, scale_is_rational, -1, -1> *scale_function = nullptr)
    {
        int spine_curve_degree = spine_curve.get_degree();
        Eigen::VectorX<T> spine_curve_knots = spine_curve.get_knots_vector();
        int ksv = spine_curve_knots.size();
        int nsect = K + 1;
        Eigen::VectorX<T> surface_v_knots;
        if (ksv <= nsect + spine_curve_degree) //(应该是确保后面的插值矩阵非奇异)
        {
            int m = nsect + spine_curve_degree - ksv + 1;
            std::vector<T> new_knots_vector(spine_curve_knots.data(), spine_curve_knots.data() + spine_curve_knots.rows());
            std::list<T> new_knots_list(new_knots_vector.begin(), new_knots_vector.end());

            for (int index = 0; index < m; ++index)
            {
                auto max_it = new_knots_list.begin();
                auto end_it = --new_knots_list.end();
                T max_interval_len = 0.0;
                for (auto it = new_knots_list.begin(); it != end_it; ++it)
                {
                    auto next_it = it;
                    ++next_it;
                    T len = *next_it - *it;
                    if (len > max_interval_len)
                    {
                        max_it = it;
                        max_interval_len = len;
                    }
                }
                T new_knot = max_interval_len / 2.0 + *max_it;
                new_knots_list.insert(++max_it, new_knot);
            }
            std::vector<T> V = std::vector<T>(new_knots_list.begin(), new_knots_list.end());
            surface_v_knots = Eigen::Map<Eigen::VectorX<T>>(V.data(), nsect + spine_curve_degree + 1);
        }
        else
        {
            surface_v_knots = std::move(spine_curve_knots);
            if (ksv > nsect + spine_curve_degree + 1)
                nsect = ksv - spine_curve_degree - 1;
        }

        std::array<T, 2> end_knots;
        spine_curve.get_ends_knots(end_knots);
        std::vector<T> params(nsect);
        params[0] = end_knots[0];
        params[nsect - 1] = end_knots[1];
        T temp_sum = end_knots[0] * spine_curve_degree;
        for (int k = 1; k < nsect - 1; ++k)
        {
            temp_sum = temp_sum - surface_v_knots[k] + surface_v_knots[k + spine_curve_degree];
            params[k] = temp_sum / spine_curve_degree;
        }
        std::vector<frame<T, 3>> frames;
        eval_frames_by_project<T, spine_is_rational>(params, spine_curve, frames);
        int profile_curve_control_points_count = profile_curve.get_control_points_count();
        std::vector<Eigen::Matrix<T, 4, Eigen::Dynamic>> temp_control_points;
        temp_control_points.resize(profile_curve_control_points_count);
        for (int index = 0; index < profile_curve_control_points_count; ++index)
        {
            temp_control_points[index].resize(4, nsect);
        }

        for (int k = 0; k < nsect; ++k)
        {
            Eigen::Matrix3<T> scale_mat;
            scale_mat.setConstant(0.0);
            if (scale_function != nullptr)
            {
                Eigen::Vector3<T> scales;
                scale_function->point_on_curve(params[k], scales);
                scale_mat(0, 0) = scales[0];
                scale_mat(1, 1) = scales[1];
                scale_mat(2, 2) = scales[2];
            }
            else
            {
                scale_mat(0, 0) = 1.0;
                scale_mat(1, 1) = 1.0;
                scale_mat(2, 2) = 1.0;
            }
            T sw = 1.0;
            if constexpr (spine_is_rational == true)
            {
                spine_curve.weight_on_curve(params[k], sw);
            }
            for (int index = 0; index < profile_curve_control_points_count; ++index)
            {
                Eigen::Vector<T, 3> p;
                profile_curve.get_control_point(index, p);
                p = frames[k].basis * scale_mat * p;
                p += frames[k].origin;
                Eigen::Vector<T, 4> homo_p;
                homo_p.block(0, 0, 3, 1) = p;
                homo_p[3] = 1.0;
                T w;
                profile_curve.get_weight(index, w);
                homo_p *= (w * sw);
                temp_control_points[index].col(k) = homo_p;
            }

        }

        Eigen::VectorX<Eigen::Matrix4X<T>> new_control_points;
        new_control_points.resize(nsect);
        for (int index = 0; index < nsect; ++index)
        {
            new_control_points[index].resize(4, profile_curve_control_points_count);
        }
        for (int index = 0; index < profile_curve_control_points_count; ++index)
        {
            nurbs_curve<T, 4, false, -1, -1> temp_nurbs;
            global_curve_interpolate<T, 4>(temp_control_points[index], spine_curve_degree, params, surface_v_knots, temp_nurbs);
            for (int i = 0; i < nsect; ++i)
            {
                Eigen::Vector4<T> p;
                temp_nurbs.get_control_point(i, p);
                new_control_points[i].col(index) = p;
            }
        }
        nurbs.set_control_points(new_control_points);
        nurbs.set_uv_knots(profile_curve.get_knots_vector(), surface_v_knots);
        nurbs.set_uv_degree(profile_curve.get_degree(), spine_curve_degree);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    /// @brief 扫掠曲面(不插值轨道曲线)(TODO: 对于无理的proflile, 使其生成的nurbs surface为无理的)
    /// @tparam T double float ...
    /// @tparam spine_is_rational 脊线是否是有理的
    /// @tparam scale_is_rational 放缩函数(nurbs)是否是有理的
    /// @param profile_curve 截面曲线(x轴为法向)
    /// @param spine_curve 轨道曲线
    /// @param v_degree 扫掠曲面v向阶数
    /// @param K 截面曲线实例的最小值(K越大, 所得到的曲面约接近理论上的扫掠曲面)
    /// @param nurbs 扫掠曲面
    /// @param scale_function 放缩函数(为nullptr表示各个方向的放缩值为1, scale的中心为原点)
    /// @return 错误码
    template<typename T, bool profile_curve_is_rational, bool spine_is_rational, bool scale_is_rational = false>
    ENUM_NURBS sweep_surface(const nurbs_curve<T, 3, profile_curve_is_rational, -1, -1> &profile_curve, const nurbs_curve<T, 3, spine_is_rational, -1, -1> &spine_curve, int v_degree, int K, 
        nurbs_surface<T, 3, -1, -1, -1, -1, true> &nurbs, const nurbs_curve<T, 3, scale_is_rational, -1, -1> *scale_function = nullptr)
    {
        //简单等分参数域
        std::array<T, 2> end_knots;
        spine_curve.get_ends_knots(end_knots);
        std::vector<T> params(K + 1);

        T step = (end_knots[1] - end_knots[0]) / K;
        params[0] = end_knots[0];
        for (int k = 1; k < K; ++k)
        {
            params[k] = step * static_cast<T> (k);
        }
        params[K] = end_knots[1];
        
        std::vector<frame<T, 3>> frames;
        eval_frames_by_project<T, spine_is_rational>(params, spine_curve, frames);
        int profile_curve_control_points_count = profile_curve.get_control_points_count();
        std::vector<Eigen::Matrix<T, 4, Eigen::Dynamic>> temp_control_points;
        temp_control_points.resize(profile_curve_control_points_count);
        for (int index = 0; index < profile_curve_control_points_count; ++index)
        {
            temp_control_points[index].resize(4, K + 1);
        }

        for (int k = 0; k <= K; ++k)
        {
            Eigen::Matrix3<T> scale_mat;
            scale_mat.setConstant(0.0);
            if (scale_function != nullptr)
            {
                Eigen::Vector3<T> scales;
                scale_function->point_on_curve(params[k], scales);
                scale_mat(0, 0) = scales[0];
                scale_mat(1, 1) = scales[1];
                scale_mat(2, 2) = scales[2];
            }
            else
            {
                scale_mat(0, 0) = 1.0;
                scale_mat(1, 1) = 1.0;
                scale_mat(2, 2) = 1.0;
            }

            for (int index = 0; index < profile_curve_control_points_count; ++index)
            {
                Eigen::Vector<T, 3> p;
                profile_curve.get_control_point(index, p);
                p = frames[k].basis * scale_mat * p;
                p += frames[k].origin;
                Eigen::Vector<T, 4> homo_p;
                homo_p.block(0, 0, 3, 1) = p;
                homo_p[3] = 1.0;
                T w;
                profile_curve.get_weight(index, w);
                homo_p *= w;
                temp_control_points[index].col(k) = homo_p;
            }

        }

        Eigen::VectorX<Eigen::Matrix4X<T>> new_control_points;
        new_control_points.resize(K + 1);
        Eigen::VectorX<T> surface_v_knots;
        for (int index = 0; index <= K; ++index)
        {
            new_control_points[index].resize(4, profile_curve_control_points_count);
        }
        for (int index = 0; index < profile_curve_control_points_count; ++index)
        {
            nurbs_curve<T, 4, false, -1, -1> temp_nurbs;
            global_curve_interpolate<T, 4>(temp_control_points[index], v_degree, params, temp_nurbs);
            if (index == 0)
                surface_v_knots = temp_nurbs.get_knots_vector();
            for (int i = 0; i <= K; ++i)
            {
                Eigen::Vector4<T> p;
                temp_nurbs.get_control_point(i, p);
                new_control_points[i].col(index) = p;
            }
        }
        nurbs.set_control_points(new_control_points);
        nurbs.set_uv_knots(profile_curve.get_knots_vector(), surface_v_knots);
        nurbs.set_uv_degree(profile_curve.get_degree(), v_degree);
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief 戈登曲面(各曲线阶数相等)
    /// @tparam T double, float...
    /// @tparam dim 曲线所在欧式空间的维数
    /// @param CK u向曲线簇(所有曲线的区间相同, 阶数相同)
    /// @param Cl v向曲线簇(所有曲线的区间相同, 阶数相同)
    /// @param u_params u向参数
    /// @param v_params v向参数
    /// @param pt 插值曲面u向的最低阶数
    /// @param qt 插值曲面u向的最低阶数
    /// @param gordon_nurbs 戈登曲面
    /// @return 错误码
    template<typename T, int dim>
    ENUM_NURBS gordon_surface(const std::vector<nurbs_curve<T, dim, false, -1, -1> *> &CK , const std::vector<nurbs_curve<T, dim, false, -1, -1> *> &Cl, const std::vector<T> &u_params,
        const std::vector<T> &v_params, int pt, int qt, nurbs_surface<T, dim, -1, -1, -1, -1, false> &gordon_nurbs)
    {
        //不检测曲线是否在规定的点相交, 假设相交
        int ck_degree = CK[0]->get_degree();
        int ck_count = CK.size();
        for (int index = 1; index < ck_count; ++index)
        {
            int temp_degree = CK[index]->get_degree();
            if (ck_degree != temp_degree)
                return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        }
        if (ck_count != static_cast<int> (v_params.size()))
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        int cl_degree = Cl[0]->get_degree();
        int cl_count = Cl.size();
        for (int index = 1; index < cl_count; ++index)
        {
            int temp_degree = Cl[index]->get_degree();
            if (cl_degree != temp_degree)
                return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        }
        if (cl_count != static_cast<int> (u_params.size()))
            return ENUM_NURBS::NURBS_PARAM_IS_INVALID;
        nurbs_surface<T, dim, -1, -1, -1, -1, false> CK_surface;
        skin_surface<T, dim>(qt, v_params, CK, CK_surface);
        nurbs_surface<T, dim, -1, -1, -1, -1, false> Cl_surface;
        skin_surface<T, dim>(pt, u_params, Cl, Cl_surface);
        Cl_surface.reverse_uv();
        int u_degree = std::max(pt, ck_degree);
        int v_degree = std::max(qt, cl_degree);
        nurbs_surface<T, dim, -1, -1, -1, -1, false> interpolate_surface;

        //交点
        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> points(ck_count);
        for (int v_index = 0; v_index < ck_count; ++v_index)
        {
            points[v_index].resize(dim, cl_count);
            for (int u_index = 0; u_index < cl_count; ++u_index)
            {
                Eigen::Vector<T, dim> p;
                CK[v_index]->point_on_curve(u_params[u_index], p);
                points[v_index].col(u_index) = std::move(p);
            }
        }
        global_surface_interpolate<T, dim>(points, pt, qt, u_params, v_params, interpolate_surface);
        
        //将interpolate_surface, CK_surface, Cl_surface的区间按照interpolate_surface的区间对齐
        std::array<T, 2> u_knots_end, v_knots_end;
        interpolate_surface.get_uv_knots_end(u_knots_end, v_knots_end);

        std::array<T, 2> temp_u_knots_end, temp_v_knots_end;
        CK_surface.get_uv_knots_end(temp_u_knots_end, temp_v_knots_end);
        T alpha = (temp_u_knots_end[1] - temp_u_knots_end[0]) / (u_knots_end[1] - u_knots_end[0]);
        T beta = temp_u_knots_end[0] - alpha * u_knots_end[0];
        nurbs_surface<T, dim, -1, -1, -1, -1, false>  temp_surface;
        CK_surface.surface_reparameter_with_linear_function(alpha, beta, ENUM_DIRECTION::U_DIRECTION, temp_surface);

        alpha = (temp_v_knots_end[1] - temp_v_knots_end[0]) / (v_knots_end[1] - v_knots_end[0]);
        beta = temp_v_knots_end[0] - alpha * v_knots_end[0];
        temp_surface.surface_reparameter_with_linear_function(alpha, beta, ENUM_DIRECTION::V_DIRECTION, CK_surface);

        Cl_surface.get_uv_knots_end(temp_u_knots_end, temp_v_knots_end);
        alpha = (temp_u_knots_end[1] - temp_u_knots_end[0]) / (u_knots_end[1] - u_knots_end[0]);
        beta = temp_u_knots_end[0] - alpha * u_knots_end[0];
        nurbs_surface<T, dim, -1, -1, -1, -1, false>  new_CK_surface;
        Cl_surface.surface_reparameter_with_linear_function(alpha, beta, ENUM_DIRECTION::U_DIRECTION, temp_surface);

        alpha = (temp_v_knots_end[1] - temp_v_knots_end[0]) / (v_knots_end[1] - v_knots_end[0]);
        beta = temp_v_knots_end[0] - alpha * v_knots_end[0];
        temp_surface.surface_reparameter_with_linear_function(alpha, beta, ENUM_DIRECTION::V_DIRECTION, Cl_surface);


        //升阶
        CK_surface.degree_elevate(u_degree - ck_degree, ENUM_DIRECTION::U_DIRECTION);
        CK_surface.degree_elevate(v_degree - qt, ENUM_DIRECTION::V_DIRECTION);

        Cl_surface.degree_elevate(u_degree - pt, ENUM_DIRECTION::U_DIRECTION);
        Cl_surface.degree_elevate(v_degree - cl_degree, ENUM_DIRECTION::V_DIRECTION);

        interpolate_surface.degree_elevate(u_degree - pt, ENUM_DIRECTION::U_DIRECTION);
        interpolate_surface.degree_elevate(v_degree - qt, ENUM_DIRECTION::V_DIRECTION);

        //合并节点矢量
        Eigen::VectorX<T> CK_u_knots_vector = CK_surface.get_u_knots();
        Eigen::VectorX<T> Cl_u_knots_vector = Cl_surface.get_u_knots();
        Eigen::VectorX<T> int_u_knots_vector = interpolate_surface.get_u_knots();
        std::vector<T> knots_vector1_add;
        std::vector<T> knots_vector2_add;
        std::vector<T> u_merged_knots;
        merge_two_knots_vector<T>(u_degree, CK_u_knots_vector, Cl_u_knots_vector, u_merged_knots, knots_vector1_add, knots_vector2_add);
        Eigen::VectorX<T> temp = Eigen::Map<Eigen::VectorX<T>>(u_merged_knots.data(), u_merged_knots.size());
        merge_two_knots_vector<T>(u_degree, temp, int_u_knots_vector, u_merged_knots, knots_vector1_add, knots_vector2_add);

        Eigen::VectorX<T> CK_v_knots_vector = CK_surface.get_v_knots();
        Eigen::VectorX<T> Cl_v_knots_vector = Cl_surface.get_v_knots();
        Eigen::VectorX<T> int_v_knots_vector = interpolate_surface.get_v_knots();
        std::vector<T> v_merged_knots;
        merge_two_knots_vector<T>(v_degree, CK_v_knots_vector, Cl_v_knots_vector, v_merged_knots, knots_vector1_add, knots_vector2_add);
        Eigen::VectorX<T> temp2 = Eigen::Map<Eigen::VectorX<T>>(v_merged_knots.data(), v_merged_knots.size());
        merge_two_knots_vector<T>(v_degree, temp2, int_v_knots_vector, v_merged_knots, knots_vector1_add, knots_vector2_add);

        //加细
        std::vector<T> temp_vec;
        Eigen::VectorX<T> u_merged_knots_vector = Eigen::Map<Eigen::VectorX<T>>(u_merged_knots.data(), u_merged_knots.size());
        merge_two_knots_vector<T>(u_degree, CK_u_knots_vector, u_merged_knots_vector, temp_vec, knots_vector1_add, knots_vector2_add);
        Eigen::VectorX<T> refine_knots = Eigen::Map<Eigen::VectorX<T>>(knots_vector1_add.data(), knots_vector1_add.size());
        CK_surface.refine_knots_vector(refine_knots, ENUM_DIRECTION::U_DIRECTION);
        merge_two_knots_vector<T>(u_degree, Cl_u_knots_vector, u_merged_knots_vector, temp_vec, knots_vector1_add, knots_vector2_add);
        refine_knots = Eigen::Map<Eigen::VectorX<T>>(knots_vector1_add.data(), knots_vector1_add.size());
        Cl_surface.refine_knots_vector(refine_knots, ENUM_DIRECTION::U_DIRECTION);
        merge_two_knots_vector<T>(u_degree, int_u_knots_vector, u_merged_knots_vector, temp_vec, knots_vector1_add, knots_vector2_add);
        refine_knots = Eigen::Map<Eigen::VectorX<T>>(knots_vector1_add.data(), knots_vector1_add.size());
        interpolate_surface.refine_knots_vector(refine_knots, ENUM_DIRECTION::U_DIRECTION);

        Eigen::VectorX<T> v_merged_knots_vector = Eigen::Map<Eigen::VectorX<T>>(v_merged_knots.data(), v_merged_knots.size());
        merge_two_knots_vector<T>(v_degree, CK_v_knots_vector, v_merged_knots_vector, temp_vec, knots_vector1_add, knots_vector2_add);
        refine_knots = Eigen::Map<Eigen::VectorX<T>>(knots_vector1_add.data(), knots_vector1_add.size());
        CK_surface.refine_knots_vector(refine_knots, ENUM_DIRECTION::V_DIRECTION);
        merge_two_knots_vector<T>(v_degree, Cl_v_knots_vector, v_merged_knots_vector, temp_vec, knots_vector1_add, knots_vector2_add);
        refine_knots = Eigen::Map<Eigen::VectorX<T>>(knots_vector1_add.data(), knots_vector1_add.size());
        Cl_surface.refine_knots_vector(refine_knots, ENUM_DIRECTION::V_DIRECTION);
        merge_two_knots_vector<T>(v_degree, int_v_knots_vector, v_merged_knots_vector, temp_vec, knots_vector1_add, knots_vector2_add);
        refine_knots = Eigen::Map<Eigen::VectorX<T>>(knots_vector1_add.data(), knots_vector1_add.size());
        interpolate_surface.refine_knots_vector(refine_knots, ENUM_DIRECTION::V_DIRECTION);


        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> new_control_points = interpolate_surface.get_control_points();
        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> points1 = Cl_surface.get_control_points();
        Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> points2 = CK_surface.get_control_points();
        int v_count = new_control_points.rows();
        for (int v_index = 0; v_index < v_count; ++v_index)
        {
            new_control_points[v_index] = points1[v_index] + points2[v_index] - new_control_points[v_index];
        }
        gordon_nurbs.set_control_points(new_control_points);
        gordon_nurbs.set_uv_degree(u_degree, v_degree);
        gordon_nurbs.set_uv_knots(u_merged_knots_vector, v_merged_knots_vector);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim, bool is_rational>
    ENUM_NURBS coons_surface(const std::array<nurbs_curve<T, dim, is_rational, -1, -1> *, 2> CK, const std::array<nurbs_curve<T, dim, is_rational, -1, -1> *, 2> Cl, nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &coons_nurbs)
    {
        if constexpr (is_rational == false)
        {
            std::vector<nurbs_curve<T, dim, false, -1, -1> *> row_curve{CK[0], CK[1]};
            std::vector<nurbs_curve<T, dim, false, -1, -1> *> col_curve{Cl[0], Cl[1]};
            std::array<T, 2> row_ends_knots;
            CK[0]->get_ends_knots(row_ends_knots);
            std::vector<T> u_params{row_ends_knots[0], row_ends_knots[1]}; 
            Cl[0]->get_ends_knots(row_ends_knots);
            std::vector<T> v_params{row_ends_knots[0], row_ends_knots[1]};
            return gordon_surface<T, dim>(row_curve, col_curve, u_params, v_params, 1, 1, coons_nurbs);
        }
        else
        {
            T u0_low_weight, u1_low_weight, u0_high_weight, u1_high_weight;
            CK[0]->get_weight(0, u0_low_weight);
            CK[1]->get_weight(0, u1_low_weight);
            int CK0_control_points_count = CK[0]->get_control_points_count();
            int CK1_control_points_count = CK[0]->get_control_points_count();
            CK[0]->get_weight(CK0_control_points_count - 1, u0_high_weight);
            CK[1]->get_weight(CK1_control_points_count - 1, u1_high_weight);
            std::array<T, 2> u_ends_knots, v_ends_knots;
            CK[0]->get_ends_knots(u_ends_knots);
            Cl[0]->get_ends_knots(v_ends_knots);
            std::vector<T> u_params(u_ends_knots.begin(), u_ends_knots.end());
            std::vector<T> v_params(v_ends_knots.begin(), v_ends_knots.end());
            int v_degree = Cl[0]->get_degree();

            T v0_low_weight, v1_low_weight, v0_high_weight, v1_high_weight;
            Cl[0]->get_weight(0, v0_low_weight);
            Cl[1]->get_weight(0, v1_low_weight);
            int Cl0_control_points_count = Cl[0]->get_control_points_count();
            int Cl1_control_points_count = Cl[0]->get_control_points_count();
            Cl[0]->get_weight(Cl0_control_points_count - 1, v0_high_weight);
            Cl[1]->get_weight(Cl1_control_points_count - 1, v1_high_weight);
            T alpha, beta, gamma, delta;
            Eigen::Matrix2<T> mat;
            nurbs_curve<T, dim, true, -1, -1> new_col_curve0;
            if (std::abs(v0_low_weight / u0_low_weight - v0_high_weight / u1_low_weight) < MIN_WEIGHT<T>::value)
            {
                T x = u0_low_weight / v0_low_weight;
                auto points = Cl[0]->get_control_points();
                points *= x;
                new_col_curve0.set_control_points(points);
                new_col_curve0.set_knots_vector(Cl[0]->get_knots_vector());
                new_col_curve0.set_degree(Cl[0]->get_degree());
            }
            else
            {
                reparameter_map<T>(v_degree, v_ends_knots[0], v_ends_knots[1], v0_low_weight, v0_high_weight, u0_low_weight, u1_low_weight, alpha, beta, gamma, delta);   
                mat(0, 0) = alpha;
                mat(0, 1) = beta;
                mat(1, 0) = gamma;
                mat(1, 1) = delta;
                Cl[0]->curve_reparameter_with_linear_function(mat, new_col_curve0);
            }
            nurbs_curve<T, dim, true, -1, -1> new_col_curve1;
            if (std::abs(v1_low_weight / u0_high_weight - v1_high_weight / u1_high_weight) < MIN_WEIGHT<T>::value)
            {
                T x = u0_high_weight / v1_low_weight;
                auto points = Cl[1]->get_control_points();
                points *= x;
                new_col_curve1.set_control_points(points);
                new_col_curve1.set_knots_vector(Cl[1]->get_knots_vector());
                new_col_curve1.set_degree(Cl[1]->get_degree());
            }
            else
            {
                reparameter_map<T>(v_degree, v_ends_knots[0], v_ends_knots[1], v1_low_weight, v1_high_weight, u0_high_weight, u1_high_weight, alpha, beta, gamma, delta);
                mat(0, 0) = alpha;
                mat(0, 1) = beta;
                mat(1, 0) = gamma;
                mat(1, 1) = delta;
                
                Cl[1]->curve_reparameter_with_linear_function(mat, new_col_curve1);
            }
            
            nurbs_curve<T, dim + 1, false, -1, -1> c0;
            c0.set_control_points(CK[0]->get_control_points());
            c0.set_knots_vector(CK[0]->get_knots_vector());
            c0.set_degree(CK[0]->get_degree());

            nurbs_curve<T, dim + 1, false, -1, -1> c1;
            c1.set_control_points(CK[1]->get_control_points());
            c1.set_knots_vector(CK[1]->get_knots_vector());
            c1.set_degree(CK[1]->get_degree());
            std::vector<nurbs_curve<T, dim + 1, false, -1, -1>*> new_cks{&c0, &c1};

            nurbs_curve<T, dim + 1, false, -1, -1> r0;
            r0.set_control_points(new_col_curve0.get_control_points());
            r0.set_knots_vector(new_col_curve0.get_knots_vector());
            r0.set_degree(new_col_curve0.get_degree());

            nurbs_curve<T, dim + 1, false, -1, -1> r1;
            r1.set_control_points(new_col_curve1.get_control_points());
            r1.set_knots_vector(new_col_curve1.get_knots_vector());
            r1.set_degree(new_col_curve1.get_degree());
            std::vector<nurbs_curve<T, dim + 1, false, -1, -1>*> new_cls{&r0, &r1};
            nurbs_surface<T, dim + 1, -1, -1, -1, -1, false> temp_surface;
            gordon_surface<T, dim + 1>(new_cks, new_cls, u_params, v_params, 1, 1, temp_surface);
            coons_nurbs.set_control_points(temp_surface.get_control_points());
            coons_nurbs.set_uv_degree(temp_surface.get_u_degree(), temp_surface.get_v_degree());
            coons_nurbs.set_uv_knots(temp_surface.get_u_knots(), temp_surface.get_v_knots());
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        return ENUM_NURBS::NURBS_ERROR;
    }


}