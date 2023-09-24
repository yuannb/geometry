#pragma once

#include "nurbs_curve.h"
#include "nurbs_surface.h"
#include <numbers>

//start_angle < end_angle : 逆时针圆弧;  start_angle > end_angle : 顺时针的圆弧
template<typename T, int dim>
ENUM_NURBS create_nurbs_circle(const Eigen::Vector<T, dim> &center, const Eigen::Vector<T, dim> &u_dir, const Eigen::Vector<T, dim> &v_dir,
    T radius, T start_angle, T end_angle, nurbs_curve<T, dim, true, -1, -1> &nurbs)
{
    // if(end_angle < start_angle)
    //     return ENUM_NURBS::NURBS_ERROR;
    T theta = end_angle - start_angle;
    if (theta > M_PI * 2 + DEFAULT_ERROR)
        return ENUM_NURBS::NURBS_ERROR;
    int narcs = std::ceil(theta / M_PI_2);
    T delta_theta = theta / static_cast<T>(narcs);
    T w1 = std::cos(delta_theta / 2.0);
    Eigen::Vector<T, dim> P0 = center + radius * std::cos(start_angle) * u_dir + radius * std::sin(start_angle) * v_dir;
    Eigen::Vector<T, dim> T0 = -std::sin(start_angle) * u_dir + std::cos(start_angle) * v_dir;
    Eigen::Matrix<T, dim + 1, Eigen::Dynamic> control_points(dim + 1, 2 * narcs + 1);
    control_points.col(0).template block<dim, 1>(0, 0) = P0;
    control_points(dim, 0) = 1.0;
    int index = 0;
    T angle = start_angle;
    for (int i = 1; i <= narcs; ++i)
    {
        angle += delta_theta;
        Eigen::Vector<T, dim> P2 = center + radius * std::cos(angle) * u_dir + radius * std::sin(angle) * v_dir;
        control_points.template block<dim, 1>(0, index + 2) = P2;
        control_points(dim, index + 2) = 1.0;
        Eigen::Vector<T, dim> T2 = -std::sin(angle) * u_dir + std::cos(angle) * v_dir;

        Eigen::Matrix<T, dim, 2> mat;
        mat << T0, -T2;
        Eigen::Vector<T, dim> vec = P2 - P0;
        Eigen::Matrix<T, dim, 3> externMat;
        externMat << T0, -T2, vec;
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> externMatSvd(externMat);
        int rankMat = matSvd.rank();
        int rankExternMat = externMatSvd.rank();
        if (rankMat != 2 && rankExternMat != 2)
        {
            if (radius < DEFAULT_ERROR)
            {
                control_points(dim, index + 1) = w1;
                control_points.template block<dim, 1>(0, index + 1) = w1 * center;
                index += 2;
                if (i < narcs)
                {
                    P0 = P2;
                    T0 = T2;
                }
                continue;
            }
            else
                return ENUM_NURBS::NURBS_ERROR;
        }
        Eigen::Vector2<T> v = matSvd.solve(vec);
        Eigen::Vector<T, dim> middle_point = ((v[0] * T0 + P0) + (v[1] * T2 + P2)) / 2.0;
        control_points(dim, index + 1) = w1;
        control_points.template block<dim, 1>(0, index + 1) = w1 * middle_point;
        index += 2;
        if (i < narcs)
        {
            P0 = P2;
            T0 = T2;
        }
    }
    int len = (narcs - 1) * 2;
    Eigen::VectorX<T> knots_vector(6 + len);
    knots_vector.template block<3, 1>(0, 0).setConstant(0.0);
    knots_vector.template block<3, 1>(3 + len, 0).setConstant(1.0);
    if (narcs == 2)
    {
        knots_vector[3] = 0.5;
        knots_vector[4] = 0.5;
    }
    else if (narcs == 3)
    {
        knots_vector[3] = knots_vector[4] = 1.0 / 3.0;
        knots_vector[5] = knots_vector[6] = 2.0 * knots_vector[3];
    }
    else if (narcs == 4)
    {
        knots_vector[3] = knots_vector[4] = 0.25;
        knots_vector[5] = knots_vector[6] = 0.5;
        knots_vector[7] = knots_vector[8] = 0.75;
    }
    nurbs.set_control_points(control_points);
    nurbs.set_knots_vector(knots_vector);
    nurbs.set_degree(2);
    return ENUM_NURBS::NURBS_SUCCESS;
}

template<typename T, int dim>
ENUM_NURBS create_one_arc(const Eigen::Vector<T, dim> &P0, const Eigen::Vector<T, dim> &T0, const Eigen::Vector<T, dim> &P2, const Eigen::Vector<T, dim> &T2,
    const Eigen::Vector<T, dim> &P, Eigen::Vector<T, dim> &P1, T &w1)
{
    if ((P2 - P0).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR)
        return ENUM_NURBS::NURBS_ERROR;
    if ((P - P0).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR)
        return ENUM_NURBS::NURBS_ERROR;
    if ((P2 - P).squaredNorm() < DEFAULT_ERROR * DEFAULT_ERROR)
        return ENUM_NURBS::NURBS_ERROR;
    
    T T0_T2_length = T0.squaredNorm() * T2.squaredNorm();
    T T0_T2_cross_product_length = (T0.cross(T2)).squaredNorm();
    if (T0_T2_cross_product_length < T0_T2_length * DEFAULT_ERROR * DEFAULT_ERROR)
    {
        //起始点切向和终止点切向平行 ： 无穷远控制点
        w1 = 0.0;
        Eigen::Vector<T, dim> direction_vec = P2 - P0;
        Eigen::Matrix<T, dim, 2> mat;
        mat << direction_vec, -T0;
        Eigen::Vector<T, dim> vec = P - P0;
        Eigen::Matrix<T, dim, 3> externMat;
        externMat << direction_vec, -T0, vec;
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> externMatSvd(externMat);
        int rankMat = matSvd.rank();
        int rankExternMat = externMatSvd.rank();
        if (rankMat != 2 && rankExternMat != 2)
        {
            return ENUM_NURBS::NURBS_ERROR;
        }
        Eigen::Vector2<T> v = matSvd.solve(vec);
        if (v[0] > 1.0 || v[0] < 0.0)
            return ENUM_NURBS::NURBS_ERROR;
        T a = std::sqrt(v[0] / (1.0 - v[0]));
        T u = a / (1.0 + a);
        T b = 2.0 * u * (1.0 - u);
        b = (1.0 - b) / b;
        P1 =-v[1] * b * T0;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    // else
    Eigen::Matrix<T, dim, 2> mat;
    mat << T0, -T2;
    Eigen::Vector<T, dim> vec = P2 - P0;
    Eigen::Matrix<T, dim, 3> externMat;
    externMat << T0, -T2, vec;
    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> externMatSvd(externMat);
    int rankMat = matSvd.rank();
    int rankExternMat = externMatSvd.rank();
    if (rankMat != 2 && rankExternMat != 2)
    {
        return ENUM_NURBS::NURBS_ERROR;
    }
    Eigen::Vector2<T> v = matSvd.solve(vec);
    
    P1 = ((v[0] * T0 + P0) + (v[1] * T2 + P2)) / 2.0;
    Eigen::Vector<T, dim> v1p = P - P1;
    Eigen::Matrix<T, dim, 2> mat2;
    mat2 << v1p, -vec;
    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd2(mat2);
    v = matSvd2.solve(P0 - P1);
    if (v[1] > 1.0 || v[1] < 0.0)
        return ENUM_NURBS::NURBS_ERROR;
    T a = std::sqrt(v[1] / (1.0 - v[1]));
    T u = a / (1.0 + a);
    T num = std::pow(1.0 - u, 2) * (P - P0).dot(P1 - P) + std::pow(u, 2) * (P - P2).dot(P1 - P);
    T den = 2.0 * u * (1.0 - u) * (P1 - P).squaredNorm();
    w1 = num / den;
    return ENUM_NURBS::NURBS_SUCCESS;
}


template<typename T, int dim>
ENUM_NURBS split_arc(const Eigen::Vector<T, dim> &P0, const Eigen::Vector<T, dim> &P1, T w1, const Eigen::Vector<T, dim> &P2,
    Eigen::Vector<T, dim> &Q1, Eigen::Vector<T, dim> &S, Eigen::Vector<T, dim> &R1, T &wqr)
{
    if (w1 == 0.0)
    {
        Q1 = P0 + P1;
        R1 = P2 + P1;
        S = (Q1 + R1) / 2.0;
        wqr = std::sqrt(2.0) / 2.0;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    //else
    Q1 = (P0 + w1 * P1) / (1 + w1);
    R1 = (w1 * P1 + P2) / (1 + w1);
    S = (Q1 + R1) / 2.0;
    wqr = std::sqrt((1.0 + w1) / 2.0);
    return ENUM_NURBS::NURBS_SUCCESS;
}

template<typename T, int dim>
ENUM_NURBS create_open_conic(const Eigen::Vector<T, dim> &P0, const Eigen::Vector<T, dim> &T0, const Eigen::Vector<T, dim> &P2, const Eigen::Vector<T, dim> &T2,
    const Eigen::Vector<T, dim> &P, nurbs_curve<T, dim, true, -1, -1> &nurbs)
{
    Eigen::Vector<T, dim> P1;
    T w1;
    create_one_arc(P0, T0, P2, T2, P, P1, w1);
    if (w1 <= -1.0)
        return ENUM_NURBS::NURBS_ERROR;
    
    int nsegs = -1;
    if (w1 >= 1.0)
        nsegs = 1;
    else
    {
        T angle = angle_between_tow_vector<T, dim>(P0 - P1, P2 - P1);
        if (w1 > 0.0 && angle > M_PI / 3.0)
            nsegs = 1;
        else if (w1 < 0.0 && angle > M_PI / 2.0)
            nsegs = 4;
        else
            nsegs = 2;
    }
    int n = 2 * nsegs;
    Eigen::VectorX<T> knots_vector(4 + n);
    knots_vector.template block<3, 1>(0, 0).setConstant(0.0);
    knots_vector.template block<3, 1>(n + 1, 0).setConstant(1.0);

    Eigen::Matrix<T, dim + 1, Eigen::Dynamic> control_points(dim + 1, n + 1);
    control_points.col(0).template block<dim, 1>(0, 0) = P0;
    control_points(dim, 0) = 1.0;

    control_points.col(n).template block<dim, 1>(0, 0) = P2;
    control_points(dim, n) = 1.0;
    if (nsegs == 1)
    {
        control_points.col(1).template block<dim, 1>(0, 0) = w1 * P1;
        control_points(dim, 1) = w1;
        nurbs.set_control_points(control_points);
        nurbs.set_knots_vector(knots_vector);
        nurbs.set_degree(2);
        std::cout << control_points << std::endl;
        std::cout << knots_vector << std::endl;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    Eigen::Vector<T, dim> Q1, S, R1;
    T wqr;
    split_arc(P0, P1, w1, P2, Q1, S, R1, wqr);
    if (nsegs == 2)
    {
        control_points.col(2).template block<dim, 1>(0, 0) = S;
        control_points(dim, 2) = 1.0;
        control_points.col(1).template block<dim, 1>(0, 0) = wqr * Q1;
        control_points(dim, 1) = wqr;
        control_points.col(3).template block<dim, 1>(0, 0) = wqr * R1;
        control_points(dim, 3) = wqr;
        knots_vector[3] = knots_vector[4] = 0.5;
        nurbs.set_control_points(control_points);
        nurbs.set_knots_vector(knots_vector);
        nurbs.set_degree(2);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    //else
    control_points.col(4).template block<dim, 1>(0, 0) = S;
    control_points(dim, 4) = 1.0;

    w1 = wqr;
    Eigen::Vector<T, dim> HQ1, HS, HR1;
    split_arc(P0, Q1, w1, S, HQ1, HS, HR1, wqr);
    control_points.col(2).template block<dim, 1>(0, 0) = HS;
    control_points(dim, 2) = 1.0;

    control_points.col(1).template block<dim, 1>(0, 0) = wqr * HQ1;
    control_points(dim, 1) = wqr;
    control_points.col(3).template block<dim, 1>(0, 0) = wqr * HR1;
    control_points(dim, 3) = wqr;

    split_arc(S, R1, w1, P2, HQ1, HS, HR1, wqr);
    control_points.col(6).template block<dim, 1>(0, 0) = HS;
    control_points(dim, 6) = 1.0;

    control_points.col(5).template block<dim, 1>(0, 0) = wqr * HQ1;
    control_points(dim, 5) = wqr;
    control_points.col(7).template block<dim, 1>(0, 0) = wqr * HR1;
    control_points(dim, 7) = wqr;
    for (int i = 0; i < 2; ++i)
    {
        knots_vector[i + 3] = 0.25;
        knots_vector[i + 5] = 0.5;
        knots_vector[i + 7] = 0.75;
    }

    return ENUM_NURBS::NURBS_SUCCESS;
}


template<typename T>
struct conic_geometry_feature
{
    int type; // 0 : parabola; 1 : hyperbola; 2 : ellipse; -1 : 此变量无效
    Eigen::Vector<T, 3> focus_or_center; // type = 1时, 此变量表示焦点, 不然表示椭圆或者双曲线的中心
    Eigen::Vector<T, 3> axis_or_major; // type = 1时, 此变量表示准线的方向, 否则表示椭圆的长轴或者双曲线的横轴
    Eigen::Vector<T, 3> vertex; //type = 1时, 此变量有意义
    T major_radii; // type != 1时, 此变量有意义
    T minor_radii; // type != 1时, 此变量有意义
    conic_geometry_feature() = default;
};

//无理nurbs需要转换成有理, 或者单独写一个函数; 此函数没有测试
template<typename T>
ENUM_NURBS evaluate_geometry_feature_from_quadratic_rational_nurbs(const nurbs_curve<T, 3, true, -1, -1> &nurbs, std::vector<conic_geometry_feature<T>> &geometry_feature)
{
    int degree = nurbs.get_degree();
    if (2 != degree)
        return ENUM_NURBS::NURBS_ERROR;
    std::vector<bezier_curve<T, 3, true, -1> *> bezier_curves;
    nurbs.decompose_to_bezier(bezier_curves);
    int curves_count = bezier_curves.size();
    for (int index = 0; index < curves_count; ++index)
    {
        Eigen::Matrix<T, 3 + 1, Eigen::Dynamic> control_points;
        Eigen::Vector<T, 3 + 1> P0, P1, P2;
        bezier_curves[index]->get_control_point(0, P0);
        bezier_curves[index]->get_control_point(1, P1);
        bezier_curves[index]->get_control_point(2, P2);
        T w0, w1, w2;
        bezier_curves[index]->get_weight(0, w0);
        bezier_curves[index]->get_weight(1, w1);
        bezier_curves[index]->get_weight(2, w2);
        Eigen::Vector<T, 3> S = P0 - P1;
        Eigen::Vector<T, 3> T1 = P2 - P1;
        T alpha = S.squaredNorm();
        T gamma = T1.squaredNorm();
        T beta = S.dot(T1);
        T delta = alpha * gamma - beta * beta;  //(S.cross(T1)).squaredNorm();

        if (delta < DEFAULT_ERROR * DEFAULT_ERROR * alpha * gamma || w0 <= 0.0 || w1 <= 0.0 || w2 < DEFAULT_ERROR)
        {
            //退化二次截面曲线或者有无穷远控制点或者负权值
            conic_geometry_feature<T> cgf;
            cgf.type = -1;
            geometry_feature.push_back(std::move(cgf));
            continue;
        }
        //else

        T zat = alpha + gamma + 2.0 * beta; // (S + T).squaredNorm();
        T k = (w0 * w2) / (w1 * w1);
        if (std::abs(k - 1.0) < DEFAULT_ERROR)
        {
            //parabola
            conic_geometry_feature<T> cgf;
            cgf.type = 0;
            cgf.axis_or_major = (S + T1) / std::sqrt(zat);
            cgf.focus_or_center = P1 + (gamma * S + alpha * T1) / zat;
            cgf.vertex = P1 + std::pow((gamma + alpha) / zat, 2) * S + std::pow((alpha + beta) / zat, 2) * T1;
            geometry_feature.push_back(std::move(cgf));
            continue;
        }
        T eit = alpha + gamma - 2.0 * beta; // (S - T).squaredNorm();
        T epsilon = k / (2.0 * (k - 1.0));
        conic_geometry_feature<T> cgf;
        if (k < 1.0)
            cgf.type = 2;
        else
            cgf.type = 1;
        cgf.focus_or_center = P1 + epsilon * (S + T1);

        T b = k * eit + 4.0 * beta;
        T DELTA = std::sqrt(std::pow(b, 2) - 16.0 * delta(k - 1));
        T lamda2 = (DELTA + b) / (4 * delta);
        T lamda1 = (-DELTA + b) / (4 * delta);
        if (std::abs(lamda1 - lamda2) < DEFAULT_ERROR)
        {
            cgf.major_radii = std::sqrt(std::abs(epsilon / lamda1));
            cgf.minor_radii = cgf.major_radii;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        cgf.major_radii = std::sqrt(std::abs(epsilon / lamda1));
        cgf.minor_radii = std::sqrt(std::abs(epsilon / lamda2));
        T x_bar = k / 2.0 - gamma * lamda1;
        T y_bar = k / 2.0 - alpha * lamda1;
        if (std::abs(x_bar) > std::abs(y_bar))
            y_bar = beta * lamda1 - k / 2.0 + 1;
        else
            x_bar = beta * lamda1 - k / 2.0 + 1;
        T rou = alpha * std::pow(x_bar, 2) + 2 * beta * x_bar * y_bar + gamma * std::pow(gamma, 2);
        T x0 = x_bar / rou;
        T y0 = y_bar / rou;
        Eigen::Vector<T, 3> Q1 = P1 + (epsilon + cgf.major_radii * x0) * S + (epsilon + cgf.minor_radii * y0) * T1;
        Eigen::Vector<T, 3> Q2 = P1 + (epsilon - cgf.major_radii * x0) * S + (epsilon - cgf.minor_radii * y0) * T1;
        cgf.axis_or_major = Q2 - Q1;
    }

    return ENUM_NURBS::NURBS_SUCCESS;

}


//创建一个xy平面, 半径为r的, 圆心在原点的整圆, 结果为5次有理bezier曲线(无内部节点), 如果想要得到任意平面整圆, 只需要做一个刚体变换即可
template<typename T, int dim>
ENUM_NURBS create_bezier_circle(T r, bezier_curve<T, dim, true, -1> &bezier)
{
    Eigen::Matrix<T, dim + 1, Eigen::Dynamic> control_points(dim + 1, 6);
    control_points.setConstant(0.0);
    control_points(1, 0) = -5.0 * r;
    control_points(dim, 0) = 5.0;
    
    control_points(0, 1) = 4.0 * r;
    control_points(1, 1) = -r;
    control_points(dim, 1) = 1.0;
    
    control_points(0, 2) = 2.0 * r;
    control_points(1, 2) = r * 3.0;
    control_points(dim, 2) = 1.0;

    control_points(0, 3) = -2.0 * r;
    control_points(1, 3) = r * 3.0;
    control_points(dim, 3) = 1.0;

    control_points(0, 4) = -4.0 * r;
    control_points(1, 4) = r * -1.0;
    control_points(dim, 4) = 1.0;

    control_points(1, 5) = r * -5.0;
    control_points(dim, 5) = 5.0;
    bezier.set_control_points(control_points);
    return ENUM_NURBS::NURBS_SUCCESS;
}

//创建一个xy平面, 半径为r的, 圆心在原点的整圆, 结果为5次有理bezier曲线(无内部节点), 如果想要得到任意平面整圆, 只需要做一个刚体变换即可
template<typename T, int dim>
ENUM_NURBS create_bezier_circle(T r, nurbs_curve<T, dim, true, -1, -1> &nurbs)
{
    Eigen::Matrix<T, dim + 1, Eigen::Dynamic> control_points(dim + 1, 6);
    control_points.setConstant(0.0);
    control_points(1, 0) = -5.0 * r;
    control_points(dim, 0) = 5.0;
    
    control_points(0, 1) = 4.0 * r;
    control_points(1, 1) = -r;
    control_points(dim, 1) = 1.0;
    
    control_points(0, 2) = 2.0 * r;
    control_points(1, 2) = r * 3.0;
    control_points(dim, 2) = 1.0;

    control_points(0, 3) = -2.0 * r;
    control_points(1, 3) = r * 3.0;
    control_points(dim, 3) = 1.0;

    control_points(0, 4) = -4.0 * r;
    control_points(1, 4) = r * -1.0;
    control_points(dim, 4) = 1.0;

    control_points(1, 5) = r * -5.0;
    control_points(dim, 5) = 5.0;
    nurbs.set_control_points(control_points);
    Eigen::VectorX<T> knots_vector(12);
    knots_vector << 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1;
    nurbs.set_knots_vector(knots_vector);
    nurbs.set_degree(5);
    return ENUM_NURBS::NURBS_SUCCESS;
}


