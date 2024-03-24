#include "declare.h"
#include "nurbs_surface.h"
#include "interval_alogrithm.h"
#include <memory>
namespace tnurbs
{
    //check uniqueness and existence
    //先计算曲面的法向(目前只处理非有理的情况，有理只需要稍加修改即可)
    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS eval_normal_and_tangent_interval(const surface_type *u_surface, const surface_type *v_surface, const Box<T, 2> &interval,
                std::array<Interval<T>, (unsigned)dim> &normal_interval, std::array<Interval<T>, (unsigned)dim> &u_tangent_interval, std::array<Interval<T>, (unsigned)dim> &v_tangent_interval)
    {
        surface_type *u_sub_surface = new surface_type();
        surface_type *v_sub_surface = new surface_type();
        
        std::unique_ptr<surface_type> raii = std::unique_ptr<surface_type>(u_sub_surface);
        std::unique_ptr<surface_type> raii2 = std::unique_ptr<surface_type>(v_sub_surface);
        
        Box<T, 2> interval_copy = interval;
        u_surface->sub_divide(interval_copy, *u_sub_surface);
        v_surface->sub_divide(interval_copy, *v_sub_surface);
        
        Box<T, dim> u_box, v_box;
        u_sub_surface->get_box(u_box);
        v_sub_surface->get_box(v_box);

        //将u_box,v_box 转为interval vector
        for (int index = 0; index < dim; ++index)
        {
            u_tangent_interval[index].m_interval[0] = u_box.Min[index];
            u_tangent_interval[index].m_interval[1] = u_box.Max[index];

            v_tangent_interval[index].m_interval[0] = v_box.Min[index];
            v_tangent_interval[index].m_interval[1] = v_box.Max[index];
        }
        normal_interval = vector_product<T>(u_tangent_interval, v_tangent_interval);
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    bool eval_priori_enclosure(const surface_type *u_surface_left, const surface_type *v_surface_left,
                                        const surface_type *u_surface_right, const surface_type *v_surface_right, const Box<T, 4> &box_domain,
                                        const Box<T, 4> &initial_box, Box<T, 4> &priori_enclosure, const T &step_size, T param_space_eps)
    {
        //目前仅仅支持三维
        static_assert(3 == dim, "dim != 3 is not supported");
        std::array<Interval<T>, (unsigned)dim> left_normal_interval, right_normal_interval;
        std::array<Interval<T>, (unsigned)dim> left_u_tangent_interval, left_v_tangent_interval;
        std::array<Interval<T>, (unsigned)dim> right_u_tangent_interval, right_v_tangent_interval;

        Box<T, 2> left_param_box(initial_box.Min.template block<2, 1>(0, 0), initial_box.Max.template block<2, 1>(0, 0));
        Box<T, 2> right_param_box(initial_box.Min.template block<2, 1>(2, 0), initial_box.Max.template block<2, 1>(2, 0));
        eval_normal_and_tangent_interval<surface_type>(u_surface_left, v_surface_left, left_param_box, left_normal_interval, left_u_tangent_interval, left_v_tangent_interval);
        eval_normal_and_tangent_interval<surface_type>(u_surface_right, v_surface_right, right_param_box, right_normal_interval, right_u_tangent_interval, right_v_tangent_interval);

        Interval<T> u_normal_dot = dot_product(left_normal_interval, left_normal_interval);
        Interval<T> v_normal_dot = dot_product(right_normal_interval, right_normal_interval);

        std::array<Interval<T>, (unsigned)dim> c = vector_product(left_normal_interval, right_normal_interval);

        T len = std::sqrt(u_normal_dot.m_interval[1]) * std::sqrt(v_normal_dot.m_interval[1]);
        Interval<T> len_interval(len, len);
        std::array<Interval<T>, (unsigned)dim> temp;
         if (false == divide<T, dim>(c, len_interval, temp))
            return false;


        int count = 0;
        for (int index = 0; index < dim; ++index)
        {
            if (0.0 >= temp[index].m_interval[0] - TDEFAULT_ERROR<T>::value &&  0.0 <= temp[index].m_interval[1] + PRECISION<T>::value)
                ++count;
        }
        
        //奇异点
        if (count == dim)
            return false;

        Interval<T> c_dot = dot_product(c, c);
        c_dot.m_interval[0] = std::sqrt(c_dot.m_interval[0]);
        c_dot.m_interval[1] = std::sqrt(c_dot.m_interval[1]);
        if (divide<T, (unsigned)dim>(c, c_dot, c) == false)
            return false;


        Interval<T> alpha, beta, u, v;
        if (false == divide(dot_product(left_normal_interval, vector_product(c, left_v_tangent_interval)), u_normal_dot, alpha))
            return false;
        if (false == divide(dot_product(left_normal_interval, vector_product(left_u_tangent_interval, c)), u_normal_dot, beta))
            return false;
        if (false == divide(dot_product(right_normal_interval, vector_product(c, right_v_tangent_interval)), v_normal_dot, u))
            return false;
        if (false == divide(dot_product(right_normal_interval, vector_product(right_u_tangent_interval, c)), v_normal_dot, v))
            return false;


        Box<T, 4> box_interval, box_interval2;
        box_interval2.Min[0] = alpha.m_interval[0];
        box_interval2.Min[1] = beta.m_interval[0];
        box_interval2.Min[2] = u.m_interval[0];
        box_interval2.Min[3] = v.m_interval[0];

        box_interval2.Max[0] = alpha.m_interval[1];
        box_interval2.Max[1] = beta.m_interval[1];
        box_interval2.Max[2] = u.m_interval[1];
        box_interval2.Max[3] = v.m_interval[1];
        box_interval2.scale(step_size);
        

        Interval<T> step_box(0.0, step_size);
        alpha = mutilply(alpha, step_box);
        beta = mutilply(beta, step_box);
        u = mutilply(u, step_box);
        v = mutilply(v, step_box);


        box_interval.Min[0] = alpha.m_interval[0];
        box_interval.Min[1] = beta.m_interval[0];
        box_interval.Min[2] = u.m_interval[0];
        box_interval.Min[3] = v.m_interval[0];

        box_interval.Max[0] = alpha.m_interval[1];
        box_interval.Max[1] = beta.m_interval[1];
        box_interval.Max[2] = u.m_interval[1];
        box_interval.Max[3] = v.m_interval[1];

        Eigen::Vector4<T> param_dis = box_interval.Max - box_interval.Min;
        if (param_dis[0] + param_dis[1] > param_space_eps || param_dis[2] + param_dis[3] > param_space_eps)
            return false;

        Box<T, 4> param_box;
        initial_box.plus(box_interval, PRECISION<T>::value).intersect(box_domain, param_box);
        if (priori_enclosure.is_contain_box(param_box))
        {
            if (initial_box.plus(box_interval2, PRECISION<T>::value).intersect(box_domain, priori_enclosure))
                return true;
        }
             
        //else
        priori_enclosure = param_box;
        return false;
    }


    //计算第一基本型和第二基本型
    template<typename T, int dim>
    ENUM_NURBS eval_first_and_second_fundenmental(const Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> ders, const Eigen::Vector<T, dim> &normal, 
                                                    Eigen::Matrix2<T> &first, Eigen::Matrix2<T> &second)
    {
        static_assert(3 == dim, "3 != dim");
        first(0, 0) = ders(1, 0).dot(ders(1, 0));
        first(1, 0) = first(0, 1) = ders(1, 0).dot(ders(0, 1));
        first(1, 1) = ders(0, 1).dot(ders(0, 1));

        second(0, 0) = ders(2, 0).dot(normal);
        second(1, 0) = second(0, 1) = ders(1, 1).dot(normal);
        second(1, 1) = ders(0, 2).dot(normal);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    
    /// @brief 计算两个曲面交线在各自参数域的原像曲线的1-2阶微分
    /// @tparam surface_type 曲面类型
    /// @tparam T double, flaot...
    /// @tparam dim 维数(目前仅仅支持三维)
    /// @param left_surface 第一个曲面
    /// @param right_surface 第二个曲面
    /// @param u 第一个曲面的参数u
    /// @param v 第一个曲面的参数v
    /// @param s 第一个曲面的参数s
    /// @param t 第一个曲面的参数t
    /// @param param_ders param_ders.col(0)为原像曲线的一阶微分, param_ders.col(1)为原像曲线的二阶微分
    /// @param space_ders space_ders.col(0)为交线的点, space_ders.col(1)为交线的一阶微分, space_ders.col(1)为交先的二阶微分
    /// @return 错误码
    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS eval_preiamge_ders(const surface_type *left_surface, const surface_type *right_surface, T u, T v, T s, T t, 
                                        Eigen::Matrix<T, 4, 2> &param_ders, Eigen::Matrix<T, dim, 3> &space_ders)
    {
        static_assert(3 == dim, "3 != dim");

        Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> ders;
        space_ders.col(0) = ders(0, 0);
        left_surface->template derivative_on_surface<2, 2>(u, v, ders);
        Eigen::Vector<T, dim> normal = ders(1, 0).cross(ders(0, 1));
        normal.normalize();
        Eigen::Matrix2<T> left_first_fundenmental;
        Eigen::Matrix2<T> left_second_fundenmental;
        eval_first_and_second_fundenmental(ders, normal, left_first_fundenmental, left_second_fundenmental);

        Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> ders2;
        right_surface->template derivative_on_surface<2, 2>(s, t, ders2);
        Eigen::Vector<T, dim> normal2 = ders2(1, 0).cross(ders2(0, 1));
        normal2.normalize();
        Eigen::Matrix2<T> right_first_fundenmental;
        Eigen::Matrix2<T> right_second_fundenmental;
        eval_first_and_second_fundenmental(ders2, normal2, right_first_fundenmental, right_second_fundenmental);

        Eigen::Vector<T, dim> tangent_vec = normal.cross(normal2);

        // if (tangent_vec.squaredNorm() < TDEFAULT_ERROR<T>::value)
        // {
        //     //两曲面在交点处的法向共线; 暂时不处理; TODO
        //     return ENUM_NURBS::NURBS_ERROR;
        // }


        tangent_vec.normalize();
        space_ders.col(1) = tangent_vec;
        Eigen::Vector<T, 2> dot_ders;
        dot_ders[0] = tangent_vec.dot(ders(1, 0));
        dot_ders[1] = tangent_vec.dot(ders(0, 1));
        Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(left_first_fundenmental);
        param_ders.template block<2, 1>(0, 0) = matSvd.solve(dot_ders);
        if (matSvd.info() !=  Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        
        dot_ders[0] = tangent_vec.dot(ders2(1, 0));
        dot_ders[1] = tangent_vec.dot(ders2(0, 1));
        Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd2(right_first_fundenmental);
        param_ders.template block<2, 1>(2, 0) = matSvd2.solve(dot_ders);
        if (matSvd2.info() !=  Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;

        //计算两曲面的交线的二阶微分(弧长参数下)
        
        //1.计算 d(u, v) 的第二基本型的值
        Eigen::Vector<T, 2> second_fundenmental_value;
        second_fundenmental_value[0] = param_ders.template block<2, 1>(0, 0).col(0).transpose() * left_second_fundenmental * param_ders.template block<2, 1>(0, 0).col(0);
        second_fundenmental_value[1] = param_ders.template block<2, 1>(2, 0).transpose() * right_second_fundenmental * param_ders.template block<2, 1>(2, 0);
        Eigen::Matrix2<T> mat;
        mat(0, 0) = mat(1, 1) = 0.0;
        mat(0, 1) = mat(1, 0) = normal.dot(normal2);
        Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd3(mat);
        Eigen::Vector2<T> coeff = matSvd3.solve(second_fundenmental_value);
        if (matSvd3.info() !=  Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        
        //交线的二阶导数
        Eigen::Vector<T, dim> alpha_dd = coeff[0] * normal + coeff[1] * normal2;
        space_ders.col(2) = alpha_dd;
        Eigen::Vector<T, dim> L = std::pow(param_ders(0, 0), 2) * ders(2, 0) + 2 * param_ders(0, 0) * param_ders(1, 0) * ders(1, 1) + std::pow(param_ders(1, 0), 2) * ders(0, 2);

        Eigen::Vector<T, dim> temp = alpha_dd - L;
        dot_ders[0] = temp.dot(ders(1, 0));
        dot_ders[1] = temp.dot(ders(0, 1));
        param_ders.template block<2, 1>(0, 1) = matSvd.solve(dot_ders);
        if (matSvd.info() !=  Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        
        L  = std::pow(param_ders(2, 0), 2) * ders2(2, 0) + 2 * param_ders(2, 0) * param_ders(3, 0) * ders2(1, 1) + std::pow(param_ders(3, 0), 2) * ders2(0, 2);
        temp = alpha_dd - L;
        dot_ders[0] = temp.dot(ders2(1, 0));
        dot_ders[1] = temp.dot(ders2(0, 1));
        param_ders.template block<2, 1>(2, 1) = matSvd2.solve(dot_ders);
        if (matSvd2.info() !=  Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    

    template<typename T, int dim>
    ENUM_NURBS guess_step_size(const Eigen::Matrix<T, dim, 3> &ders, T eps, T &step_size)
    {
        T y_norm = ders.col(0).norm();
        T y_d_norm = ders.col(1).norm();
        T y_dd_norm = ders.col(2).norm();
        step_size = std::sqrt(eps * y_norm / (2.0 * y_dd_norm));
        if (y_norm < step_size * y_d_norm || step_size < eps)
        {
            step_size = eps * y_d_norm / (2.0 * y_dd_norm);
        } 
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    //TODO: 参数域边界需要处理
    template<typename surface_type, bool direction, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS eval_step_size_and_yBox(const surface_type *left_surf, const surface_type *right_surf, 
                                        const surface_type *left_u_tangent_surface, const surface_type *left_v_tangent_surface, 
                                        const surface_type *right_u_tangent_surface, const surface_type *right_v_tangent_surface, const Box<T, 4> & box_domain,
                                        const Box<T, 4> &initial_box, Box<T, 4> &priori_enclosure, Eigen::Vector<T, 4> &next_initinal_param, T param_space_eps, T &step_size)
    {
        Eigen::Vector<T, 4> param = (initial_box.Min + initial_box.Max) / 2.0;

        Eigen::Matrix<T, 4, 2> params_ders;
        Eigen::Matrix<T, dim, 3> space_ders;
        eval_preiamge_ders(left_surf, right_surf, param[0], param[1], param[2], param[3], params_ders, space_ders);

        //此处需要修改
        T temp = direction == true ? (2.0 * param_space_eps) / params_ders.col(0).norm() : (-2.0 * param_space_eps) / params_ders.col(0).norm();
        if (std::abs(temp) < std::abs(step_size))
            step_size = temp;

        //TODO : 10为最大循环次数, 以后整理一下
        priori_enclosure = initial_box;
        int index = 0;
        for (;index < 10; ++index)
        {
            //需要修改一点点
            if (true == eval_priori_enclosure(left_u_tangent_surface, left_v_tangent_surface, right_u_tangent_surface, right_v_tangent_surface, box_domain, initial_box, priori_enclosure, step_size, param_space_eps))
            {
                break;
            }
                
            step_size /= 2.0;
            if (std::abs(step_size) < TDEFAULT_ERROR<T>::value)
                return ENUM_NURBS::NURBS_ERROR;
        }
        if (index == 10)
            return ENUM_NURBS::NURBS_ERROR;
        
        next_initinal_param = params_ders.col(0) * step_size + param;
        for (int j = 0; j < 4; ++j)
        {
            if (next_initinal_param[j] < 0.0)
                next_initinal_param[j] = 0.0;
            else if (next_initinal_param[j] > 1.0)
                next_initinal_param[j] = 1.0;
        }


        return ENUM_NURBS::NURBS_SUCCESS;
    }



    //intersect points
    //surf_surf_int 负责释放内存
    template<typename T, int dim>
    struct surf_surf_intersect_point
    {
        surf_surf_intersect_point<T, dim> *m_next;
        Eigen::Vector<T, dim> m_point;
        Eigen::Vector2<T> m_left_uv;
        Eigen::Vector2<T> m_right_uv;
        
        //是否相切
        bool m_is_tangent;

        //是否时孤立交点
        bool m_is_isolate;

        //奇异点必然为起始点或者终止点
        bool m_is_singular;
        std::vector<surf_surf_intersect_point*> m_connect_points;

        //TODO:几阶接触
    };

    //面面相交的结果
    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    struct surf_surf_int
    {
        surface_type* m_left_surf;
        surface_type* m_right_surf;
        std::vector<surf_surf_intersect_point<T, dim>*> m_intersect_points;

        //每个交点链表的交点的个数
        std::vector<int> m_intersect_points_count;
        ~surf_surf_int()
        {
            delete m_left_surf;
            delete m_right_surf;

            for (auto it = m_intersect_points.begin(); it != m_intersect_points.end(); ++it)
            {
                surf_surf_intersect_point<T, dim> *current = *it;
                while (current)
                {
                    surf_surf_intersect_point<T, dim> *next = current->m_next;
                    delete current;
                    current = next;
                }
                *it = nullptr;
            }
        }
    };


    //迭代加细交点(周期性曲面待处理)
    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS intersect_point_iteration(const surface_type *left_surf, const surface_type *right_surf, 
                                        const surface_type *left_u_tangent_surface, const surface_type *left_v_tangent_surface, 
                                        const surface_type *right_u_tangent_surface, const surface_type *right_v_tangent_surface,
                                        const Box<T, 4> &domian, Eigen::Vector<T, 4> current_param, Eigen::Vector<T, 4> &intersect_param, T eps)
    {
        // bool left_surf_u_closed = left_surf->is_u_closed();
        // bool left_surf_v_closed = left_surf->is_v_closed();
        // bool right_surf_u_closed = right_surf->is_u_closed();
        // bool right_surf_v_closed = right_surf->is_v_closed();

        Eigen::Vector<T, dim> left_point, right_point;
        left_surf->point_on_surface(current_param[0], current_param[1], left_point);
        right_surf->point_on_surface(current_param[2], current_param[3], right_point);
        Eigen::Vector<T, dim> vec = right_point - left_point;
        T min_distance = vec.squaredNorm();
        intersect_param = current_param;

        //迭代次数需要修改
        for (int loop_index = 0; loop_index < 10 * SURFACE_ITERATE_DEEP; ++loop_index)
        {
            Eigen::Matrix<T, dim, 4> mat;
            Eigen::Vector<T, dim> temp;
            left_u_tangent_surface->point_on_surface(current_param[0], current_param[1], temp);
            mat.col(0) = temp;
            left_v_tangent_surface->point_on_surface(current_param[0], current_param[1], temp);
            mat.col(1) = temp;
            right_u_tangent_surface->point_on_surface(current_param[2], current_param[3], temp);
            mat.col(2) = -1.0 * temp;
            right_v_tangent_surface->point_on_surface(current_param[2], current_param[3], temp);
            mat.col(3) = -1.0 * temp;

            Eigen::JacobiSVD<Eigen::Matrix<T, dim, 4>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
            Eigen::Vector<T, 4> delta = matSvd.solve(vec);
            if (matSvd.info() !=  Eigen::Success)
                return ENUM_NURBS::NURBS_ERROR;
            
            Eigen::Vector<T, 4> next_param = current_param + delta;
            for (int index = 0; index < 4; ++index)
            {
                if (next_param[index] < domian.Min[index])
                    next_param[index] = domian.Min[index];
                else if (next_param[index] > domian.Max[index])
                    next_param[index] = domian.Max[index];
            }


            // bool is_closed_flag = cur->is_closed();
            // if (left_surf_u_closed)
            // {
            //     if (next_u < min)
            //         next_u = max - (min - next_u);
            //     else if (next_u > max)
            //         next_u = min + (next_u - max);
            // }
            // else
            // {
            //     if (next_u < min)
            //         next_u = min;
            //     else if (next_u > max)
            //         next_u = max;
            // }

            current_param = next_param;


            left_surf->point_on_surface(current_param[0], current_param[1], left_point);
            right_surf->point_on_surface(current_param[2], current_param[3], right_point);
            vec = right_point - left_point;
            T distance = vec.squaredNorm();
            if (distance < min_distance)
            {
                min_distance = distance;
                intersect_param = current_param;

                bool flag = true;
                for (int index = 0; index < dim; ++index)
                {
                    if (vec[index] > PRECISION<T>::value)
                    {
                        flag = false;
                        break;
                    }
                }
                if (flag == true)
                    return ENUM_NURBS::NURBS_SUCCESS;
            }
        }

        return ENUM_NURBS::NURBS_ERROR;
    }


    //TODO: surface_type参数前面加const
    template<typename surface_type, bool direction = true, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS ODESolver(Box<T, 4> initial_box,  surface_type *left, surface_type *right, surf_surf_int<surface_type> *intersection, T eps = TDEFAULT_ERROR<T>::value)
    {
        surface_type *left_u_tangent_surface = new surface_type();
        surface_type *left_v_tangent_surface = new surface_type();
        surface_type *right_u_tangent_surface = new surface_type();
        surface_type *right_v_tangent_surface = new surface_type();
        std::unique_ptr<surface_type> raii = std::unique_ptr<surface_type>(left_u_tangent_surface);
        std::unique_ptr<surface_type> raii2 = std::unique_ptr<surface_type>(left_v_tangent_surface);
        std::unique_ptr<surface_type> raii3 = std::unique_ptr<surface_type>(right_u_tangent_surface);
        std::unique_ptr<surface_type> raii4 = std::unique_ptr<surface_type>(right_v_tangent_surface);

        left->tangent_u_surface(*left_u_tangent_surface);
        left->tangent_v_surface(*left_v_tangent_surface);
        right->tangent_u_surface(*right_u_tangent_surface);
        right->tangent_v_surface(*right_v_tangent_surface);

        Box<T, dim> left_u_tangent_surface_box, left_v_tangent_surface_box, right_u_tangent_surface_box, right_v_tangent_surface_box;
        left_u_tangent_surface->get_box(left_u_tangent_surface_box);
        left_v_tangent_surface->get_box(left_v_tangent_surface_box);
        right_u_tangent_surface->get_box(right_u_tangent_surface_box);
        right_v_tangent_surface->get_box(right_v_tangent_surface_box);

        Eigen::Vector<T, dim> origin;
        origin.setConstant(0.0);
        T l1 = left_u_tangent_surface_box.eval_maximal_distance(origin);
        T l2 = left_v_tangent_surface_box.eval_maximal_distance(origin);
        T l3 = right_u_tangent_surface_box.eval_maximal_distance(origin);
        T l4 = right_v_tangent_surface_box.eval_maximal_distance(origin);
        T l = std::max({l1, l2, l3, l4});
        T param_space_eps = eps / l;
        if (param_space_eps <= PRECISION<T>::value)
            return ENUM_NURBS::NURBS_ERROR;

        Box<T, 4> priori_enclosure;
        Eigen::Vector<T, 4> box_size;
        T box_size_num = std::min(TDEFAULT_ERROR<T>::value, param_space_eps / 4.0);
        box_size.setConstant(box_size_num);

        surf_surf_intersect_point<T, dim> *current_intersection = intersection->m_intersect_points.back();
        Box<T, 4> product_box = left->get_interval().product_box(right->get_interval());

        //加细测试用例, 后面需要删除
        Eigen::Vector<T, 4> intersect_param2;
        Eigen::Vector<T, 4> x(0.02, 0.66630158097819048, 0.71793367976331557, 0.97999999999999987);
        intersect_point_iteration<surface_type>(left, right, left_u_tangent_surface, left_v_tangent_surface, right_u_tangent_surface, right_v_tangent_surface, 
                                        initial_box, x, intersect_param2, eps);


        initial_box.Min = intersect_param2 - box_size;
        initial_box.Max = intersect_param2 + box_size;

        Box<T, 4> x_box2;
        initial_box.intersect(product_box, x_box2);

        int index = 0;

        //需要修改
        T step_size = direction == true ? MAX_ITERATE_DEEP : -MAX_ITERATE_DEEP;
        // T step_size = direction == true ? eps : -eps;
        for (; index < MAXINTERSETORPOINTNUMBER; ++index)
        {
            Eigen::Vector<T, 4> next_init_param;
            ENUM_NURBS flag = eval_step_size_and_yBox<surface_type, direction>(left, right, left_u_tangent_surface, left_v_tangent_surface, right_u_tangent_surface, right_v_tangent_surface, product_box, 
                                                initial_box, priori_enclosure, next_init_param, param_space_eps, step_size);

            if (flag != ENUM_NURBS::NURBS_SUCCESS)
            {
                intersection->m_intersect_points_count.push_back(index);
                return ENUM_NURBS::NURBS_SUCCESS;
            }

            Eigen::Vector<T, 4> intersect_param;
            flag = intersect_point_iteration<surface_type>(left, right, left_u_tangent_surface, left_v_tangent_surface, right_u_tangent_surface, right_v_tangent_surface, 
                                        priori_enclosure, next_init_param, intersect_param, eps);
            if (flag != ENUM_NURBS::NURBS_SUCCESS)
            {
                step_size /= 2.0;
                index--;
                continue;
            }
            step_size = direction == true ? MAX_ITERATE_DEEP : -MAX_ITERATE_DEEP;
            // step_size = direction == true ? eps : -eps;
            surf_surf_intersect_point<T, dim> *new_intersect = new surf_surf_intersect_point<T, dim>();
            left->point_on_surface(intersect_param[0], intersect_param[1], new_intersect->m_point);
            new_intersect->m_left_uv = intersect_param.block(0, 0, 2, 1);
            new_intersect->m_right_uv = intersect_param.block(2, 0, 2, 1);

            current_intersection->m_next = new_intersect;
            current_intersection = new_intersect;

            initial_box.Min = intersect_param - box_size;
            initial_box.Max = intersect_param + box_size;

            Box<T, 4> x_box;
            initial_box.intersect(product_box, x_box);
            initial_box = x_box;

        }
        if (index == MAXINTERSETORPOINTNUMBER)
        {
            //TODO:析构内存
            return ENUM_NURBS::NURBS_ERROR;
        }
        intersection->m_intersect_points_count.push_back(index);
        return ENUM_NURBS::NURBS_SUCCESS;
    }



}