#include "declare.h"
#include "nurbs_surface.h"
#include "interval_alogrithm.h"
#include <memory>
#include "bezier_curve_int.h"
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
    

    //intersect points
    //surf_surf_int 负责释放内存
    template<typename T, int dim>
    struct surf_surf_intersect_point
    {
        surf_surf_intersect_point<T, dim> *m_next;
        Eigen::Vector<T, dim> m_point;
        Eigen::Vector4<T> m_uv;
        
        //是否相切.TODO:几阶接触
        bool m_is_tangent;

        //是否时孤立交点
        bool m_is_isolate;

        //奇异点必然为起始点或者终止点
        //奇异点的意思是在此点的切平面上存在至少两个方向v1, v2, 使得两个相交的平面在这两个方向上任意阶偏导数相等(弧长参数下), 奇异点必然相切
        bool m_is_singular;
        std::vector<surf_surf_intersect_point*> m_connect_points;

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

    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    class trace_nurbs_surface
    {
    private:
        surface_type* m_left_surface;
        surface_type* m_right_surface;
        
        //ders of left surface
        surface_type* m_left_u_tangent_surface;
        surface_type* m_left_v_tangent_surface;
        
        //ders of right surface
        surface_type* m_right_u_tangent_surface;
        surface_type* m_right_v_tangent_surface;

        T m_param_space_eps;

        Eigen::Vector<T, 4> m_box_size;

        Box<T, 4> m_product_box;

        Interval<T> m_alpha;
        Interval<T> m_beta;
        Interval<T> m_u;
        Interval<T> m_v;

    public:
        trace_nurbs_surface() = delete;
        trace_nurbs_surface(surface_type* left_surface, surface_type* right_surface) : m_left_surface(left_surface), m_right_surface(right_surface) { }
        ~trace_nurbs_surface() 
        { 
            //delete m_left_surface;
            //delete m_right_surface;
            delete m_left_u_tangent_surface;
            delete m_left_v_tangent_surface;
            delete m_right_u_tangent_surface;
            delete m_right_v_tangent_surface;
        }
        ENUM_NURBS init(T space_eps)
        {
            m_left_u_tangent_surface = new surface_type();
            m_left_v_tangent_surface = new surface_type();
            m_right_u_tangent_surface = new surface_type();
            m_right_v_tangent_surface = new surface_type();

            m_left_surface->tangent_u_surface(*m_left_u_tangent_surface);
            m_left_surface->tangent_v_surface(*m_left_v_tangent_surface);
            m_right_surface->tangent_u_surface(*m_right_u_tangent_surface);
            m_right_surface->tangent_v_surface(*m_right_v_tangent_surface);

            Box<T, dim> left_u_tangent_surface_box, left_v_tangent_surface_box, right_u_tangent_surface_box, right_v_tangent_surface_box;
            m_left_u_tangent_surface->get_box(left_u_tangent_surface_box);
            m_left_v_tangent_surface->get_box(left_v_tangent_surface_box);
            m_right_u_tangent_surface->get_box(right_u_tangent_surface_box);
            m_right_v_tangent_surface->get_box(right_v_tangent_surface_box);

            Eigen::Vector<T, dim> origin;
            origin.setConstant(0.0);
            T l1 = left_u_tangent_surface_box.eval_maximal_distance(origin);
            T l2 = left_v_tangent_surface_box.eval_maximal_distance(origin);
            T l3 = right_u_tangent_surface_box.eval_maximal_distance(origin);
            T l4 = right_v_tangent_surface_box.eval_maximal_distance(origin);
            T l = std::max({ l1, l2, l3, l4 });
            m_param_space_eps = space_eps / l;
            if (m_param_space_eps <= PRECISION<T>::value)
                return ENUM_NURBS::NURBS_ERROR;


            m_product_box = m_left_surface->get_interval().product_box(m_right_surface->get_interval());

            T box_size_num = std::min(SPACE_PARAM_TIMES<T>::value * PRECISION<T>::value, m_param_space_eps / 16.0);
            m_box_size.setConstant(box_size_num);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        bool eval_priori_enclosure_inner(const Box<T, 4>& initial_box, T step_size, Box<T, 4>& priori_enclosure, T param_space_eps)
        {
            //目前仅仅支持三维
            static_assert(3 == dim, "dim != 3 is not supported");

            //计算交线的切向, 两个曲面非相切的时候; TODO:相切的时候
            std::array<Interval<T>, (unsigned)dim> left_normal_interval, right_normal_interval;
            std::array<Interval<T>, (unsigned)dim> left_u_tangent_interval, left_v_tangent_interval;
            std::array<Interval<T>, (unsigned)dim> right_u_tangent_interval, right_v_tangent_interval;

            Box<T, 2> left_param_box(priori_enclosure.Min.template block<2, 1>(0, 0), priori_enclosure.Max.template block<2, 1>(0, 0));
            Box<T, 2> right_param_box(priori_enclosure.Min.template block<2, 1>(2, 0), priori_enclosure.Max.template block<2, 1>(2, 0));
            eval_normal_and_tangent_interval<surface_type>(m_left_u_tangent_surface, m_left_v_tangent_surface, left_param_box, left_normal_interval, left_u_tangent_interval, left_v_tangent_interval);
            eval_normal_and_tangent_interval<surface_type>(m_right_u_tangent_surface, m_right_v_tangent_surface, right_param_box, right_normal_interval, right_u_tangent_interval, right_v_tangent_interval);

            Interval<T> u_normal_dot = dot_product(left_normal_interval, left_normal_interval);
            Interval<T> v_normal_dot = dot_product(right_normal_interval, right_normal_interval);
            std::array<Interval<T>, (unsigned)dim> c = vector_product(left_normal_interval, right_normal_interval);

            Interval<T> c_dot = dot_product(c, c);

            //正常不应该在此return false
            if (c_dot.m_interval[0] <= TDEFAULT_ERROR<T>::value * TDEFAULT_ERROR<T>::value)
                return false;

            c_dot.m_interval[0] = std::sqrt(c_dot.m_interval[0]);
            c_dot.m_interval[1] = std::sqrt(c_dot.m_interval[1]);
            if (divide<T, (unsigned)dim>(c, c_dot, c) == false)
                return false;

            if (false == divide(dot_product(left_normal_interval, vector_product(c, left_v_tangent_interval)), u_normal_dot, m_alpha))
                return false;
            if (false == divide(dot_product(left_normal_interval, vector_product(left_u_tangent_interval, c)), u_normal_dot, m_beta))
                return false;
            if (false == divide(dot_product(right_normal_interval, vector_product(c, right_v_tangent_interval)), v_normal_dot, m_u))
                return false;
            if (false == divide(dot_product(right_normal_interval, vector_product(right_u_tangent_interval, c)), v_normal_dot, m_v))
                return false;

            Interval<T> alpha, beta, u, v;
            Interval<T> step_box(0.0, step_size);
            alpha = mutilply(m_alpha, step_box);
            beta = mutilply(m_beta, step_box);
            u = mutilply(m_u, step_box);
            v = mutilply(m_v, step_box);

            Box<T, 4> box_interval;
            box_interval.Min[0] = alpha.m_interval[0];
            box_interval.Min[1] = beta.m_interval[0];
            box_interval.Min[2] = u.m_interval[0];
            box_interval.Min[3] = v.m_interval[0];

            box_interval.Max[0] = alpha.m_interval[1];
            box_interval.Max[1] = beta.m_interval[1];
            box_interval.Max[2] = u.m_interval[1];
            box_interval.Max[3] = v.m_interval[1];

            Eigen::Vector4<T> param_dis = box_interval.Max - box_interval.Min;
            if (param_dis[0] + param_dis[1] > m_param_space_eps || param_dis[2] + param_dis[3] > m_param_space_eps)
                return false;

            Box<T, 4> temp_priori_enclosure;
            if (true == initial_box.plus(box_interval, PRECISION<T>::value).intersect(m_product_box, temp_priori_enclosure))
            {
                priori_enclosure = temp_priori_enclosure;
                return true;
            }
            return false;
        }

        //TODO:删掉
        bool check_priori_enclosure(const Box<T, 4>& next_priori_enclosure, Box<T, 4>& priori_enclosure) 
        {
            if (priori_enclosure.is_contain_box(next_priori_enclosure))
            {
                return true;
            }
            return false;
        }

        bool eval_priori_enclosure(const Box<T, 4>& initial_box, T &step_size, Box<T, 4>& priori_enclosure, T param_space_eps) 
        {
            if (false == eval_priori_enclosure_inner(initial_box, step_size, priori_enclosure, param_space_eps))
            {
                return false;
            }
            Box<T, 4> next_priori_enclosure = priori_enclosure;
            step_size /= 2.0;
            if (false == eval_priori_enclosure_inner(initial_box, step_size, next_priori_enclosure, param_space_eps))
            {
                return false;
            }
            if (check_priori_enclosure(next_priori_enclosure, priori_enclosure))
            {
                Box<T, 4> box_interval;
                box_interval.Min[0] = m_alpha.m_interval[0];
                box_interval.Min[1] = m_beta.m_interval[0];
                box_interval.Min[2] = m_u.m_interval[0];
                box_interval.Min[3] = m_v.m_interval[0];

                box_interval.Max[0] = m_alpha.m_interval[1];
                box_interval.Max[1] = m_beta.m_interval[1];
                box_interval.Max[2] = m_u.m_interval[1];
                box_interval.Max[3] = m_v.m_interval[1];
                box_interval.scale(step_size);

                Box<T, 4> temp_priori_enclosure;
                if (true == initial_box.plus(box_interval, PRECISION<T>::value).intersect(m_product_box, temp_priori_enclosure))
                {
                    priori_enclosure = temp_priori_enclosure;
                    return true;
                }
                return false;
            }
            return false;


        }

        //TODO: 参数域边界需要处理
        template<bool direction>
        ENUM_NURBS eval_step_size_and_yBox(const Box<T, 4>& initial_box, const Eigen::Vector<T, 4> &param, Box<T, 4>& priori_enclosure, Eigen::Vector<T, 4>& next_initinal_param, T param_space_eps)
        {
            Eigen::Matrix<T, 4, 2> params_ders;
            params_ders.setConstant(0.0);
            Eigen::Matrix<T, dim, 3> space_ders;
            space_ders.setConstant(0.0);
            eval_preiamge_ders<surface_type, typename surface_type::Type, surface_type::dimension>(m_left_surface, m_right_surface, param[0], param[1], param[2], param[3], params_ders, space_ders);

            T step_size = direction == true ? param_space_eps / (2.0 * params_ders.col(0).norm()) : param_space_eps / (-2.0 * params_ders.col(0).norm());

            if (false == eval_priori_enclosure(initial_box, step_size, priori_enclosure, param_space_eps))
                return ENUM_NURBS::NURBS_ERROR;

            next_initinal_param = priori_enclosure.get_middle_point();
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        //迭代加细交点(周期性曲面待处理)
        ENUM_NURBS intersect_point_iteration(const Box<T, 4>& domian, Eigen::Vector<T, 4> current_param, Eigen::Vector<T, 4>& intersect_param) 
        {
            // bool left_surf_u_closed = left_surf->is_u_closed();
            // bool left_surf_v_closed = left_surf->is_v_closed();
            // bool right_surf_u_closed = right_surf->is_u_closed();
            // bool right_surf_v_closed = right_surf->is_v_closed();

            Eigen::Vector<T, dim> left_point, right_point;
            m_left_surface->point_on_surface(current_param[0], current_param[1], left_point);
            m_right_surface->point_on_surface(current_param[2], current_param[3], right_point);
            Eigen::Vector<T, dim> vec = right_point - left_point;
            T min_distance = vec.squaredNorm();
            intersect_param = current_param;

            //迭代次数需要修改
            for (int loop_index = 0; loop_index < SURFACE_ITERATE_DEEP; ++loop_index)
            {
                Eigen::Matrix<T, dim, 4> mat;
                Eigen::Vector<T, dim> temp;
                m_left_u_tangent_surface->point_on_surface(current_param[0], current_param[1], temp);
                mat.col(0) = temp;
                m_left_v_tangent_surface->point_on_surface(current_param[0], current_param[1], temp);
                mat.col(1) = temp;
                m_right_u_tangent_surface->point_on_surface(current_param[2], current_param[3], temp);
                mat.col(2) = -1.0 * temp;
                m_right_v_tangent_surface->point_on_surface(current_param[2], current_param[3], temp);
                mat.col(3) = -1.0 * temp;

                Eigen::JacobiSVD<Eigen::Matrix<T, dim, 4>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
                Eigen::Vector<T, 4> delta = matSvd.solve(vec);
                if (matSvd.info() != Eigen::Success)
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


                m_left_surface->point_on_surface(current_param[0], current_param[1], left_point);
                m_right_surface->point_on_surface(current_param[2], current_param[3], right_point);
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

        template<bool direction>
        ENUM_NURBS trace_point(Box<T, 4> &initial_box, const Eigen::Vector<T, 4> &param ,surf_surf_intersect_point<T, dim>& next_point, Box<T, 4>& next_initial_box)
        {
            //TODO:应在调用函数外部确保
            if (false == initial_box.intersect(m_product_box, initial_box))
                return ENUM_NURBS::NURBS_ERROR;
            int count1 = 0, count2 = 0;
            Eigen::Vector<T, 4> next_init_param;
            Box<T, 4> priori_enclosure;
            T param_space_eps = 2.0 * m_param_space_eps;
            Eigen::Vector<T, 4> intersect_param;
            priori_enclosure = initial_box;
            while (true)
            {
                if (ENUM_NURBS::NURBS_SUCCESS == eval_step_size_and_yBox<direction>(initial_box, param, priori_enclosure, next_init_param, param_space_eps))
                {
                    if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration(priori_enclosure, next_init_param, intersect_param))
                    {
                        break;
                    }
                    count2++;
                }
                count1++;
                if (count2 >= 2)
                {
                    std::cout << "count1: " << count1 << std::endl;
                    std::cout << "count2: " << count1 << std::endl;
                }
                
                param_space_eps /= 2.0;
                if (param_space_eps < TDEFAULT_ERROR<T>::value)
                {
                    return ENUM_NURBS::NURBS_ERROR;
                }
            }

            m_left_surface->point_on_surface(intersect_param[0], intersect_param[1], next_point.m_point);
            next_point.m_uv = intersect_param;

            next_initial_box.Min = intersect_param - m_box_size;
            next_initial_box.Max = intersect_param + m_box_size;
            
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        template<bool direction>
        ENUM_NURBS trace_curve_by_point(const Box<T, 4>& initial_box, surf_surf_int<surface_type>* intersection)
        {
            Box<T, 4> current_box = initial_box;
            Eigen::Vector<T, 4> initial_param = initial_box.get_middle_point();
            Box<T, 4> next_initial_box;

            surf_surf_intersect_point<T, dim>* current_intersection = new surf_surf_intersect_point<T, dim>();
            Eigen::Vector<T, dim> initial_point;
            m_left_surface->point_on_surface(initial_param[0], initial_param[1], initial_point);
            current_intersection->m_point = initial_point;
            current_intersection->m_uv = initial_param;
            intersection->m_intersect_points.push_back(current_intersection);

            int index = 0;

            for (; index < MAXINTERSETORPOINTNUMBER; ++index)
            {
                std::cout << index << std::endl;
                surf_surf_intersect_point<T, dim>* next_point = new surf_surf_intersect_point<T, dim>();
                if (trace_point<direction>(current_box, current_intersection->m_uv, *next_point, next_initial_box) != ENUM_NURBS::NURBS_SUCCESS)
                {
                    intersection->m_intersect_points_count.push_back(index);
                    return ENUM_NURBS::NURBS_SUCCESS;
                }
                current_box = next_initial_box;

                if ((current_intersection->m_uv - next_point->m_uv).norm() < 0.1 * TDEFAULT_ERROR<T>::value)
                {
                    intersection->m_intersect_points_count.push_back(index);
                    return ENUM_NURBS::NURBS_SUCCESS;
                }

                current_intersection->m_next = next_point;
                current_intersection = next_point;

            }
            if (index == MAXINTERSETORPOINTNUMBER)
            {
                //TODO:析构内存
                return ENUM_NURBS::NURBS_ERROR;
            }
            intersection->m_intersect_points_count.push_back(index);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        //暂时不处理相交为环的情况(目前只支持bezier情况)
        ENUM_NURBS surafces_intersection(surf_surf_int<surface_type>** intersection)
        {
            if (*intersection != nullptr)
            {
                return ENUM_NURBS::NURBS_ERROR;
            }
            *intersection = new surf_surf_int<surface_type>();

            //先取等参线和曲面求交(非闭合)
            //using curve_type = template surface_type::iso_curve_type;
            typename surface_type::iso_curve_type curve1;
            typename surface_type::iso_curve_type curve2;
            typename surface_type::iso_curve_type curve3;
            typename surface_type::iso_curve_type curve4;

            m_right_surface->get_isoparameter_curve<ENUM_DIRECTION::U_DIRECTION>(0.0, curve1);
            m_right_surface->get_isoparameter_curve<ENUM_DIRECTION::U_DIRECTION>(1.0, curve2);
            m_right_surface->get_isoparameter_curve<ENUM_DIRECTION::V_DIRECTION>(0.0, curve3);
            m_right_surface->get_isoparameter_curve<ENUM_DIRECTION::V_DIRECTION>(1.0, curve4);

            std::vector<Eigen::Vector<T, 3>> int_params1 = bezier_curve_int_bezier_surface(curve1, *m_left_surface);
            std::vector<Eigen::Vector<T, 3>> int_params2 = bezier_curve_int_bezier_surface(curve2, *m_left_surface);
            std::vector<Eigen::Vector<T, 3>> int_params3 = bezier_curve_int_bezier_surface(curve3, *m_left_surface);
            std::vector<Eigen::Vector<T, 3>> int_params4 = bezier_curve_int_bezier_surface(curve4, *m_left_surface);

            //TODO:去重, 加细; 另一个曲面的四条边界
            std::vector<Box<T, 4>> int_params;
            for (auto& param : int_params1)
            {
                Eigen::Vector<T, 4> initial_param = Eigen::Vector<T, 4>{ param[1], param[2], 0.0, param[0] };
                Box<T, 4> box;
                create_box(initial_param, TDEFAULT_ERROR<T>::value, box);
                int_params.push_back(box);
            }
            for (auto& param : int_params2)
            {
                Eigen::Vector<T, 4> initial_param = Eigen::Vector<T, 4>{ param[1], param[2], 1.0, param[0]};
                Box<T, 4> box;
                create_box(initial_param, TDEFAULT_ERROR<T>::value, box);
                int_params.push_back(box);
            }
            for (auto& param : int_params3)
            {
                Eigen::Vector<T, 4> initial_param = Eigen::Vector<T, 4>{ param[1], param[2], param[0], 0.0 };
                Box<T, 4> box;
                create_box(initial_param, TDEFAULT_ERROR<T>::value, box);
                int_params.push_back(box);
            }
            for (auto& param : int_params4)
            {
                Eigen::Vector<T, 4> initial_param = Eigen::Vector<T, 4>{ param[1], param[2], param[0], 1.0 };
                Box<T, 4> box;
                create_box(initial_param, TDEFAULT_ERROR<T>::value, box);
                int_params.push_back(box);
            }

            for (auto &box : int_params)
            {
                trace_curve_by_point<true>(box, *intersection);
                trace_curve_by_point<false>(box, *intersection);
            }
            return ENUM_NURBS::NURBS_SUCCESS;

        }


};

}