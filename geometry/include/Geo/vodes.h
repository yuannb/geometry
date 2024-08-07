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

    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS eval_point_interval(const surface_type* surf, const Box<T, 2>& interval, std::array<Interval<T>, (unsigned)dim>& box_interval)
    {
        surface_type* sub_surface = new surface_type();
        std::unique_ptr<surface_type> raii = std::unique_ptr<surface_type>(sub_surface);

        Box<T, 2> interval_copy = interval;
        surf->sub_divide(interval_copy, *sub_surface);

        Box<T, dim> sub_box;
        sub_surface->get_box(sub_box);

        //将u_box,v_box 转为interval vector
        for (int index = 0; index < dim; ++index)
        {
            box_interval[index].m_interval[0] = sub_box.Min[index];
            box_interval[index].m_interval[1] = sub_box.Max[index];
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS eval_point_interval(const surface_type* surf, const Box<T, 2>& interval, Box<T, dim> &box)
    {
        surface_type* sub_surface = new surface_type();
        std::unique_ptr<surface_type> raii = std::unique_ptr<surface_type>(sub_surface);

        Box<T, 2> interval_copy = interval;
        surf->sub_divide(interval_copy, *sub_surface);

        sub_surface->get_box(box);
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
                                        Eigen::Matrix<T, 4, 2> &param_ders, Eigen::Matrix<T, dim, 3> &space_ders, T &x1, T &x2)
    {
        static_assert(3 == dim, "3 != dim");

        Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> ders;
        left_surface->template derivative_on_surface<2, 2>(u, v, ders);
        space_ders.col(0) = ders(0, 0);
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
        x1 = param_ders.template block<2, 1>(0, 0).transpose() * left_first_fundenmental * param_ders.template block<2, 1>(0, 0);
        dot_ders[0] = tangent_vec.dot(ders2(1, 0));
        dot_ders[1] = tangent_vec.dot(ders2(0, 1));
        Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd2(right_first_fundenmental);
        param_ders.template block<2, 1>(2, 0) = matSvd2.solve(dot_ders);
        if (matSvd2.info() !=  Eigen::Success)
            return ENUM_NURBS::NURBS_ERROR;
        x2 = param_ders.template block<2, 1>(2, 0).transpose() * right_first_fundenmental * param_ders.template block<2, 1>(2, 0);
        //计算两曲面的交线的二阶微分(弧长参数下)
        
        //1.计算 d(u, v) 的第二基本型的值
        Eigen::Vector<T, 2> second_fundenmental_value;
        second_fundenmental_value[0] = param_ders.template block<2, 1>(0, 0).col(0).transpose() * left_second_fundenmental * param_ders.template block<2, 1>(0, 0).col(0);
        second_fundenmental_value[1] = param_ders.template block<2, 1>(2, 0).transpose() * right_second_fundenmental * param_ders.template block<2, 1>(2, 0);
        Eigen::Matrix2<T> mat;
        mat(0, 0) = mat(1, 1) = 1.0;
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
    

    //intersect points(以后有空再改吧)
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
        std::vector<surface_type*> m_left_dxx;
        
        //ders of right surface
        surface_type* m_right_u_tangent_surface;
        surface_type* m_right_v_tangent_surface;
        std::vector<surface_type*> m_right_dxx;

        surface_type* m_left_normal_surface;
        surface_type* m_right_normal_surface;

        T m_param_space_eps;
        T m_angle_eps;

        Eigen::Vector<T, 4> m_box_size;

        Box<T, 4> m_product_box;

        Interval<T> m_alpha;
        Interval<T> m_beta;
        Interval<T> m_u;
        Interval<T> m_v;

        Eigen::Vector<T, 4> m_current_param;

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
        ENUM_NURBS init(T space_eps, T angle_eps = TDEFAULTANGLETOL<T>::value)
        {
            m_left_u_tangent_surface = new surface_type();
            m_left_v_tangent_surface = new surface_type();
            m_right_u_tangent_surface = new surface_type();
            m_right_v_tangent_surface = new surface_type();
            m_left_normal_surface = new surface_type();
            m_right_normal_surface = new surface_type();
            for (int index = 0; index < 3; ++index)
            {
                m_left_dxx.push_back(new surface_type());
                m_right_dxx.push_back(new surface_type());
            }

            m_left_surface->tangent_u_surface(*m_left_u_tangent_surface);
            m_left_surface->tangent_v_surface(*m_left_v_tangent_surface);
            m_right_surface->tangent_u_surface(*m_right_u_tangent_surface);
            m_right_surface->tangent_v_surface(*m_right_v_tangent_surface);

            Box<T, dim> left_u_tangent_surface_box, left_v_tangent_surface_box, right_u_tangent_surface_box, right_v_tangent_surface_box;
            m_left_u_tangent_surface->get_box(left_u_tangent_surface_box);
            m_left_v_tangent_surface->get_box(left_v_tangent_surface_box);
            m_right_u_tangent_surface->get_box(right_u_tangent_surface_box);
            m_right_v_tangent_surface->get_box(right_v_tangent_surface_box);

            //Xuu
            m_left_u_tangent_surface->tangent_u_surface(*(m_left_dxx[0]));
            //Xuv = Xvu
            m_left_u_tangent_surface->tangent_v_surface(*(m_left_dxx[1]));
            //Xvv
            m_left_v_tangent_surface->tangent_v_surface(*(m_left_dxx[2]));


            //Yuu
            m_right_u_tangent_surface->tangent_u_surface(*(m_right_dxx[0]));
            //Yuv = Yvu
            m_right_u_tangent_surface->tangent_v_surface(*(m_right_dxx[1]));
            //Yvv
            m_right_v_tangent_surface->tangent_v_surface(*(m_right_dxx[2]));

            m_left_surface->eval_normal_surface(*m_left_normal_surface);
            m_right_surface->eval_normal_surface(*m_right_normal_surface);

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

            m_angle_eps = angle_eps;
            m_product_box = m_left_surface->get_interval().product_box(m_right_surface->get_interval());

            T box_size_num = std::min(SPACE_PARAM_TIMES<T>::value * PRECISION<T>::value, m_param_space_eps / 16.0);
            m_box_size.setConstant(box_size_num);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        // type = 0, 1, 2
        ENUM_NURBS eval_tangent(const Eigen::Vector<T, 4> &initial_param, Eigen::Vector<T, dim> &tangent, 
            Eigen::Vector<T, 4>& param_tangent,  int type)
        {
            Eigen::Vector<T, dim> Xu;
            m_left_u_tangent_surface->point_on_surface(initial_param[0], initial_param[1], Xu);
            Eigen::Vector<T, dim> Xv;
            m_left_v_tangent_surface->point_on_surface(initial_param[0], initial_param[1], Xv);

            Eigen::Vector<T, dim> Yu;
            m_right_u_tangent_surface->point_on_surface(initial_param[2], initial_param[3], Yu);
            Eigen::Vector<T, dim> Yv;
            m_right_v_tangent_surface->point_on_surface(initial_param[2], initial_param[3], Yv);
            Eigen::Vector<T, dim> Nx = Xu.cross(Xv);
            T Nx_len = Nx.norm();
            Nx.normalize();
            Eigen::Vector<T, dim> Ny = Yu.cross(Yv);
            T Ny_len = Ny.norm();
            Ny.normalize();
            if (type == 0)
            {
                tangent = Nx.cross(Ny);
                tangent.normalize();
                //return ENUM_NURBS::NURBS_SUCCESS;
            }
            else
            {
                Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_ders;
                m_left_surface->template derivative_on_surface<2>(initial_param[0], initial_param[1], left_ders);
                Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_ders;
                m_right_surface->template derivative_on_surface<2>(initial_param[2], initial_param[3], right_ders);

                Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
                left_normal.normalize();

                Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
                right_normal.normalize();
                T denominator = 1.0 / right_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal);
                T a11 = left_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal) * denominator;
                T a12 = left_ders(0, 1).cross(right_ders(0, 1)).dot(right_normal) * denominator;
                T a21 = right_ders(1, 0).cross(left_ders(1, 0)).dot(right_normal) * denominator;
                T a22 = right_ders(1, 0).cross(left_ders(0, 1)).dot(right_normal) * denominator;
                T left_L = left_ders(2, 0).dot(left_normal);
                T left_M = left_ders(1, 1).dot(left_normal);
                T left_N = left_ders(0, 2).dot(left_normal);

                T right_L = right_ders(2, 0).dot(right_normal);
                T right_M = right_ders(1, 1).dot(right_normal);
                T right_N = right_ders(0, 2).dot(right_normal);

                T b11 = a11 * a11 * right_L;
                b11 += 2 * a11 * a21 * right_M;
                b11 += a21 * a21 * right_N;
                b11 -= left_L;

                T b12 = a11 * a12 * right_L;
                b12 += (a11 * a22 + a21 * a12) * right_M;
                b12 += a21 * a22 * right_N;
                b12 -= left_M;

                T b22 = a12 * a12 * right_L;
                b22 += 2 * a12 * a22 * right_M;
                b22 += a22 * a22 * right_N;
                b22 -= left_N;
                //double delta = b12 * b12 - b11 * b22;
                if (type == 1)
                {
                    T nu = b12 / b11;
                    tangent = left_ders(0, 1) - nu * left_ders(1, 0);
                    tangent.normalize();
                    //return ENUM_NURBS::NURBS_SUCCESS;
                }
                else
                {
                    T mu = b12 / b22;
                    tangent = left_ders(1, 0) - mu * left_ders(0, 1);
                    tangent.normalize();
                    //return ENUM_NURBS::NURBS_SUCCESS;
                }
            }
                
            param_tangent[0] = tangent.cross(Xv).dot(Nx) / Nx_len;
            param_tangent[1] = Xu.cross(tangent).dot(Nx) / Nx_len;
            param_tangent[2] = tangent.cross(Yv).dot(Ny) / Ny_len;
            param_tangent[3] = Yu.cross(tangent).dot(Ny) / Ny_len;
            return ENUM_NURBS::NURBS_SUCCESS; 
        }
            

        //2阶微分的box太大，要么能够缩小2阶微分的box，要么只使用一阶微分的box
        //TODO: 相切的情形的一阶微分的box
        template<unsigned order>
        ENUM_NURBS eval_preiamge_and_space_ders(const Box<T, 4> &param_box, std::vector<std::array<Box<T, 4>, order>> &param_ders, 
            std::vector<std::array<Box<T, dim>, order>>&space_ders)
        {
            //目前仅仅支持三维
            static_assert(3 == dim, "dim != 3 is not supported");
            static_assert(order == 1 || order == 2, "order != 1 && order != 2");

            //计算交线的切向, 两个曲面非相切的时候;
            Box<T, dim> left_normal_box, right_normal_box;
            Box<T, dim> left_u_tangent_box, left_v_tangent_box;
            Box<T, dim> right_u_tangent_box, right_v_tangent_box;

            Box<T, 2> left_param_box(param_box.Min.template block<2, 1>(0, 0), param_box.Max.template block<2, 1>(0, 0));
            Box<T, 2> right_param_box(param_box.Min.template block<2, 1>(2, 0), param_box.Max.template block<2, 1>(2, 0));
            eval_point_interval(m_left_u_tangent_surface, left_param_box, left_u_tangent_box);
            eval_point_interval(m_left_v_tangent_surface, left_param_box, left_v_tangent_box);
            eval_point_interval(m_right_u_tangent_surface, right_param_box, right_u_tangent_box);
            eval_point_interval(m_right_v_tangent_surface, right_param_box, right_v_tangent_box);

            //使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
            eval_point_interval(m_left_normal_surface, left_param_box, left_normal_box);
            eval_point_interval(m_right_normal_surface, right_param_box, right_normal_box);

            Eigen::Matrix2<Box<T, dim>> left_dxx;
            eval_point_interval(m_left_dxx[0], left_param_box, left_dxx(0, 0));
            eval_point_interval(m_left_dxx[1], left_param_box, left_dxx(1, 0));
            left_dxx(0, 1) = left_dxx(1, 0);
            eval_point_interval(m_left_dxx[2], left_param_box, left_dxx(1, 1));

            Eigen::Matrix2<Box<T, dim>> right_dxx;
            eval_point_interval(m_right_dxx[0], right_param_box, right_dxx(0, 0));
            eval_point_interval(m_right_dxx[1], right_param_box, right_dxx(1, 0));
            right_dxx(0, 1) = right_dxx(1, 0);
            eval_point_interval(m_right_dxx[2], right_param_box, right_dxx(1, 1));

            Box<T, dim> c;
            c = interval_algorithm::cross(left_normal_box, right_normal_box);
            Box<T, 1> left_normal_len2 = interval_algorithm::dot(left_normal_box, left_normal_box);
            Box<T, 1> right_normal_len2 = interval_algorithm::dot(right_normal_box, right_normal_box);

            Box<T, dim> unit_c;
            ENUM_NURBS error_code = interval_algorithm::normalized(c, unit_c);
            if (error_code != ENUM_NURBS::NURBS_SUCCESS)
                return error_code;

            space_ders.resize(1);
            space_ders[0][0] = c;
            param_ders.resize(1);

            std::vector<Box<T, 1>> alpha_d(2);
            std::vector<Box<T, 1>> beta_d(2);
            if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(unit_c, left_v_tangent_box)), left_normal_len2, alpha_d[0]))
                return ENUM_NURBS::NURBS_ERROR;

            if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(left_u_tangent_box, unit_c)), left_normal_len2, alpha_d[1]))
                return ENUM_NURBS::NURBS_ERROR;

            if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(unit_c, right_v_tangent_box)), right_normal_len2, beta_d[0]))
                return ENUM_NURBS::NURBS_ERROR;

            if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(right_u_tangent_box, unit_c)), right_normal_len2, beta_d[1]))
                return ENUM_NURBS::NURBS_ERROR;
            param_ders[0][0].set_index_interval(0, alpha_d[0]);
            param_ders[0][0].set_index_interval(1, alpha_d[1]);
            param_ders[0][0].set_index_interval(2, beta_d[0]);
            param_ders[0][0].set_index_interval(3, beta_d[1]);
            if constexpr (order == 2)
            {
                Box<T, dim> left_L;
                left_L.Min.setConstant(0.0);
                left_L.Max.setConstant(0.0);

                Box<T, dim> right_L;
                right_L.Min.setConstant(0.0);
                right_L.Max.setConstant(0.0);
                for (int i = 0; i < 2; ++i)
                {
                    for (int j = 0; j < 2; ++j)
                    {
                        left_L = left_L + left_dxx(i, j) * (alpha_d[i] * alpha_d[j]);
                        right_L = right_L + right_dxx(i, j) * (beta_d[i] * beta_d[j]);
                    }
                }

                Eigen::Matrix2<Box<T, 1>> left_first_fundenmental;
                left_first_fundenmental(0, 0) = interval_algorithm::dot(left_u_tangent_box, left_u_tangent_box);
                left_first_fundenmental(1, 0) = left_first_fundenmental(0, 1) = interval_algorithm::dot(left_u_tangent_box, left_v_tangent_box);
                left_first_fundenmental(1, 1) = interval_algorithm::dot(left_v_tangent_box, left_v_tangent_box);
                Box<T, 1> left_first_det = left_first_fundenmental.determinant();

                Eigen::Matrix2<Box<T, 1>> right_first_fundenmental;
                right_first_fundenmental(0, 0) = interval_algorithm::dot(right_u_tangent_box, right_u_tangent_box);
                right_first_fundenmental(1, 0) = right_first_fundenmental(0, 1) = interval_algorithm::dot(right_u_tangent_box, right_v_tangent_box);
                right_first_fundenmental(1, 1) = interval_algorithm::dot(right_v_tangent_box, right_v_tangent_box);
                Box<T, 1> right_first_det = right_first_fundenmental.determinant();


                Eigen::Vector<Box<T, 1>, 2> second_fundenmental_value;
                second_fundenmental_value[0] = interval_algorithm::dot(left_L, left_normal_box);
                second_fundenmental_value[1] = interval_algorithm::dot(right_L, right_normal_box);
                Eigen::Matrix2<Box<T, 1>> mat;
                mat(0, 0) = left_normal_len2;
                mat(1, 1) = right_normal_len2;
                mat(1, 0) = mat(0, 1) = interval_algorithm::dot(left_normal_box, right_normal_box);
                Box<T, 1> mat_det = mat.determinant();

                Eigen::Matrix2<Box<T, 1>> b1 = mat;
                b1.col(0) = second_fundenmental_value;
                Box<T, 1> coeff1;
                if (false == interval_algorithm::divide(b1.determinant(), mat_det, coeff1))
                    return ENUM_NURBS::NURBS_ERROR;

                Eigen::Matrix2<Box<T, 1>> b2 = mat;
                b2.col(1) = second_fundenmental_value;
                Box<T, 1> coeff2;
                if (false == interval_algorithm::divide(b2.determinant(), mat_det, coeff2))
                    return ENUM_NURBS::NURBS_ERROR;

                Box<T, dim> c_dd = left_normal_box * coeff1 + right_normal_box * coeff2;
                //交线的二阶导数
                space_ders[0][1] = c_dd;

                Box<T, dim> temp = c_dd - left_L;
                Eigen::Vector2<Box<T, 1>> dot_ders;
                dot_ders[0] = interval_algorithm::dot(temp, left_u_tangent_box);
                dot_ders[1] = interval_algorithm::dot(temp, left_v_tangent_box);

                b1 = left_first_fundenmental;
                b1.col(0) = dot_ders;
                if (false == interval_algorithm::divide(b1.determinant(), left_first_det, coeff1))
                    return ENUM_NURBS::NURBS_ERROR;
                b2 = left_first_fundenmental;
                b2.col(1) = dot_ders;
                if (false == interval_algorithm::divide(b2.determinant(), left_first_det, coeff2))
                    return ENUM_NURBS::NURBS_ERROR;

                param_ders[0][1].set_index_interval(0, coeff1);
                param_ders[0][1].set_index_interval(1, coeff2);

                temp = c_dd - right_L;
                dot_ders[0] = interval_algorithm::dot(temp, right_u_tangent_box);
                dot_ders[1] = interval_algorithm::dot(temp, right_v_tangent_box);

                b1 = right_first_fundenmental;
                b1.col(0) = dot_ders;
                if (false == interval_algorithm::divide(b1.determinant(), right_first_det, coeff1))
                    return ENUM_NURBS::NURBS_ERROR;
                b2 = right_first_fundenmental;
                b2.col(1) = dot_ders;
                if (false == interval_algorithm::divide(b2.determinant(), right_first_det, coeff2))
                    return ENUM_NURBS::NURBS_ERROR;

                param_ders[0][1].set_index_interval(2, coeff1);
                param_ders[0][1].set_index_interval(3, coeff2);
            }
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
        template<int order>
        ENUM_NURBS eval_preiamge_and_space_ders(const Eigen::Vector<T, 4> &param, std::vector<Eigen::Matrix<T, 4, order>>& param_ders, std::vector<Eigen::Matrix<T, dim, order>>& space_ders, int type)
        {
            static_assert(3 == dim, "3 != dim");

            param_ders.clear();
            space_ders.clear();
            Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_ders;
            m_left_surface->template derivative_on_surface<2, 2>(param[0], param[1], left_ders);
            Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
            left_normal.normalize();
            Eigen::Matrix2<T> left_first_fundenmental;
            Eigen::Matrix2<T> left_second_fundenmental;
            eval_first_and_second_fundenmental(left_ders, left_normal, left_first_fundenmental, left_second_fundenmental);

            Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_ders;
            m_right_surface->template derivative_on_surface<2, 2>(param[2], param[3], right_ders);
            Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
            right_normal.normalize();
            Eigen::Matrix2<T> right_first_fundenmental;
            Eigen::Matrix2<T> right_second_fundenmental;
            eval_first_and_second_fundenmental(right_ders, right_normal, right_first_fundenmental, right_second_fundenmental);

            Eigen::Vector<T, dim> tangent;
            if (type == 0)
            {
                tangent = left_normal.cross(right_normal);
                Eigen::Matrix<T, dim, order> ders;
                ders.col(0) = tangent.normalized();
                space_ders.push_back(ders);
            }
            else
            {
                //Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_ders;
                //m_left_surface->template derivative_on_surface<2>(initial_param[0], initial_param[1], left_ders);
                //Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_ders;
                //m_right_surface->template derivative_on_surface<2>(initial_param[2], initial_param[3], right_ders);

                //Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
                //left_normal.normalize();

                //Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
                //right_normal.normalize();
                int flag = right_normal.dot(left_normal) > 0 ? 1 : -1;
                T denominator = 1.0 / right_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal);
                T a11 = left_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal) * denominator;
                T a12 = left_ders(0, 1).cross(right_ders(0, 1)).dot(right_normal) * denominator;
                T a21 = right_ders(1, 0).cross(left_ders(1, 0)).dot(right_normal) * denominator;
                T a22 = right_ders(1, 0).cross(left_ders(0, 1)).dot(right_normal) * denominator;
                T left_L = left_ders(2, 0).dot(left_normal) * flag;
                T left_M = left_ders(1, 1).dot(left_normal) * flag;
                T left_N = left_ders(0, 2).dot(left_normal) * flag;

                T right_L = right_ders(2, 0).dot(right_normal);
                T right_M = right_ders(1, 1).dot(right_normal);
                T right_N = right_ders(0, 2).dot(right_normal);

                T b11 = a11 * a11 * right_L;
                b11 += 2 * a11 * a21 * right_M;
                b11 += a21 * a21 * right_N;
                b11 -= left_L;

                T b12 = a11 * a12 * right_L;
                b12 += (a11 * a22 + a21 * a12) * right_M;
                b12 += a21 * a22 * right_N;
                b12 -= left_M;

                T b22 = a12 * a12 * right_L;
                b22 += 2 * a12 * a22 * right_M;
                b22 += a22 * a22 * right_N;
                b22 -= left_N;
                double delta = 4.0 * (b12 * b12 - b11 * b22);
                bool is_branch = false;
                if (type == 3)
                {
                    if (delta < -KNOTS_EPS<T>::value)
                    {
                        return ENUM_NURBS::NURBS_ISOLATED_TANGENTIAL_POINT;
                    }
                    if (delta > KNOTS_EPS<T>::value)
                    {
                        is_branch = true;
                    }
                    if (std::abs(b11) < KNOTS_EPS<T>::value && std::abs(b12) < KNOTS_EPS<T>::value
                        && std::abs(b22) < KNOTS_EPS<T>::value)
                    {
                        return ENUM_NURBS::NURBS_HIGH_ORDER_TANGENTIAL;
                    }
                }

                if (type == 1 || (type == 3 && std::abs(b11) >= std::abs(b22)))
                {
                    std::vector<T> coeff;
                    if (is_branch == true)
                    {
                        T delta_r = std::sqrt(delta);
                        coeff.push_back(2 * b12 + delta_r);
                        coeff.push_back(2 * b12 - delta_r);
                    }
                    else
                    {
                        coeff.push_back(2 * b12);
                    }
                    for (const T& coef : coeff)
                    {
                        T nu = coef / (2 * b11);
                        tangent = left_ders(0, 1) - nu * left_ders(1, 0);
                        Eigen::Matrix<T, dim, order> ders;
                        ders.col(0) = tangent.normalized();
                        space_ders.push_back(ders);
                    }
                }
                if (type == 2 || (type == 3 && std::abs(b11) < std::abs(b22)))
                {
                    std::vector<T> coeff;
                    if (is_branch == true)
                    {
                        T delta_r = std::sqrt(delta);
                        coeff.push_back(2 * b12 + delta_r);
                        coeff.push_back(2 * b12 - delta_r);
                    }
                    else
                    {
                        coeff.push_back(2 * b12);
                    }
                    for (const T& coef : coeff)
                    {
                        T mu = coef / (2 * b22);
                        tangent = left_ders(1, 0) - mu * left_ders(0, 1);
                        Eigen::Matrix<T, dim, order> ders;
                        ders.col(0) = tangent.normalized();
                        space_ders.push_back(ders);
                    }
                    //return ENUM_NURBS::NURBS_SUCCESS;
                }
            }

            Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(left_first_fundenmental);
            Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd2(right_first_fundenmental);
            Eigen::Matrix2<T> mat;
            mat(0, 0) = mat(1, 1) = 1.0;
            mat(0, 1) = mat(1, 0) = left_normal.dot(right_normal);
            Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd3(mat);
            
            for (Eigen::Matrix<T, dim, order> &space_der : space_ders)
            {
                Eigen::Vector<T, dim> tangent_vec = space_der.col(0);
                Eigen::Vector<T, 2> dot_ders;
                dot_ders[0] = tangent_vec.dot(left_ders(1, 0));
                dot_ders[1] = tangent_vec.dot(left_ders(0, 1));

                Eigen::Matrix<T, 4, order> param_der;
                param_der.template block<2, 1>(0, 0) = matSvd.solve(dot_ders);
                if (matSvd.info() != Eigen::Success)
                    return ENUM_NURBS::NURBS_ERROR;
                dot_ders[0] = tangent_vec.dot(right_ders(1, 0));
                dot_ders[1] = tangent_vec.dot(right_ders(0, 1));

                param_der.template block<2, 1>(2, 0) = matSvd2.solve(dot_ders);
                if (matSvd2.info() != Eigen::Success)
                    return ENUM_NURBS::NURBS_ERROR;
                //计算两曲面的交线的二阶微分(弧长参数下)

                if constexpr (order == 2)
                {
                    //1.计算 d(u, v) 的第二基本型的值
                    Eigen::Vector<T, 2> second_fundenmental_value;
                    second_fundenmental_value[0] = param_der.template block<2, 1>(0, 0).col(0).transpose() * left_second_fundenmental * param_der.template block<2, 1>(0, 0).col(0);
                    second_fundenmental_value[1] = param_der.template block<2, 1>(2, 0).transpose() * right_second_fundenmental * param_der.template block<2, 1>(2, 0);
                    
                    Eigen::Vector2<T> coeff = matSvd3.solve(second_fundenmental_value);
                    if (matSvd3.info() != Eigen::Success)
                        return ENUM_NURBS::NURBS_ERROR;

                    //交线的二阶导数
                    Eigen::Vector<T, dim> alpha_dd = coeff[0] * left_normal + coeff[1] * right_normal;
                    space_der.col(1) = alpha_dd;
                    Eigen::Vector<T, dim> L = std::pow(param_der(0, 0), 2) * left_ders(2, 0) + 2 * param_der(0, 0) * param_der(1, 0) * left_ders(1, 1) + std::pow(param_der(1, 0), 2) * left_ders(0, 2);

                    Eigen::Vector<T, dim> temp = alpha_dd - L;
                    dot_ders[0] = temp.dot(left_ders(1, 0));
                    dot_ders[1] = temp.dot(left_ders(0, 1));
                    param_der.template block<2, 1>(0, 1) = matSvd.solve(dot_ders);
                    if (matSvd.info() != Eigen::Success)
                        return ENUM_NURBS::NURBS_ERROR;

                    L = std::pow(param_der(2, 0), 2) * right_ders(2, 0) + 2 * param_der(2, 0) * param_der(3, 0) * right_ders(1, 1) + std::pow(param_der(3, 0), 2) * right_ders(0, 2);
                    temp = alpha_dd - L;
                    dot_ders[0] = temp.dot(right_ders(1, 0));
                    dot_ders[1] = temp.dot(right_ders(0, 1));
                    param_der.template block<2, 1>(2, 1) = matSvd2.solve(dot_ders);
                    if (matSvd2.info() != Eigen::Success)
                        return ENUM_NURBS::NURBS_ERROR;
                }
                param_ders.push_back(param_der);
            }

            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS estimate_next_param(const Eigen::Vector<T, 4>& initial_param, const Box<T, 4> domain, T &step, Eigen::Vector<T, 4>& next_param, int type)
        {
            while (std::abs(step)> TDEFAULT_ERROR<T>::value)
            {
                //domain = Box<T, 4>(Eigen::Vector<T, 4>(0, 0, 0, 0), Eigen::Vector<T, 4>(1, 1, 1, 1));
                Eigen::Vector<T, dim> tangent;
                Eigen::Vector<T, 4> k1;
                //Eigen::Vector<T, 4> k11;
                //Eigen::Vector<T, 4> k21;
                //Eigen::Vector<T, 4> k31;
                //Eigen::Vector<T, 4> k41;

                std::vector<Eigen::Matrix<T, 4, 1>> param_ders;
                std::vector<Eigen::Matrix<T, dim, 1>> space_ders;
                eval_preiamge_and_space_ders<1>(initial_param, param_ders, space_ders, type);
                k1 = step * param_ders[0];
                //eval_tangent(initial_param, tangent, k11, type);
                //k11 -= param_ders[0];
                //T dis1 = k11.norm();
                //k1 *= step;

                Eigen::Vector<T, 4> mid_param = initial_param + 0.5 * k1;
                if (domain.is_contain_point(mid_param) == false)
                {
                    step /= 1.2;
                    continue;
                    //return ENUM_NURBS::NURBS_ERROR;
                }
                Eigen::Vector<T, 4> k2;
                eval_preiamge_and_space_ders<1>(mid_param, param_ders, space_ders, type);
                k2 = step * param_ders[0];
                //eval_tangent(mid_param, tangent, k21, type);
                //k21 -= param_ders[0];
                //T dis2 = k21.norm();
                //k2 *= step;
                mid_param = initial_param + 0.5 * k2;
                if (domain.is_contain_point(mid_param) == false)
                {
                    step /= 1.2;
                    continue;
                    //return ENUM_NURBS::NURBS_ERROR;
                }
                Eigen::Vector<T, 4> k3;
                eval_preiamge_and_space_ders<1>(mid_param, param_ders, space_ders, type);
                k3 = step * param_ders[0];
                //eval_tangent(mid_param, tangent, k31, type);
                //k31 -= param_ders[0];
                //T dis3 = k31.norm();
                //k3 *= step;
                mid_param = initial_param + k3;
                if (domain.is_contain_point(mid_param) == false)
                {
                    step /= 1.2;
                    continue;
                    //return ENUM_NURBS::NURBS_ERROR;
                }
                Eigen::Vector<T, 4> k4;
                eval_preiamge_and_space_ders<1>(mid_param, param_ders, space_ders, type);
                k4 = step * param_ders[0];
                //eval_tangent(mid_param, tangent, k41, type);
                //k41 -= param_ders[0];
                //T dis4 = k41.norm();
                //k4 *= step;

                next_param = initial_param + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
                if (domain.is_contain_point(next_param) == false)
                {
                    step /= 1.2;
                    continue;
                    //return ENUM_NURBS::NURBS_ERROR;
                }
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            return ENUM_NURBS::NURBS_ERROR;
        }

        ENUM_NURBS eval_priori_enclosure_inner2(const Box<T, 4>& initial_box, const Box<T, 4>& param_box, T step_size, Box<T, 4>& priori_enclosure, T& angle_diff, int& type, bool transversal = true)
        {
            std::vector<std::array<Box<T, 4>, 1>> param_ders;
            std::vector<std::array<Box<T, dim>, 1>> space_ders;
            ENUM_NURBS error_code = eval_preiamge_and_space_ders<1>(param_box, param_ders, space_ders);
            if (ENUM_NURBS::NURBS_SUCCESS != error_code)
            {
                return error_code;
            }

            Interval<T> step_box(0, step_size);

            Box<T, 4> box_interval = param_ders[0][0] * step_box;
            Box<T, 4> temp_priori_enclosure = initial_box + box_interval;

            if (true == m_product_box.is_contain_box(temp_priori_enclosure, 5.0 * PRECISION<T>::value))
            {
                temp_priori_enclosure.intersect(m_product_box, priori_enclosure);
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            return ENUM_NURBS::NURBS_ERROR;


        }

        ENUM_NURBS eval_priori_enclosure_inner(const Box<T, 4>& initial_box, const Box<T, 4>& param_box, T step_size, Box<T, 4>& priori_enclosure, T &angle_diff, int &type, bool transversal = true)
        {
            //目前仅仅支持三维
            static_assert(3 == dim, "dim != 3 is not supported");

            //计算交线的切向, 两个曲面非相切的时候; TODO:相切的时候
            std::array<Interval<T>, (unsigned)dim> left_normal_interval, right_normal_interval;
            std::array<Interval<T>, (unsigned)dim> left_normal_interval_t, right_normal_interval_t;
            std::array<Interval<T>, (unsigned)dim> left_u_tangent_interval, left_v_tangent_interval;
            std::array<Interval<T>, (unsigned)dim> right_u_tangent_interval, right_v_tangent_interval;

            Box<T, 2> left_param_box(param_box.Min.template block<2, 1>(0, 0), param_box.Max.template block<2, 1>(0, 0));
            Box<T, 2> right_param_box(param_box.Min.template block<2, 1>(2, 0), param_box.Max.template block<2, 1>(2, 0));
            eval_point_interval(m_left_u_tangent_surface, left_param_box, left_u_tangent_interval);
            eval_point_interval(m_left_v_tangent_surface, left_param_box, left_v_tangent_interval);
            eval_point_interval(m_right_u_tangent_surface, right_param_box, right_u_tangent_interval);
            eval_point_interval(m_right_v_tangent_surface, right_param_box, right_v_tangent_interval);
            eval_point_interval(m_left_normal_surface, left_param_box, left_normal_interval);
            eval_point_interval(m_right_normal_surface, right_param_box, right_normal_interval);


            //eval_normal_and_tangent_interval<surface_type>(m_left_u_tangent_surface, m_left_v_tangent_surface, left_param_box, left_normal_interval, left_u_tangent_interval, left_v_tangent_interval);
            //eval_normal_and_tangent_interval<surface_type>(m_right_u_tangent_surface, m_right_v_tangent_surface, right_param_box, right_normal_interval, right_u_tangent_interval, right_v_tangent_interval);
            Box<T, dim> u_normal_box, v_normal_box, left_u_tangent_box, left_v_tangent_box, right_u_tangent_box, right_v_tangent_box;
            for (int index = 0; index < dim; ++index)
            {
                u_normal_box.Min[index] = left_normal_interval[index].get_low();
                u_normal_box.Max[index] = left_normal_interval[index].get_high();
                v_normal_box.Min[index] = right_normal_interval[index].get_low();
                v_normal_box.Max[index] = right_normal_interval[index].get_high();

                left_u_tangent_box.Min[index] = left_u_tangent_interval[index].get_low();
                left_u_tangent_box.Max[index] = left_u_tangent_interval[index].get_high();

                left_v_tangent_box.Min[index] = left_v_tangent_interval[index].get_low();
                left_v_tangent_box.Max[index] = left_v_tangent_interval[index].get_high();

                right_u_tangent_box.Min[index] = right_u_tangent_interval[index].get_low();
                right_u_tangent_box.Max[index] = right_u_tangent_interval[index].get_high();

                right_v_tangent_box.Min[index] = right_v_tangent_interval[index].get_low();
                right_v_tangent_box.Max[index] = right_v_tangent_interval[index].get_high();
            }

            T u_minLen = u_normal_box.eval_minimal_distance(Eigen::Vector3d(0, 0, 0));
            T u_maxLen = u_normal_box.eval_maximal_distance(Eigen::Vector3d(0, 0, 0));

            T v_minLen = v_normal_box.eval_minimal_distance(Eigen::Vector3d(0, 0, 0));
            T v_maxLen = v_normal_box.eval_maximal_distance(Eigen::Vector3d(0, 0, 0));
            Interval<T> u_normal_dot(u_minLen * u_minLen, u_maxLen * u_maxLen);
            Interval<T> v_normal_dot(v_minLen * v_minLen, v_maxLen * v_maxLen);

            // Transversal Intersection
            std::array<Interval<T>, (unsigned)dim> c;
            if (transversal == true)
            {
                c = vector_product(left_normal_interval, right_normal_interval);
                type = 0;
            }
            else
            {
                //Eigen::Vector<T, 4> initial_param = initial_box.get_middle_point();
                //Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_ders;
                //m_left_surface->template derivative_on_surface<2>(initial_param[0], initial_param[1], left_ders);
                //Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_ders;
                //m_right_surface->template derivative_on_surface<2>(initial_param[2], initial_param[3], right_ders);

                //Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
                //left_normal.normalize();

                //Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
                //right_normal.normalize();
                //T de = right_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal);
                //T a11_t = left_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal) / de;
                //T a12_t = left_ders(0, 1).cross(right_ders(0, 1)).dot(right_normal) / de;
                //T a21_t = right_ders(1, 0).cross(left_ders(1, 0)).dot(right_normal) / de;
                //T a22_t = right_ders(1, 0).cross(left_ders(0, 1)).dot(right_normal) / de;
                //T left_L_t = left_ders(2, 0).dot(left_normal);
                //T left_M_t = left_ders(1, 1).dot(left_normal);
                //T left_N_t = left_ders(0, 2).dot(left_normal);

                //T right_L_t = right_ders(2, 0).dot(right_normal);
                //T right_M_t = right_ders(1, 1).dot(right_normal);
                //T right_N_t = right_ders(0, 2).dot(right_normal);

                //T b11_t = a11_t * a11_t * right_L_t;
                //b11_t += 2 * a11_t * a21_t * right_M_t;
                //b11_t += a21_t * a21_t * right_N_t;
                //b11_t -= left_L_t;

                //T b12_t = a11_t * a12_t * right_L_t;
                //b12_t += (a11_t * a22_t + a21_t * a12_t) * right_M_t;
                //b12_t += a21_t * a22_t * right_N_t;
                //b12_t -= left_M_t;

                //T b22_t = a12_t * a12_t * right_L_t;
                //b22_t += 2 * a12_t * a22_t * right_M_t;
                //b22_t += a22_t * a22_t * right_N_t;
                //b22_t -= left_N_t;
                //double delta_t = b12_t * b12_t - b11_t * b22_t;


                // Tangential Intersection
                Interval<T> noraml_dot = dot_product(left_normal_interval, right_normal_interval);
                if (noraml_dot.contain(0.0, PRECISION<T>::value) == true)
                {
                    return ENUM_NURBS::NURBS_ERROR;
                }
                int flag = noraml_dot.m_interval[0] > 0 ? 1 : -1;
                Box<T, dim> left_normal_box, right_normal_box;
                std::array<Interval<T>, (unsigned)dim> left_normal_interval_t = left_normal_interval;
                std::array<Interval<T>, (unsigned)dim> right_normal_interval_t = right_normal_interval;
                for (int index = 0; index < dim; ++index)
                {
                    left_normal_box.Min[index] = left_normal_interval_t[index].get_low();
                    left_normal_box.Max[index] = left_normal_interval_t[index].get_high();

                    right_normal_box.Min[index] = right_normal_interval_t[index].get_low();
                    right_normal_box.Max[index] = right_normal_interval_t[index].get_high();
                }

                T minLen = left_normal_box.eval_minimal_distance(Eigen::Vector3d(0, 0, 0));
                T maxLen = left_normal_box.eval_maximal_distance(Eigen::Vector3d(0, 0, 0));

                Interval<T> left_normal_dot(minLen, maxLen);

                if (left_normal_dot.m_interval[0] <= TDEFAULT_ERROR<T>::value * TDEFAULT_ERROR<T>::value)
                    return ENUM_NURBS::NURBS_ERROR;

                if (divide<T, (unsigned)dim>(left_normal_interval_t, left_normal_dot, left_normal_interval_t) == false)
                    return ENUM_NURBS::NURBS_ERROR;

                minLen = right_normal_box.eval_minimal_distance(Eigen::Vector3d(0, 0, 0));
                maxLen = right_normal_box.eval_maximal_distance(Eigen::Vector3d(0, 0, 0));

                Interval<T> right_normal_dot(minLen, maxLen);
                if (right_normal_dot.m_interval[0] <= TDEFAULT_ERROR<T>::value * TDEFAULT_ERROR<T>::value)
                    return ENUM_NURBS::NURBS_ERROR;

                if (divide<T, (unsigned)dim>(right_normal_interval_t, right_normal_dot, right_normal_interval_t) == false)
                    return ENUM_NURBS::NURBS_ERROR;


                // TODO : 缩小下面量的范围，估计更准确一些
                //Interval<T> denominator = dot_product(right_normal_interval, right_normal_interval);
                Interval<T> denominator_r;
                divide(Interval<T>(1.0, 1.0), right_normal_dot, denominator_r);
                Interval<T> a11 = mutilply<T>(dot_product<T>(vector_product<T>(left_u_tangent_interval, right_v_tangent_interval), right_normal_interval_t), denominator_r);
                Interval<T> a12 = mutilply<T>(dot_product<T>(vector_product<T>(left_v_tangent_interval, right_v_tangent_interval), right_normal_interval_t), denominator_r);
                Interval<T> a21 = mutilply<T>(dot_product<T>(vector_product<T>(right_u_tangent_interval, left_u_tangent_interval), right_normal_interval_t), denominator_r);
                Interval<T> a22 = mutilply<T>(dot_product<T>(vector_product<T>(right_u_tangent_interval, left_v_tangent_interval), right_normal_interval_t), denominator_r);
                std::array<Interval<T>, (unsigned)dim> Xuu_box, Xuv_box, Xvv_box, Yuu_box, Yuv_box, Yvv_box;
                eval_point_interval<surface_type>(m_left_dxx[0], left_param_box, Xuu_box);
                eval_point_interval<surface_type>(m_left_dxx[1], left_param_box, Xuv_box);
                eval_point_interval<surface_type>(m_left_dxx[2], left_param_box, Xvv_box);

                eval_point_interval<surface_type>(m_right_dxx[0], right_param_box, Yuu_box);
                eval_point_interval<surface_type>(m_right_dxx[1], right_param_box, Yuv_box);
                eval_point_interval<surface_type>(m_right_dxx[2], right_param_box, Yvv_box);
                Interval<T> left_L = dot_product<T>(Xuu_box, left_normal_interval_t);
                Interval<T> left_M = dot_product<T>(Xuv_box, left_normal_interval_t);
                Interval<T> left_N = dot_product<T>(Xvv_box, left_normal_interval_t);
                if (flag == -1)
                {
                    Interval<T> Flag = Interval<T>(flag, flag);
                    left_L = mutilply(left_L, Flag);
                    left_M = mutilply(left_M, Flag);
                    left_N = mutilply(left_N, Flag);
                }

                Interval<T> right_L = dot_product<T>(Yuu_box, right_normal_interval_t);
                Interval<T> right_M = dot_product<T>(Yuv_box, right_normal_interval_t);
                Interval<T> right_N = dot_product<T>(Yvv_box, right_normal_interval_t);

                Interval<T> temp = mutilply<T>(a11, a11);
                temp.m_interval[0] = 0 > temp.m_interval[0] ? 0 : temp.m_interval[0];
                Interval<T> b11 = mutilply(temp, right_L);
                b11 = plus<T>(b11, mutilply<T>(Interval<T>(2.0, 2.0), mutilply<T>(mutilply<T>(a11, a21), right_M)));

                temp = mutilply<T>(a21, a21);
                temp.m_interval[0] = 0 > temp.m_interval[0] ? 0 : temp.m_interval[0];
                b11 = plus<T>(b11, mutilply(temp, right_N));
                b11 = minus<T>(b11, left_L);

                Interval<T> b12 = mutilply(mutilply<T>(a11, a12), right_L);
                b12 = plus<T>(b12, mutilply<T>(right_M, plus<T>(mutilply<T>(a11, a22), mutilply<T>(a21, a12))));
                b12 = plus<T>(b12, mutilply(mutilply<T>(a21, a22), right_N));
                b12 = minus<T>(b12, left_M);

                temp = mutilply<T>(a12, a12);
                temp.m_interval[0] = 0 > temp.m_interval[0] ? 0 : temp.m_interval[0];
                Interval<T> b22 = mutilply(temp, right_L);
                b22 = plus<T>(b22, mutilply<T>(Interval<T>(2.0, 2.0), mutilply<T>(mutilply<T>(a12, a22), right_M)));

                temp = mutilply<T>(a22, a22);
                temp.m_interval[0] = 0 > temp.m_interval[0] ? 0 : temp.m_interval[0];
                b22 = plus<T>(b22, mutilply(temp, right_N));
                b22 = minus<T>(b22, left_N);

                temp = mutilply(b12, b12);
                temp.m_interval[0] = 0 > temp.m_interval[0] ? 0 : temp.m_interval[0];
                Interval<T> delta = minus(temp, mutilply(b11, b22));

                if (delta.get_low() > 0 || delta.get_high() < 0)
                {
                    return ENUM_NURBS::NURBS_ERROR;
                }
                if (b11.contain(0, 0) && b12.contain(0, 0) && b22.contain(0, 0))
                {
                    return ENUM_NURBS::NURBS_ERROR;
                }
                if (b11.contain(0, 0) == false)
                {
                    Interval<T> nu;
                    divide(b12, b11, nu);
                    c = mutilply<T, dim>(left_u_tangent_interval, nu);
                    c = minus<T, dim>(left_v_tangent_interval, c);
                    type = 1;
                }
                else if (b22.contain(0, 0) == false)
                {
                    Interval<T> mu;
                    divide(b12, b22, mu);
                    c = mutilply<T, dim>(left_v_tangent_interval, mu);
                    c = minus<T, dim>(left_u_tangent_interval, c);
                    type = 2;
                }
                //else
                //{
                //    return ENUM_NURBS::NURBS_ERROR
                //}
            }


            Box<T, dim> c_box;
            for (int index = 0; index < dim; ++index)
            {
                c_box.Min[index] = c[index].get_low();
                c_box.Max[index] = c[index].get_high();
            }
            Box<T, dim> un_box, vn_box;
            for (int index = 0; index < dim; ++index)
            {
                un_box.Min[index] = left_normal_interval[index].get_low();
                un_box.Max[index] = left_normal_interval[index].get_high();

                vn_box.Min[index] = right_normal_interval[index].get_low();
                vn_box.Max[index] = right_normal_interval[index].get_high();
            }
            cone<T, dim> un_cone = point_box(un_box, Eigen::Vector3d(0, 0, 0));
            cone<T, dim> vn_cone = point_box(vn_box, Eigen::Vector3d(0, 0, 0));
            angle_diff = std::max(un_cone.m_angle * 2.0, vn_cone.m_angle * 2.0);

            cone<T, dim> left_u_tangent_cone = point_box(left_u_tangent_box, Eigen::Vector3d(0, 0, 0));
            cone<T, dim> left_v_tangent_cone = point_box(left_v_tangent_box, Eigen::Vector3d(0, 0, 0));
            T temp_angle_diff = std::max(left_u_tangent_cone.m_angle * 2.0, left_v_tangent_cone.m_angle * 2.0);

            cone<T, dim> right_u_tangent_cone = point_box(right_u_tangent_box, Eigen::Vector3d(0, 0, 0));
            cone<T, dim> right_v_tangent_cone = point_box(right_v_tangent_box, Eigen::Vector3d(0, 0, 0));
            T temp_angle_diff2 = std::max(right_u_tangent_cone.m_angle * 2.0, right_v_tangent_cone.m_angle * 2.0);
            angle_diff = std::min(std::max(temp_angle_diff, temp_angle_diff2), angle_diff);

            T minLen = c_box.eval_minimal_distance(Eigen::Vector3d(0, 0, 0));
            T maxLen = c_box.eval_maximal_distance(Eigen::Vector3d(0, 0, 0));

            //Interval<T> c_dot = dot_product(c, c);
            Interval<T> c_dot(minLen, maxLen);

            cone<T, dim> c_cone = point_box(c_box, Eigen::Vector3d(0, 0, 0));
            angle_diff = std::min(c_cone.m_angle * 2.0, angle_diff);

            //正常不应该在此return false
            if (c_dot.m_interval[0] <= TDEFAULT_ERROR<T>::value * TDEFAULT_ERROR<T>::value)
                return ENUM_NURBS::NURBS_ERROR;

            if (divide<T, (unsigned)dim>(c, c_dot, c) == false)
                return ENUM_NURBS::NURBS_ERROR;

            if (false == divide(dot_product(left_normal_interval, vector_product(c, left_v_tangent_interval)), u_normal_dot, m_alpha))
                return ENUM_NURBS::NURBS_ERROR;
            if (false == divide(dot_product(left_normal_interval, vector_product(left_u_tangent_interval, c)), u_normal_dot, m_beta))
                return ENUM_NURBS::NURBS_ERROR;
            if (false == divide(dot_product(right_normal_interval, vector_product(c, right_v_tangent_interval)), v_normal_dot, m_u))
                return ENUM_NURBS::NURBS_ERROR;
            if (false == divide(dot_product(right_normal_interval, vector_product(right_u_tangent_interval, c)), v_normal_dot, m_v))
                return ENUM_NURBS::NURBS_ERROR;

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

            Box<T, 4> temp_priori_enclosure = initial_box.plus(box_interval, PRECISION<T>::value);
            if (true == m_product_box.is_contain_box(temp_priori_enclosure, 5.0 * PRECISION<T>::value))
            {
                temp_priori_enclosure.intersect(m_product_box, priori_enclosure);
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            return ENUM_NURBS::NURBS_ERROR;
        }

        bool eval_priori_enclosure(const Box<T, 4>& initial_box, const T& bigger_step_size, T& smalll_step_size, Box<T, 4>& priori_enclosure, T& angle_diff, int &type, bool transversal = true)
        {
            int loopCount = 8;
            angle_diff = 4.0;
            Box<T, 4> initial_box_t = initial_box;
            Box<T, 4> priori_enclosure_t = priori_enclosure;
            T angle_diff_t = angle_diff;
            //ENUM_NURBS code = eval_priori_enclosure_inner2(initial_box, initial_box, bigger_step_size, priori_enclosure_t, angle_diff_t, type, transversal);
            //code = eval_priori_enclosure_inner(initial_box, initial_box, bigger_step_size, priori_enclosure, angle_diff, type, transversal);
            if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner2(initial_box, initial_box, bigger_step_size, priori_enclosure, angle_diff, type, transversal))
            {
                angle_diff = 4.0;
                return false;
            }
            while (loopCount > 0)
            {
                loopCount -= 1;
                //if (loopCount == 4)
                //{
                //    Eigen::Vector<T, 4> &min = priori_enclosure.Min;
                //    Eigen::Vector<T, 4> &max = priori_enclosure.Max;

                //    for (int index = 0; index < 2; ++index)
                //    {
                //        T startIndex = 2 * index;
                //        if ((max[startIndex] - min[startIndex]) > 10 * (max[startIndex + 1] - min[startIndex + 1]))
                //        {
                //            T mid = (max[1 + startIndex] + min[1 + startIndex]) / 2.0;
                //            T len = (max[startIndex] - min[startIndex]) / 10;
                //            max[1 + startIndex] = mid + len;
                //            min[1 + startIndex] = mid - len;
                //        }
                //        else if ((max[1 + startIndex] - min[1 + startIndex]) > 10 * (max[startIndex] - min[startIndex]))
                //        {
                //            T mid = (max[startIndex] + min[startIndex]) / 2.0;
                //            T len = (max[1 + startIndex] - min[1 + startIndex]) / 10;
                //            max[startIndex] = mid + len;
                //            min[startIndex] = mid - len;
                //        }
                //    }
                //    priori_enclosure.intersect(m_product_box, priori_enclosure);
                //}
                Box<T, 4> next_priori_enclosure;
                Box<T, 4> next_priori_enclosure_t;
                //code = eval_priori_enclosure_inner2(initial_box, priori_enclosure_t, smalll_step_size, priori_enclosure_t, angle_diff_t, type, transversal);
                //if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, priori_enclosure, smalll_step_size, next_priori_enclosure, angle_diff, type, transversal))
                if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner2(initial_box, priori_enclosure, smalll_step_size, next_priori_enclosure, angle_diff, type, transversal))
                {
                    return false;
                }
                if (priori_enclosure.is_contain_box(next_priori_enclosure))
                {
                    //code = eval_priori_enclosure_inner2(initial_box, next_priori_enclosure_t, smalll_step_size, priori_enclosure_t, angle_diff_t, type, transversal);
                    /*if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, next_priori_enclosure, smalll_step_size, priori_enclosure, angle_diff, type, transversal))*/
                    if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, next_priori_enclosure, smalll_step_size, priori_enclosure, angle_diff, type, transversal))
                    {
                        return false;
                    }
                    T costheta = std::cos(m_angle_eps);
                    if (angle_diff < m_angle_eps * 1000)
                    {
                        //Box<T, 4> box_interval;
                        //box_interval.Min[0] = m_alpha.m_interval[0];
                        //box_interval.Min[1] = m_beta.m_interval[0];
                        //box_interval.Min[2] = m_u.m_interval[0];
                        //box_interval.Min[3] = m_v.m_interval[0];

                        //box_interval.Max[0] = m_alpha.m_interval[1];
                        //box_interval.Max[1] = m_beta.m_interval[1];
                        //box_interval.Max[2] = m_u.m_interval[1];
                        //box_interval.Max[3] = m_v.m_interval[1];
                        //box_interval.scale(step_size);

                        //Box<T, 4> temp_priori_enclosure;
                        //if (true == initial_box.plus(box_interval, PRECISION<T>::value).intersect(m_product_box, temp_priori_enclosure))
                        //{
                        //    priori_enclosure = temp_priori_enclosure;
                        //    return true;
                        //}
                        return true;
                    }
                }
                //else
                Box<T, 4> temp_box = next_priori_enclosure;
                //priori_enclosure = ;
                Box<T, 4> temp_box_t = next_priori_enclosure;
                //code = eval_priori_enclosure_inner2(initial_box, temp_box_t, smalll_step_size, next_priori_enclosure_t, angle_diff_t, type, transversal);
                /*if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, temp_box, smalll_step_size, next_priori_enclosure, angle_diff, type, transversal))*/
                if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner2(initial_box, temp_box, smalll_step_size, next_priori_enclosure, angle_diff, type, transversal))
                {
                    return false;
                }
                if (temp_box.is_contain_box(next_priori_enclosure))
                {
                    //code = eval_priori_enclosure_inner2(initial_box, temp_box_t, smalll_step_size, priori_enclosure_t, angle_diff_t, type, transversal);
                    //if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, next_priori_enclosure, smalll_step_size, priori_enclosure, angle_diff, type, transversal))
                    if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner2(initial_box, next_priori_enclosure, smalll_step_size, priori_enclosure, angle_diff, type, transversal))
                    {
                        return false;
                    }
                    if (angle_diff < m_angle_eps * 100)
                    {
                        //Box<T, 4> box_interval;
                        //box_interval.Min[0] = m_alpha.m_interval[0];
                        //box_interval.Min[1] = m_beta.m_interval[0];
                        //box_interval.Min[2] = m_u.m_interval[0];
                        //box_interval.Min[3] = m_v.m_interval[0];

                        //box_interval.Max[0] = m_alpha.m_interval[1];
                        //box_interval.Max[1] = m_beta.m_interval[1];
                        //box_interval.Max[2] = m_u.m_interval[1];
                        //box_interval.Max[3] = m_v.m_interval[1];
                        //box_interval.scale(step_size);

                        //Box<T, 4> temp_priori_enclosure;
                        //if (true == initial_box.plus(box_interval, PRECISION<T>::value).intersect(m_product_box, temp_priori_enclosure))
                        //{
                        //    priori_enclosure = temp_priori_enclosure;
                        //    return true;
                        //}
                        return true;
                    }
                }
                //
                priori_enclosure = next_priori_enclosure;
                smalll_step_size /= 1.2;
                if (std::abs(smalll_step_size) < TDEFAULT_ERROR<T>::value)
                {
                    return false;
                }
            }
            return false;
        }

        //迭代加细交点(周期性曲面待处理)
        ENUM_NURBS intersect_point_iteration(const Box<T, 4>& domian, Eigen::Vector<T, 4> current_param, Eigen::Vector<T, 4>& intersect_param) 
        {
            // bool left_surf_u_closed = left_surf->is_u_closed();
            // bool left_surf_v_closed = left_surf->is_v_closed();
            // bool right_surf_u_closed = right_surf->is_u_closed();
            // bool right_surf_v_closed = right_surf->is_v_closed();

            Eigen::Vector<T, dim> left_point, right_point;
            Eigen::Vector<T, dim> vec;
            T min_distance = 100000;
            intersect_param = current_param;

            //迭代次数需要修改
            for (int loop_index = 0; loop_index < 2 * SURFACE_ITERATE_DEEP; ++loop_index)
            {
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
                        if (std::abs(vec[index]) > PRECISION<T>::value)
                        {
                            flag = false;
                            break;
                        }
                    }
                    if (flag == true)
                        return ENUM_NURBS::NURBS_SUCCESS;
                }

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



            }

            return ENUM_NURBS::NURBS_ERROR;
        }

        //迭代加细交点(周期性曲面待处理)
        ENUM_NURBS intersect_singular_point_iteration(const Box<T, 4>& domian, Eigen::Vector<T, 4> current_param, Eigen::Vector<T, 4>& intersect_param)
        {
            intersect_param = current_param;
            T min_distance = 10000000;
            //迭代次数需要修改
            for (int loop_index = 0; loop_index < SURFACE_ITERATE_DEEP; ++loop_index)
            {
                Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_derivatives, right_derivatives;
                m_left_surface->template derivative_on_surface<2>(current_param[0], current_param[1], left_derivatives);
                m_right_surface->template derivative_on_surface<2>(current_param[2], current_param[3], right_derivatives);
                Eigen::Vector<T, dim + 1> vec;
                vec.template block<dim, 1>(0, 0) =  right_derivatives(0, 0) - left_derivatives(0, 0);
                

                Eigen::Matrix<T, dim + 1, 4> mat;
                mat.template block<dim, 1>(0, 0) = left_derivatives(1, 0) ;
                mat.template block<dim, 1>(0, 1) = left_derivatives(0, 1) ;
                mat.template block<dim, 1>(0, 2) = -right_derivatives(1, 0);
                mat.template block<dim, 1>(0, 3) = -right_derivatives(0, 1);

                Eigen::Vector<T, dim> dsdt = right_derivatives(1, 0).cross(right_derivatives(0, 1));
                Eigen::Vector<T, dim> dudv = left_derivatives(1, 0).cross(left_derivatives(0, 1));
                Eigen::Vector<T, dim> dudvdsdt = dudv.cross(dsdt);
                Eigen::Vector<T, dim> dtdss = right_derivatives(0, 1).cross(right_derivatives(2, 0));
                Eigen::Vector<T, dim> dsdst = right_derivatives(1, 0).cross(right_derivatives(1, 1));
                Eigen::Vector<T, dim> dtdst = right_derivatives(0, 1).cross(right_derivatives(1, 1));
                Eigen::Vector<T, dim> dsdtt = right_derivatives(1, 0).cross(right_derivatives(0, 2));

                Eigen::Vector<T, dim> dvduu = left_derivatives(0, 1).cross(left_derivatives(2, 0));
                Eigen::Vector<T, dim> duduv = left_derivatives(1, 0).cross(left_derivatives(1, 1));
                Eigen::Vector<T, dim> dvduv = left_derivatives(0, 1).cross(left_derivatives(1, 1));
                Eigen::Vector<T, dim> dudvv = left_derivatives(1, 0).cross(left_derivatives(0, 2));

                double coeff = 2.0 / (dudv.squaredNorm() * dsdt.squaredNorm());
                //double coeff = 1.0;
                mat(dim, 0) = coeff * (-dvduu + duduv).cross(dsdt).dot(dudvdsdt);
                mat(dim, 1) = coeff * (-dvduv + dudvv).cross(dsdt).dot(dudvdsdt);
                mat(dim, 2) = coeff * dudv.cross(-dtdss + dsdst).dot(dudvdsdt);
                mat(dim, 3) = coeff * dudv.cross(-dtdst + dsdtt).dot(dudvdsdt);
                
                vec[dim] = -coeff * dudvdsdt.squaredNorm();
                T distance = vec.squaredNorm();
                if (distance < min_distance)
                {
                    min_distance = distance;
                    intersect_param = current_param;

                    bool flag = true;
                    for (int index = 0; index <= dim; ++index)
                    {
                        if (std::abs(vec[index]) > PRECISION<T>::value)
                        {
                            flag = false;
                            break;
                        }
                    }
                    if (flag == true)
                        return ENUM_NURBS::NURBS_SUCCESS;
                }

                Eigen::JacobiSVD<Eigen::Matrix<T, dim + 1, 4>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
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

                current_param = next_param;

            }

            return ENUM_NURBS::NURBS_ERROR;
        }

        template<bool direction>
        ENUM_NURBS trace_point(Box<T, 4> &initial_box, const Eigen::Vector<T, 4> &param ,surf_surf_intersect_point<T, dim>& next_point, Box<T, 4>& next_initial_box, bool &stop, bool transversal = true)
        {
            stop = false;
            //TODO:应在调用函数外部确保
            if (false == initial_box.intersect(m_product_box, initial_box))
                return ENUM_NURBS::NURBS_ERROR;
            int count1 = 0, count2 = 0;
            Eigen::Vector<T, 4> next_init_param;
            Box<T, 4> priori_enclosure;
            T param_space_eps = 2.0 * m_param_space_eps;
            Eigen::Vector<T, 4> intersect_param;
            priori_enclosure = initial_box;

           /* Eigen::Matrix<T, 4, 2> params_ders;
            params_ders.setConstant(0.0);
            Eigen::Matrix<T, dim, 3> space_ders;
            space_ders.setConstant(0.0);
            T x1, x2;
            eval_preiamge_ders<surface_type, typename surface_type::Type, surface_type::dimension>(m_left_surface, m_right_surface, param[0], param[1], param[2], param[3], params_ders, space_ders, x1, x2);*/
            std::vector<Eigen::Matrix<T, 4, 2>> param_ders;
            std::vector<Eigen::Matrix<T, dim, 2>> space_ders;
            eval_preiamge_and_space_ders<2>(param, param_ders, space_ders, transversal == true ? 0 : 3);

            // TODO : 处理branch
            T curvature = space_ders[0].col(1).norm();
            //T curvature = 0;
            T arc_length = 2 * m_angle_eps;
            if (curvature < TDEFAULT_ERROR<T>::value)
            {
                arc_length = TINFINITE<T>::value;
            }
            else
            {
                arc_length /= curvature;
            }
            T param_space_length = arc_length;

            // TODO : 处理step太大超过参数域
            //T abs_step = 0.5 * param_ders[0].col(0).norm() / (2.0 * param_ders[0].col(1).norm());
            //T bigger_step_size = direction == false ? abs_step : -abs_step;
            T bigger_step_size = direction == false ? param_space_length  : -param_space_length;
            T small_step = bigger_step_size;
            int type;
            T angle_diff;
            while (false == eval_priori_enclosure(initial_box, bigger_step_size, small_step, priori_enclosure, angle_diff, type, transversal))
            {
                if (std::abs(bigger_step_size) < 1e-4)
                {
                    stop = true;

                    return ENUM_NURBS::NURBS_SUCCESS;

                }
                if (angle_diff < 0.1)
                {
                    //if (angle_diff < 0.01)
                    //{
                    //    return ENUM_NURBS::NURBS_ERROR;
                    //}
                    // 计算奇异点
                    Eigen::Vector<T, 4> pre_param;
                    if (estimate_next_param(param, priori_enclosure, small_step, next_init_param, type) != ENUM_NURBS::NURBS_SUCCESS)
                    {
                        return ENUM_NURBS::NURBS_ERROR;
                    }

                    // try find singular point
                    if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration(m_product_box, next_init_param, pre_param))
                    {
                        if (ENUM_NURBS::NURBS_SUCCESS == intersect_singular_point_iteration(m_product_box, pre_param, intersect_param))
                        {
                            if (priori_enclosure.is_contain_point(intersect_param) == true)
                            {
                                // check is singular point
                                    // TODO: 计算出奇异点，从此点出沿着切向发先出发一个微小的距离到达p，从p点开始追踪(存在唯一的追踪)，过程中可能链接到存在的曲线。追踪完之后，继续追踪因为奇异点导致停止的追踪过程
                                std::vector<Eigen::Matrix<T, 4, 1>> param_ders_temp;
                                std::vector<Eigen::Matrix<T, dim, 1>> space_ders_temp;
                                eval_preiamge_and_space_ders<1>(intersect_param, param_ders_temp, space_ders_temp, 3);
                                if (space_ders_temp.size() == 2)
                                {
                                    stop = true;
                                    m_left_surface->point_on_surface(intersect_param[0], intersect_param[1], next_point.m_point);
                                    next_point.m_uv = intersect_param;

                                    next_initial_box.Min = intersect_param - m_box_size;
                                    next_initial_box.Max = intersect_param + m_box_size;

                                    return ENUM_NURBS::NURBS_SUCCESS;
                                    //break;
                                }
                            }      
                        }
                    }  
                }
                bigger_step_size /= 1.5;
                small_step = bigger_step_size;
            }

            //T temp_bigger = bigger_step_size / 1.2;
            //T temp_small = small_step * 1.2;
            //while (true)
            //{
            //    Box<T, 4> temp_priori_enclosure;
            //    if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, initial_box, temp_bigger, temp_priori_enclosure, angle_diff, type, transversal))
            //    {
            //        break;
            //    }
            //    Box<T, 4> next_priori_enclosure;
            //    if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, temp_priori_enclosure, temp_small, next_priori_enclosure, angle_diff, type, transversal))
            //    {
            //        break;
            //    }
            //    if (temp_priori_enclosure.is_contain_box(next_priori_enclosure))
            //    {
            //        priori_enclosure = temp_priori_enclosure;
            //        small_step = temp_small;
            //        bigger_step_size = temp_bigger;
            //        temp_small *= 1.2;
            //        temp_bigger /= 1.2;
            //    }
            //    else
            //    {
            //        temp_bigger /= 1.2;
            //        if (std::abs(temp_bigger) < std::abs(temp_small))
            //        {
            //            break;
            //        }
            //    }
            //}
            

            do
            {
                if (ENUM_NURBS::NURBS_SUCCESS != estimate_next_param(param, priori_enclosure, small_step, next_init_param, type))
                {
                    return ENUM_NURBS::NURBS_ERROR;
                }
                if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration(m_product_box, next_init_param, intersect_param))
                {
                    if (transversal == false)
                    {
                        next_init_param = intersect_param;
                        if (ENUM_NURBS::NURBS_SUCCESS == intersect_singular_point_iteration(m_product_box, next_init_param, intersect_param))
                        {
                            if (priori_enclosure.is_contain_point(intersect_param) == true)
                            {
                                break;
                            }
                        }
                    }
                    else
                    {
                        if (priori_enclosure.is_contain_point(intersect_param) == true)
                        {
                            break;
                        }
                    }
                }
                small_step /= 1.2;

            } while (true);

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
            Eigen::Vector<T, 4> current_param = initial_box.get_middle_point();
            Box<T, 4> next_initial_box;

            Eigen::Vector<T, 4> initial_param;
            surf_surf_intersect_point<T, dim>* current_intersection = new surf_surf_intersect_point<T, dim>();
            Eigen::Vector<T, dim> initial_point;

            if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration(Box<T, 4>(Eigen::Vector4<T>(0, 0, 0, 0), Eigen::Vector4<T>(1, 1, 1, 1)), current_param, initial_param))
            {
                current_param = initial_param;
                //if (ENUM_NURBS::NURBS_SUCCESS != intersect_singular_point_iteration(Box<T, 4>(Eigen::Vector4<T>(0, 0, 0, 0), Eigen::Vector4<T>(1, 1, 1, 1)), current_param, initial_param))
                //{
                //    initial_param = current_param;
                //}
            }

            //m_current_param = initial_param;
            m_left_surface->point_on_surface(initial_param[0], initial_param[1], initial_point);
            current_intersection->m_point = initial_point;
            current_intersection->m_uv = initial_param;
            intersection->m_intersect_points.push_back(current_intersection);

            int index = 1;
            Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> left_ders;
            m_left_surface->template derivative_on_surface<1>(initial_param[0], initial_param[1], left_ders);
            Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_ders;
            m_right_surface->template derivative_on_surface<1>(initial_param[2], initial_param[3], right_ders);

            Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
            left_normal.normalize();

            Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
            right_normal.normalize();
            bool transversal = true;
            double cosTheta = left_normal.dot(right_normal);

            if (cosTheta > 1.0 - KNOTS_EPS<T>::value)
            {
                transversal = false;
            }

            for (; index < MAXINTERSETORPOINTNUMBER; ++index)
            {
                std::cout << index << std::endl;
                surf_surf_intersect_point<T, dim>* next_point = new surf_surf_intersect_point<T, dim>();
                bool stop = false;
                if (trace_point<direction>(current_box, current_intersection->m_uv, *next_point, next_initial_box, stop, transversal) != ENUM_NURBS::NURBS_SUCCESS)
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

                if (stop == true)
                {
                    intersection->m_intersect_points_count.push_back(index + 1);
                    return ENUM_NURBS::NURBS_SUCCESS;
                }

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

            
            std::vector<Eigen::Vector<T, 3>> int_params2 = bezier_curve_int_bezier_surface(curve2, *m_left_surface);
            std::vector<Eigen::Vector<T, 3>> int_params1 = bezier_curve_int_bezier_surface(curve1, *m_left_surface);
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
            m_current_param = Eigen::Vector<T, 4>(0.5, 1, 0.5, 1);
            int index = 0;
            for (auto &box : int_params)
            {
                //if (index != 3)
                //{
                //    ++index;
                //    continue;
                //}
                trace_curve_by_point<true>(box, *intersection);
                trace_curve_by_point<false>(box, *intersection);
            }
            return ENUM_NURBS::NURBS_SUCCESS;

        }


    };

}