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
    
    

    enum POINT_TYPE
    {
        UNKNOWN = -1,
        BOUNDARY = 0,
        SINGULAR = 1,
        LOOP = 2
    };

    //intersect points(以后有空再改吧)
    template<typename T, int dim>
    struct surf_surf_intersect_point
    {
        //surf_surf_intersect_point<T, dim> *m_next;
        Eigen::Vector<T, dim> m_point;
        Eigen::Vector4<T> m_uv;
        Box<T, 4> m_priori_enclosure;
    };

    template<typename T, int dim>
    struct surfs_int_points_chat
    {
        std::vector<surf_surf_intersect_point<T, dim>> m_inter_points;
        bool m_is_transversal{true};
        bool m_is_positive_direction{ true };
        POINT_TYPE start_point_type{ POINT_TYPE::UNKNOWN };
        POINT_TYPE end_point_type{ POINT_TYPE::UNKNOWN };
    };


    //面面相交的结果
    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    struct surf_surf_int
    {        
        surface_type* m_left_surf;
        surface_type* m_right_surf;
        std::vector<std::vector<surf_surf_intersect_point<T, dim>>> m_intersect_points;
        std::vector<surfs_int_points_chat<T, dim>> m_int_chats;

        //每个交点链表的交点的个数
        std::vector<int> m_intersect_points_count;
        ~surf_surf_int()
        {
            //if (m_left_surf)
            //    delete m_left_surf;
            //if (m_right_surf)
            //    delete m_right_surf;
        }
    };


    // TODO: 整理和完善代码
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

        T m_angle_eps;
        Box<T, 4> m_product_box;

        //std::vector<std::vector<surf_surf_intersect_point<T, dim>>> m_boundary_intpoints;
        std::vector<surfs_int_points_chat<T, dim>> m_boundary_chats;
        
        //'0'表示负向追踪，其余反向追踪; std::vector<bool>真sb
        std::vector<char> m_boundary_point_trace_dir;
        std::vector<char> m_singular_point_trace_dir;
        //'1'表示transversal
        //std::vector<char> m_boundary_transversal;

        //奇异点必然为起始点或者终止点
        //奇异点的意思是在此点的切平面上存在至少两个方向v1, v2, 使得两个相交的平面在这两个方向上任意阶偏导数相等(弧长参数下), 奇异点必然相切
        //确切的说这里存的是奇异点的前一个点
        //std::vector<std::vector<surf_surf_intersect_point<T, dim>>> m_singular_points;
        std::vector<surfs_int_points_chat<T, dim>> m_singular_chats;
        std::vector<Box<T, 4>> m_singular_boxes;

        //闭曲线的追踪点
        std::vector<surfs_int_points_chat<T, dim>> m_loop_chats;

        ////闭曲线的在参数域的切向量角度的积分，简单闭曲线的这个值应该为2pi或者-2pi
        //std::vector<T> m_loop_angles;

    public:
        surf_surf_int<surface_type> m_result;

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
            //m_result = new surf_surf_int<surface_type>();
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

            m_angle_eps = angle_eps;
            m_product_box = m_left_surface->get_interval().product_box(m_right_surface->get_interval());

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
        //validated solution of initial value可以计算交线在此部分的包围盒,此包围盒应该会有其他的一些用处
        template<unsigned order>
        ENUM_NURBS eval_preiamge_and_space_ders(const Box<T, 4> &param_box, std::array<Box<T, 4>, order> &param_ders, 
            std::array<Box<T, dim>, order> &space_ders, bool& may_arrived_singular)
        {
            //目前仅仅支持三维
            static_assert(3 == dim, "dim != 3 is not supported");
            static_assert(order == 1 || order == 2, "order != 1 && order != 2");
            may_arrived_singular = false;
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
            {
                may_arrived_singular = true;
                return error_code;
            }
            //space_ders.resize(1);
            space_ders[0] = c;
            //param_ders.resize(1);

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
            param_ders[0].set_index_interval(0, alpha_d[0]);
            param_ders[0].set_index_interval(1, alpha_d[1]);
            param_ders[0].set_index_interval(2, beta_d[0]);
            param_ders[0].set_index_interval(3, beta_d[1]);
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
                space_ders[1] = c_dd;

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

                param_ders[1].set_index_interval(0, coeff1);
                param_ders[1].set_index_interval(1, coeff2);

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

                param_ders[1].set_index_interval(2, coeff1);
                param_ders[1].set_index_interval(3, coeff2);
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        //相切的情形的一阶微分的box
        //validated solution of initial value可以计算交线在此部分的包围盒,此包围盒应该会有其他的一些用处
        //目前只能使用两个surface的box内部normal的变化作为一个步长选择的标准
        template<unsigned order>
        ENUM_NURBS eval_preiamge_and_space_ders_tangent(const Box<T, 4>& param_box, std::array<Box<T, 4>, order>& param_ders,
            std::array<Box<T, dim>, order>& space_ders, int& type)
        {
            //目前仅仅支持三维
            static_assert(3 == dim, "dim != 3 is not supported");
            static_assert(order == 1, "order != 1");
            Box<T, dim> left_normal_box, right_normal_box;
            Box<T, dim> left_u_tangent_box, left_v_tangent_box;
            Box<T, dim> right_u_tangent_box, right_v_tangent_box;

            Box<T, 2> left_param_box(param_box.Min.template block<2, 1>(0, 0), param_box.Max.template block<2, 1>(0, 0));
            Box<T, 2> right_param_box(param_box.Min.template block<2, 1>(2, 0), param_box.Max.template block<2, 1>(2, 0));
            eval_point_interval(m_left_u_tangent_surface, left_param_box, left_u_tangent_box);
            eval_point_interval(m_left_v_tangent_surface, left_param_box, left_v_tangent_box);
            eval_point_interval(m_right_u_tangent_surface, right_param_box, right_u_tangent_box);
            eval_point_interval(m_right_v_tangent_surface, right_param_box, right_v_tangent_box);

            eval_point_interval(m_left_normal_surface, left_param_box, left_normal_box);
            eval_point_interval(m_right_normal_surface, right_param_box, right_normal_box);

            Box<T, 1> left_normal_len2 = interval_algorithm::dot(left_normal_box, left_normal_box);
            Box<T, 1> right_normal_len2 = interval_algorithm::dot(right_normal_box, right_normal_box);

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

            Box<T, 1> normal_dot = interval_algorithm::dot(left_normal_box, right_normal_box);
            Eigen::Vector<T, 1> zero_num;
            zero_num[0] = 0;
            if (normal_dot.is_contain_point(zero_num, TDEFAULT_ERROR<T>::value) == true)
            {
                return ENUM_NURBS::NURBS_ERROR;
            }
            //int flag = normal_dot.Min[0] > 0 ? 1 : -1;
            Box<T, dim> conormal;
            ENUM_NURBS code = interval_algorithm::normalized(right_normal_box, conormal);
            if (ENUM_NURBS::NURBS_SUCCESS != code)
            {
                return code;
            }
            
            Eigen::Vector<T, dim> origin;
            origin.setConstant(0.0);
            T right_normal_min_length = right_normal_box.eval_minimal_distance(origin);
            T right_normal_max_length = right_normal_box.eval_minimal_distance(origin);
            Box<T, 1> denominator;
            denominator.Min[0] = right_normal_min_length;
            denominator.Max[0] = right_normal_max_length;
            Box<T, 1> denominator_r, unitInterval;
            unitInterval.Min[0] = 1.0;
            unitInterval.Max[0] = 1.0;
            code = interval_algorithm::divide(unitInterval, denominator, denominator_r);
            if (ENUM_NURBS::NURBS_SUCCESS != code)
            {
                return code;
            }

            Box<T, 1> a11 = interval_algorithm::dot(interval_algorithm::cross(left_u_tangent_box, right_v_tangent_box), conormal) * denominator_r;
            Box<T, 1> a12 = interval_algorithm::dot(interval_algorithm::cross(left_v_tangent_box, right_v_tangent_box), conormal) * denominator_r;
            Box<T, 1> a21 = interval_algorithm::dot(interval_algorithm::cross(right_u_tangent_box, left_u_tangent_box), conormal) * denominator_r;
            Box<T, 1> a22 = interval_algorithm::dot(interval_algorithm::cross(right_u_tangent_box, left_v_tangent_box), conormal) * denominator_r;
            Box<T, dim> Xuu_box, Xuv_box, Xvv_box, Yuu_box, Yuv_box, Yvv_box;
            eval_point_interval(m_left_dxx[0], left_param_box, Xuu_box);
            eval_point_interval(m_left_dxx[1], left_param_box, Xuv_box);
            eval_point_interval(m_left_dxx[2], left_param_box, Xvv_box);

            eval_point_interval(m_right_dxx[0], right_param_box, Yuu_box);
            eval_point_interval(m_right_dxx[1], right_param_box, Yuv_box);
            eval_point_interval(m_right_dxx[2], right_param_box, Yvv_box);
            
            Box<T, 1> left_L = interval_algorithm::dot(Xuu_box, conormal);
            Box<T, 1> left_M = interval_algorithm::dot(Xuv_box, conormal);
            Box<T, 1> left_N = interval_algorithm::dot(Xvv_box, conormal);

            Box<T, 1> right_L = interval_algorithm::dot(Yuu_box, conormal);
            Box<T, 1> right_M = interval_algorithm::dot(Yuv_box, conormal);
            Box<T, 1> right_N = interval_algorithm::dot(Yvv_box, conormal);

            Box<T, 1> b11 = a11 * a11;
            b11.Min[0] = std::max(0.0, b11.Min[0]);
            b11 = b11 * right_L;
            b11 = b11 + (a11 * a21 * right_M) * 2.0;
            Box<T, 1> temp = a21 * a21;
            temp.Min[0] = std::max(0.0, temp.Min[0]);
            b11 = b11 + temp * right_N;
            b11 = b11 - left_L;

            Box<T, 1> b12 = a11 * a12 * right_L;
            b12 = b12 + (a11 * a22 + a21 * a12) * right_M;
            b12 = b12 + a21 * a22 * right_N;
            b12 = b12 - left_M;

            Box<T, 1> b22 = a12 * a12;
            b22.Min[0] = std::max(0.0, b22.Min[0]);
            b22 = b22 * right_L;
            b22 = b22 + (a12 * a22 * right_M) * 2.0;
            temp = a22 * a22;
            temp.Min[0] = std::max(0.0, temp.Min[0]);
            b22 = b22 + temp * right_N;
            b22 = b22 - left_N;
            
            bool flag1 = b11.is_contain_point(zero_num, TDEFAULT_ERROR<T>::value);
            //bool flag2 = b12.is_contain_point(zero_num, TDEFAULT_ERROR<T>::value);
            bool flag2 = b22.is_contain_point(zero_num, TDEFAULT_ERROR<T>::value);
            if (flag1 && flag2)
            {
                return ENUM_NURBS::NURBS_ERROR;
            }
            T b11_extream = std::min(std::abs(b11.Min[0]), std::abs(b11.Max[0]));
            T b22_extream = std::min(std::abs(b22.Min[0]), std::abs(b22.Max[0]));
            
            if (flag2 == true || (b11_extream > b22_extream && flag1 == false))
            {
                Box<T, 1> nu;
                interval_algorithm::divide(b12, b11, nu);
                Box<T, dim> tangent = left_u_tangent_box * nu + left_v_tangent_box;
                code = interval_algorithm::normalized(tangent, space_ders[0]);
                if (ENUM_NURBS::NURBS_SUCCESS != code)
                {
                    return code;
                }
                type = 1;
            }
            else
            {
                Box<T, 1> mu;
                interval_algorithm::divide(b12, b22, mu);
                Box<T, dim> tangent = left_v_tangent_box * mu + left_u_tangent_box;
                code = interval_algorithm::normalized(tangent, space_ders[0]);
                if (ENUM_NURBS::NURBS_SUCCESS != code)
                {
                    return code;
                }
                type = 2;
            }
            std::vector<Box<T, 1>> alpha_d(2);
            std::vector<Box<T, 1>> beta_d(2);
            if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(space_ders[0], left_v_tangent_box)), left_normal_len2, alpha_d[0]))
                return ENUM_NURBS::NURBS_ERROR;

            if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(left_u_tangent_box, space_ders[0])), left_normal_len2, alpha_d[1]))
                return ENUM_NURBS::NURBS_ERROR;

            if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(space_ders[0], right_v_tangent_box)), right_normal_len2, beta_d[0]))
                return ENUM_NURBS::NURBS_ERROR;

            if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(right_u_tangent_box, space_ders[0])), right_normal_len2, beta_d[1]))
                return ENUM_NURBS::NURBS_ERROR;
            param_ders[0].set_index_interval(0, alpha_d[0]);
            param_ders[0].set_index_interval(1, alpha_d[1]);
            param_ders[0].set_index_interval(2, beta_d[0]);
            param_ders[0].set_index_interval(3, beta_d[1]);

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
                //int flag = right_normal.dot(left_normal) > 0 ? 1 : -1;
                T denominator = 1.0 / right_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal);
                T a11 = left_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal) * denominator;
                T a12 = left_ders(0, 1).cross(right_ders(0, 1)).dot(right_normal) * denominator;
                T a21 = right_ders(1, 0).cross(left_ders(1, 0)).dot(right_normal) * denominator;
                T a22 = right_ders(1, 0).cross(left_ders(0, 1)).dot(right_normal) * denominator;
                T left_L = left_ders(2, 0).dot(right_normal);
                T left_M = left_ders(1, 1).dot(right_normal);
                T left_N = left_ders(0, 2).dot(right_normal);

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
                    // 下面delta的值需要修改
                    if (delta < -TDEFAULT_ERROR<T>::value)
                    {
                        T try_value = -b12 / b11;
                        T value = b11 * (try_value * try_value) + 2 * b12 * try_value + b22;
                        if (std::abs(value) < TDEFAULT_ERROR<T>::value)
                        {
                            delta = 0.0;
                        }
                        else
                        {
                            return ENUM_NURBS::NURBS_ISOLATED_TANGENTIAL_POINT;
                        }
                    }
                    if (delta > TDEFAULT_ERROR<T>::value)
                    {
                        is_branch = true;
                    }
                    if (std::abs(b11) < TDEFAULT_ERROR<T>::value && std::abs(b12) < TDEFAULT_ERROR<T>::value
                        && std::abs(b22) < TDEFAULT_ERROR<T>::value)
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
            int loop_count = 8;
            while (loop_count-- > 0)
            {
                //domain = Box<T, 4>(Eigen::Vector<T, 4>(0, 0, 0, 0), Eigen::Vector<T, 4>(1, 1, 1, 1));
                //Eigen::Vector<T, dim> tangent;
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


        bool may_contain_param(const Eigen::Vector2<T>& current_param, const Eigen::Vector2<T>& next_param, const Eigen::Vector2<T>& test_param)
        {
            Eigen::Vector2<T> param_vec = next_param - current_param;
            T vec_len = param_vec.norm();
            param_vec /= vec_len;
            //param_vec.normalize();
            Eigen::Vector2<T> test_vec = test_param - current_param;
            T project_len = test_vec.dot(param_vec);
            if (project_len > -PRECISION<T>::value && project_len < vec_len + PRECISION<T>::value)
            {
                return true;
            }
            return false;
        }


        bool is_point_on_arc(Eigen::Vector4<T>& current_param, Eigen::Vector2<T>& next_param, const Eigen::Vector2<T>& test_param, const Box<T, 4>& priori_enclosure, int type, bool is_positive)
        {
            Eigen::Vector<T, dim> test_point;
            m_right_surface->point_on_surface(test_param[0], test_param[1], test_point);
            if (may_contain_param(current_param.template block<2, 1>(2, 0), next_param, test_param))
            {
                for (int iter_index = 0; iter_index < 2 * SURFACE_ITERATE_DEEP; ++iter_index)
                {
                    Eigen::Vector<T, dim> start_point, end_point;

                    m_right_surface->point_on_surface(current_param[2], current_param[3], start_point);
                    m_right_surface->point_on_surface(next_param[0], next_param[1], end_point);

                    Eigen::Vector<T, dim> chord_vec = end_point - start_point;
                    T chord_length = chord_vec.norm();
                    chord_vec /= chord_length;
                    Eigen::Vector<T, dim> test_vec = test_point - start_point;
                    T test_length = test_vec.dot(chord_vec);
                    if (is_positive == false)
                    {
                        test_length = -test_length;
                    }
                    Eigen::Vector4<T> next_init_param, intersect_param;
                    estimate_next_param(current_param, priori_enclosure, test_length, next_init_param, type);
                    intersect_point_iteration(m_product_box, next_init_param, intersect_param);
                    Eigen::Vector<T, dim> intersec_point;
                    
                    m_right_surface->point_on_surface(intersect_param[2], intersect_param[3], intersec_point);
                    double dis = (intersec_point - test_point).norm();
                    if (priori_enclosure.is_contain_point(intersect_param) && (intersec_point - test_point).norm() < 1e-8)
                    {
                        return true;
                    }
                    Eigen::Vector2<T> mid_param = intersect_param.template block<2, 1>(2, 0);
                  
                    bool flag1 = may_contain_param(current_param.template block<2, 1>(2, 0), mid_param, test_param);
                    bool flag2 = may_contain_param(mid_param, next_param, test_param);

                    if (flag1 == true && flag2 == true)
                    {
                        std::cout << "error" << std::endl;
                    }
                    if (flag1 == true)
                    {
                        next_param = mid_param;
                    }
                    else if (flag2 == true)
                    {
                        current_param.template block<2, 1>(2, 0) = mid_param;
                        m_left_surface->find_point_on_surface(intersec_point, current_param[0], current_param[1]);
                    }
                    else
                    {
                        break;
                    }
                }
            }
            return false;
        }

        bool is_closed(surfs_int_points_chat<T, dim>& int_chat, const surf_surf_intersect_point<T, dim>& int_point, int type)
        {
            // TODO: 重写
            if (int_chat.m_inter_points.size() > 3 && int_point.m_priori_enclosure.is_contain_point(int_chat.m_inter_points[0].m_uv))
            {
                Eigen::Vector2<T> test_param = int_chat.m_inter_points[0].m_uv.template block<2, 1>(2, 0);
                Eigen::Vector2<T> next_param= int_point.m_uv.template block<2, 1>(2, 0);
                Eigen::Vector4<T> current_param = int_chat.m_inter_points.back().m_uv;
                return is_point_on_arc(current_param, next_param, test_param, int_point.m_priori_enclosure, type, int_chat.m_is_positive_direction);
            }
            return false;
        }

        bool is_point_on_intcurve(const Eigen::Vector4<T>& param, bool is_transversal, int type)
        {
            //Eigen::Vector2<T> test_param = param.template block<2, 1>(0, 0);
            Eigen::Vector2<T> test_param = param.template block<2, 1>(2, 0);
            //Eigen::Vector<T, dim> test_point;
            ////m_left_surface->point_on_surface(test_param[0], test_param[1], test_point);
            //m_right_surface->point_on_surface(test_param[0], test_param[1], test_point);
            for (surfs_int_points_chat<T, dim>& int_chat : m_result.m_int_chats)
            {
                if (is_transversal != int_chat.m_is_transversal)
                {
                    continue;
                }
                int int_points_count = int_chat.m_inter_points.size();
                for (int index = 0; index < int_points_count - 1; ++index)
                {
                    const Box<T, 4>& priori_enclosure = int_chat.m_inter_points[index].m_priori_enclosure;
                    if (priori_enclosure.is_contain_point(param))
                    {
                      /*  Eigen::Vector2<T> current_param = int_chat.m_inter_points[index].m_uv.template block<2, 1>(0, 0);*/
                        
                        
                        Eigen::Vector4<T> current_param = int_chat.m_inter_points[index].m_uv;
                        //Eigen::Vector2<T> next_param = int_chat.m_inter_points[index + 1].m_uv.template block<2, 1>(0, 0);
                        Eigen::Vector2<T> next_param = int_chat.m_inter_points[index + 1].m_uv.template block<2, 1>(2, 0);
                        
                        //Eigen::Vector4<T> next_param = int_chat.m_inter_points[index + 1].m_uv;
                        if (is_point_on_arc(current_param, next_param, test_param, priori_enclosure, type, int_chat.m_is_positive_direction) == true)
                        {
                            return true;
                        }


                        //if (may_contain_param(current_param.template block<2, 1>(0, 0), next_param, test_param))
                        //if (may_contain_param(current_param.template block<2, 1>(2, 0), next_param, test_param))
                        //{
                        //    for (int iter_index = 0; iter_index < 2 * SURFACE_ITERATE_DEEP; ++iter_index)
                        //    {
                        //        Eigen::Vector<T, dim> start_point, end_point;
                        //        //
                        //        //m_left_surface->point_on_surface(current_param[0], current_param[1], start_point);
                        //        //m_left_surface->point_on_surface(next_param[0], next_param[1], end_point);

                        //        m_right_surface->point_on_surface(current_param[2], current_param[3], start_point);
                        //        m_right_surface->point_on_surface(next_param[0], next_param[1], end_point);

                        //        Eigen::Vector<T, dim> chord_vec = end_point - start_point;
                        //        T chord_length = chord_vec.norm();
                        //        chord_vec /= chord_length;
                        //        Eigen::Vector<T, dim> test_vec = test_point - start_point;
                        //        T test_length = test_vec.dot(chord_vec);
                        //        if (int_chat.m_is_positive_direction == false)
                        //        {
                        //            test_length = -test_length;
                        //        }
                        //        Eigen::Vector4<T> next_init_param, intersect_param;
                        //        estimate_next_param(current_param, priori_enclosure, test_length, next_init_param, type);
                        //        intersect_point_iteration(m_product_box, next_init_param, intersect_param);
                        //        Eigen::Vector<T, dim> intersec_point;
                        //        //m_left_surface->point_on_surface(intersect_param[0], intersect_param[1], intersec_point);
                        //        m_right_surface->point_on_surface(intersect_param[2], intersect_param[3], intersec_point);
                        //        double dis = (intersec_point - test_point).norm();
                        //        if (priori_enclosure.is_contain_point(intersect_param) && (intersec_point - test_point).norm() < 1e-8)
                        //        {
                        //            return true;
                        //        }
                        //        //Eigen::Vector2<T> mid_param = intersect_param.template block<2, 1>(0, 0);
                        //        Eigen::Vector2<T> mid_param = intersect_param.template block<2, 1>(2, 0);
                        //        //Eigen::Vector4<T> mid_param = intersect_param;
                        //        //bool flag1 = may_contain_param(current_param.template block<2, 1>(0, 0), mid_param, test_param);
                        //        bool flag1 = may_contain_param(current_param.template block<2, 1>(2, 0), mid_param, test_param);
                        //        bool flag2 = may_contain_param(mid_param, next_param, test_param);

                        //        if (flag1 == true && flag2 == true)
                        //        {
                        //            std::cout << "error" << std::endl;
                        //        }
                        //        if (flag1 == true)
                        //        {
                        //            //current_param = mid_param;
                        //            next_param = mid_param;
                        //        }
                        //        else if (flag2 == true)
                        //        {
                        //            //current_param.template block<2, 1>(0, 0) = mid_param;
                        //            //m_right_surface->find_point_on_surface(intersec_point, current_param[2], current_param[3]);
                        //            
                        //            current_param.template block<2, 1>(2, 0) = mid_param;
                        //            m_left_surface->find_point_on_surface(intersec_point, current_param[0], current_param[1]);
                        //            
                        //            //next_param = mid_param;
                        //        }
                        //        else
                        //        {
                        //            break;
                        //        }
                        //    }
                        //}
                    
                    }
                }
            }
            return false;
        }
        
        
        ENUM_NURBS eval_priori_enclosure_inner(const Box<T, 4>& initial_box, const Box<T, 4>& param_box, T step_size, Box<T, 4>& priori_enclosure, int& type, bool& may_arrived_singular, bool transversal = true)
        {
            may_arrived_singular = false;
            std::array<Box<T, 4>, 1> param_ders;
            std::array<Box<T, dim>, 1> space_ders;
            ENUM_NURBS error_code;
            if (transversal == true)
            {
                error_code = eval_preiamge_and_space_ders<1>(param_box, param_ders, space_ders, may_arrived_singular);
                type = 0;
            }
            else
            {
                error_code = eval_preiamge_and_space_ders_tangent<1>(param_box, param_ders, space_ders, type);
            }
            if (ENUM_NURBS::NURBS_SUCCESS != error_code)
            {
                return error_code;
            }

            Interval<T> step_box(0, step_size);

            Box<T, 4> box_interval = param_ders[0] * step_box;
            Box<T, 4> temp_priori_enclosure = initial_box + box_interval;

            if (true == m_product_box.is_contain_box(temp_priori_enclosure, 5.0 * PRECISION<T>::value))
            {
                temp_priori_enclosure.intersect(m_product_box, priori_enclosure);
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            return ENUM_NURBS::NURBS_ERROR;


        }

        bool eval_priori_enclosure(const Box<T, 4>& initial_box, const T& bigger_step_size, T& smalll_step_size, Box<T, 4>& priori_enclosure, T& angle_diff, int &type, bool& may_arrived_singular, bool transversal = true)
        {
            may_arrived_singular = false;
            int loopCount = 8;
            angle_diff = 4.0;
            Box<T, 4> initial_box_t = initial_box;
            T angle_diff_t = angle_diff;
            if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, initial_box, bigger_step_size, priori_enclosure, type, may_arrived_singular, transversal))
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
                if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, priori_enclosure, smalll_step_size, next_priori_enclosure, type, may_arrived_singular, transversal))
                {
                    Box<T, 2> left_param_box(priori_enclosure.Min.template block<2, 1>(0, 0), priori_enclosure.Max.template block<2, 1>(0, 0));
                    Box<T, 2> right_param_box(priori_enclosure.Min.template block<2, 1>(2, 0), priori_enclosure.Max.template block<2, 1>(2, 0));
                    Box<T, dim> left_normal_box, right_normal_box;
                    eval_point_interval(m_left_normal_surface, left_param_box, left_normal_box);
                    eval_point_interval(m_right_normal_surface, right_param_box, right_normal_box);
                    Eigen::Vector<T, dim> origin;
                    origin.setConstant(0.0);
                    cone<T, dim> normal_cone = point_box(left_normal_box, origin);
                    angle_diff = normal_cone.m_angle;
                    normal_cone = point_box(right_normal_box, origin);
                    angle_diff = std::max(angle_diff, normal_cone.m_angle);

                    return false;
                }
                if (priori_enclosure.is_contain_box(next_priori_enclosure))
                {
                    if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, next_priori_enclosure, smalll_step_size, priori_enclosure, type, may_arrived_singular, transversal))
                    {
                        return false;
                    }
                    return true;
                }
                //else
                Box<T, 4> temp_box = next_priori_enclosure;
                if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, temp_box, smalll_step_size, next_priori_enclosure, type, may_arrived_singular, transversal))
                {
                    priori_enclosure = temp_box;
                    Box<T, 2> left_param_box(priori_enclosure.Min.template block<2, 1>(0, 0), priori_enclosure.Max.template block<2, 1>(0, 0));
                    Box<T, 2> right_param_box(priori_enclosure.Min.template block<2, 1>(2, 0), priori_enclosure.Max.template block<2, 1>(2, 0));
                    Box<T, dim> left_normal_box, right_normal_box;
                    eval_point_interval(m_left_normal_surface, left_param_box, left_normal_box);
                    eval_point_interval(m_right_normal_surface, right_param_box, right_normal_box);
                    Eigen::Vector<T, dim> origin;
                    origin.setConstant(0.0);
                    cone<T, dim> normal_cone = point_box(left_normal_box, origin);
                    angle_diff = normal_cone.m_angle;
                    normal_cone = point_box(right_normal_box, origin);
                    angle_diff = std::max(angle_diff, normal_cone.m_angle);
                    
                    return false;
                }
                if (temp_box.is_contain_box(next_priori_enclosure))
                {
                    if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, next_priori_enclosure, smalll_step_size, priori_enclosure, type, may_arrived_singular, transversal))
                    {
                        return false;
                    }
                    return true;
                }
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

        ENUM_NURBS trace_point(surfs_int_points_chat<T, dim>& current_intpoints, bool& is_stop, bool& is_arrived_singular, T box_len = TDEFAULT_ERROR<T>::value * 0.01)
        {
            //bool direction;
            bool transversal = current_intpoints.m_is_transversal;
            is_stop = false;
            is_arrived_singular = false;
            surf_surf_intersect_point<T, dim>& current_intpoint = current_intpoints.m_inter_points.back();
            Eigen::Vector<T, 4> param = current_intpoints.m_inter_points.back().m_uv;
            Box<T, 4> initial_box;
            create_box(param, box_len, initial_box);
            initial_box.intersect(m_product_box, initial_box);

            std::vector<Eigen::Matrix<T, 4, 2>> param_ders;
            std::vector<Eigen::Matrix<T, dim, 2>> space_ders;
            eval_preiamge_and_space_ders<2>(param, param_ders, space_ders, transversal == true ? 0 : 3);

            // TODO : 处理branch
            T curvature = space_ders[0].col(1).norm();
            T arc_length = 2 * m_angle_eps;
            if (curvature < TDEFAULT_ERROR<T>::value)
            {
                arc_length = TINFINITE<T>::value;
            }
            else
            {
                arc_length /= curvature;
            }

            // TODO : 处理step太大超过参数域
            T bigger_step_size = current_intpoints.m_is_positive_direction == false ? -arc_length : arc_length;
            T small_step = bigger_step_size;
            int type;
            T angle_diff;
            Box<T, 4> priori_enclosure;
            Eigen::Vector<T, 4> intersect_param;
            Eigen::Vector<T, 4> next_init_param;
            bool may_arrived_singular = false;
            while (false == eval_priori_enclosure(initial_box, bigger_step_size, small_step, priori_enclosure, angle_diff, type, may_arrived_singular, transversal))
            {
                // TODO: delete
                if (std::abs(bigger_step_size) < 1e-4)
                {
                    return ENUM_NURBS::NURBS_ERROR;

                }
                if (angle_diff < 0.1 || may_arrived_singular == true)
                {
                    bool flag = true;
                    for (const Box<T, 4> &singular_box : m_singular_boxes)
                    {
                        Box<T, 4> int_box;
                        if (singular_box.intersect(priori_enclosure, int_box) == true)
                        {
                            bigger_step_size /= 1.5;
                            small_step = bigger_step_size;
                            flag = false;
                            break;
                        }
                    }
                    if (flag == false)
                    {
                        continue;
                    }
                    // 计算奇异点
                    if (estimate_next_param(param, priori_enclosure, small_step, next_init_param, type) != ENUM_NURBS::NURBS_SUCCESS)
                    {
                        bigger_step_size /= 1.5;
                        small_step = bigger_step_size;
                        continue;
                    }

                    // try find singular point
                    Eigen::Vector<T, 4> pre_param;
                    if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration(m_product_box, next_init_param, pre_param))
                    {
                        if (ENUM_NURBS::NURBS_SUCCESS == intersect_singular_point_iteration(m_product_box, pre_param, intersect_param))
                        {
                            if (priori_enclosure.is_contain_point(intersect_param) == true)
                            {
                                // TODO: 整理
                                std::vector<Eigen::Matrix<T, 4, 1>> param_ders_temp;
                                std::vector<Eigen::Matrix<T, dim, 1>> space_ders_temp;
                                eval_preiamge_and_space_ders<1>(intersect_param, param_ders_temp, space_ders_temp, 3);
                                if (space_ders_temp.size() == 2)
                                {
                                    double dot_product = space_ders_temp[0].col(0).dot(space_ders_temp[1].col(0));
                                    double angele = std::acos(dot_product);
                                    Box<T, 4> singular_box;
                                    create_box(intersect_param, TDEFAULT_ERROR<T>::value * 0.01, singular_box);
                                    m_singular_boxes.push_back(singular_box);
                                    surf_surf_intersect_point<T, dim> singluar_point;// = new surf_surf_intersect_point<T, dim>();
                                    singluar_point.m_uv = intersect_param;
                                    m_left_surface->point_on_surface(intersect_param[0], intersect_param[1], singluar_point.m_point);
                                    std::vector<Eigen::Vector4<T>> next_params(4);
                                    for (int index = 0; index < 2; ++index)
                                    {
                                        Eigen::Vector4<T> param_tangent = param_ders_temp[index].col(0).normalized();
                                        Eigen::Vector4<T> initial_param1 = intersect_param + param_tangent * TDEFAULT_ERROR<T>::value;
                                        Eigen::Vector4<T> intial_param2 = intersect_param - param_tangent * TDEFAULT_ERROR<T>::value;
                                        //TOOD:处理错误, 追踪两个点
                                        intersect_point_iteration(m_product_box, initial_param1, next_params[2 * index]);
                                        intersect_point_iteration(m_product_box, intial_param2, next_params[2 * index + 1]);
                                    }
                                    T min_len = 100000;
                                    for (const Eigen::Vector4<T>& corner_param : next_params)
                                    {
                                        min_len = std::min(min_len, (corner_param - intersect_param).norm());
                                    }
                                    for (int i = 0; i < 3; ++i)
                                    {
                                        for (int j = i + 1; j < 4; ++j)
                                        {
                                            min_len = std::min(min_len, (next_params[i] - next_params[j]).norm());
                                        }
                                    }
                                    assert(min_len > 1e-6);

                                    for (int index = 0; index < 4; ++index)
                                    {
                                        //surfs_int_points_chat<T, dim> singular_chat;
                                        m_singular_chats.push_back(surfs_int_points_chat<T, dim>());
                                        surfs_int_points_chat<T, dim>& end_chat = m_singular_chats.back();
                                        // TODO: 不要写死
                                        end_chat.m_is_transversal = true;
                                        std::vector<surf_surf_intersect_point<T, dim>>& int_point = end_chat.m_inter_points;
                                        int_point.resize(2);
                                        int_point[0] = singluar_point;
                                        int_point[1].m_uv = next_params[index];
                                        m_left_surface->point_on_surface(next_params[index][0], next_params[index][1], int_point[1].m_point);

                                        //m_singular_points.push_back(int_point);

                                        //TODO:修改下面这行代码，因为这话代码的direction是写死的
                                        m_singular_point_trace_dir.push_back(index >= 2 ? '1' : '0');
                                        end_chat.m_is_positive_direction = index >= 2;
                                    }
                                    is_arrived_singular = true;
                                    return ENUM_NURBS::NURBS_SUCCESS;
                                }
                            }
                        }
                    }
                }
                bigger_step_size /= 1.5;
                small_step = bigger_step_size;
            }

            //if (false == transversal)
            //{
            //    while (angle_diff > m_angle_eps)
            //    {
            //        // TODO: delete
            //        if (std::abs(small_step) < 1e-4)
            //        {
            //            return ENUM_NURBS::NURBS_ERROR;

            //        }

            //        small_step /= 1.2;
            //        Box<T, 4> param_box = priori_enclosure;
            //        if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, param_box, small_step, priori_enclosure, type, transversal))
            //        {
            //            return ENUM_NURBS::NURBS_ERROR;
            //        }
            //        Box<T, 2> left_param_box(priori_enclosure.Min.template block<2, 1>(0, 0), priori_enclosure.Max.template block<2, 1>(0, 0));
            //        Box<T, 2> right_param_box(priori_enclosure.Min.template block<2, 1>(2, 0), priori_enclosure.Max.template block<2, 1>(2, 0));
            //        Box<T, dim> left_normal_box, right_normal_box;
            //        eval_point_interval(m_left_normal_surface, left_param_box, left_normal_box);
            //        eval_point_interval(m_right_normal_surface, right_param_box, right_normal_box);
            //        Eigen::Vector<T, dim> origin;
            //        origin.setConstant(0.0);
            //        cone<T, dim> normal_cone = point_box(left_normal_box, origin);
            //        angle_diff = normal_cone.m_angle;
            //        normal_cone = point_box(right_normal_box, origin);
            //        angle_diff = std::max(angle_diff, normal_cone.m_angle);
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


            for (int index = 0; index < m_boundary_chats.size(); ++index)
            {
                surf_surf_intersect_point<T, dim>& boundary_point = m_boundary_chats[index].m_inter_points.back();
                
                if (priori_enclosure.is_contain_point(boundary_point.m_uv))
                {
                    surf_surf_intersect_point<T, dim> current_point = current_intpoints.m_inter_points.back();
                    Eigen::Vector2<T> next_param = intersect_param.template block<2, 1>(2, 0);
                    Eigen::Vector2<T> test_param = boundary_point.m_uv.template block<2, 1>(2, 0);
                    if (is_point_on_arc(current_point.m_uv, next_param, test_param, priori_enclosure, type, current_intpoints.m_is_positive_direction))
                    {
                        is_stop = true;
                        current_intpoint.m_priori_enclosure = priori_enclosure;
                        // TODO: 处理priori_enclosure
                        int points_count = m_boundary_chats[index].m_inter_points.size();
                        for (int ipoint_index = 0; ipoint_index < points_count - 1; ++ipoint_index)
                        {
                            m_boundary_chats[index].m_inter_points[ipoint_index].m_priori_enclosure = m_boundary_chats[index].m_inter_points[ipoint_index + 1].m_priori_enclosure;
                        }

                        current_intpoints.m_inter_points.insert(current_intpoints.m_inter_points.end(), m_boundary_chats[index].m_inter_points.rbegin(), m_boundary_chats[index].m_inter_points.rend());

                        m_boundary_chats.erase(m_boundary_chats.begin() + index);
                        m_boundary_point_trace_dir.erase(m_boundary_point_trace_dir.begin() + index);
                        return ENUM_NURBS::NURBS_SUCCESS;
                    }
                }
            }
            for (int index = 0; index < m_singular_chats.size(); ++index)
            {
                surf_surf_intersect_point<T, dim>& singular_point = m_singular_chats[index].m_inter_points.back();
                if (priori_enclosure.is_contain_point(singular_point.m_uv))
                {
                    surf_surf_intersect_point<T, dim> current_point = current_intpoints.m_inter_points.back();
                    Eigen::Vector2<T> next_param = intersect_param.template block<2, 1>(2, 0);
                    Eigen::Vector2<T> test_param = singular_point.m_uv.template block<2, 1>(2, 0);
                    if (is_point_on_arc(current_point.m_uv, next_param, test_param, priori_enclosure, type, current_intpoints.m_is_positive_direction))
                    {
                        is_stop = true;
                        current_intpoint.m_priori_enclosure = priori_enclosure;
                        // TODO: 处理priori_enclosure
                        int points_count = m_singular_chats[index].m_inter_points.size();
                        for (int index = 0; index < points_count - 1; ++index)
                        {
                            m_singular_chats[index].m_inter_points[index].m_priori_enclosure = m_singular_chats[index + 1].m_inter_points[index].m_priori_enclosure;
                        }

                        current_intpoints.m_inter_points.insert(current_intpoints.m_inter_points.end(), m_singular_chats[index].m_inter_points.rbegin(), m_singular_chats[index].m_inter_points.rend());
                        m_singular_chats.erase(m_singular_chats.begin() + index);
                        m_singular_point_trace_dir.erase(m_singular_point_trace_dir.begin() + index);
                        return ENUM_NURBS::NURBS_SUCCESS;
                    }
                }
            }


            surf_surf_intersect_point<T, dim> next_intpoint;
            next_intpoint.m_priori_enclosure = priori_enclosure;
            next_intpoint.m_uv = intersect_param;
            m_left_surface->point_on_surface(intersect_param[0], intersect_param[1], next_intpoint.m_point);
            if (current_intpoints.start_point_type == POINT_TYPE::LOOP)
            {
                if (is_closed(current_intpoints, next_intpoint, type))
                {
                    is_stop = true;
                    next_intpoint = current_intpoints.m_inter_points.front();
                    current_intpoints.end_point_type == POINT_TYPE::LOOP;
                }
            }

            current_intpoints.m_inter_points.back().m_priori_enclosure = priori_enclosure;
            current_intpoints.m_inter_points.push_back(next_intpoint);
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ENUM_NURBS trace_curve_by_point(surfs_int_points_chat<T, dim>& current_intpoints, bool& is_stop)
        {
            int index = 0;
            for (; index < MAXINTERSETORPOINTNUMBER; ++index)
            {
                std::cout << index << std::endl;
                bool arrived_bound = false;
                bool is_arrived_singular = false;
                if (trace_point(current_intpoints, is_stop, is_arrived_singular) != ENUM_NURBS::NURBS_SUCCESS)
                {
                    //TODO:析构内存
                    return ENUM_NURBS::NURBS_ERROR;
                }
                if (is_arrived_singular == true || is_stop == true)
                {
                    break;
                }
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }


        //暂时不处理相交为环的情况(目前只支持bezier情况)
        ENUM_NURBS surafces_intersection2()
        {
            //先取等参线和曲面求交(非闭合)
            //using curve_type = template surface_type::iso_curve_type;
            typename surface_type::iso_curve_type curve1;
            typename surface_type::iso_curve_type curve2;
            typename surface_type::iso_curve_type curve3;
            typename surface_type::iso_curve_type curve4;

            Eigen::Vector3d point;
            m_right_surface->point_on_surface(0.5, 0.5, point);
            
            m_right_surface->get_isoparameter_curve<ENUM_DIRECTION::U_DIRECTION>(0.0, curve1);
            m_right_surface->get_isoparameter_curve<ENUM_DIRECTION::U_DIRECTION>(1.0, curve2);
            m_right_surface->get_isoparameter_curve<ENUM_DIRECTION::V_DIRECTION>(0.0, curve3);
            m_right_surface->get_isoparameter_curve<ENUM_DIRECTION::V_DIRECTION>(1.0, curve4);

           
            std::vector<Eigen::Vector<T, 3>> int_params2 = bezier_curve_int_bezier_surface(curve2, *m_left_surface);
            std::vector<Eigen::Vector<T, 3>> int_params1 = bezier_curve_int_bezier_surface(curve1, *m_left_surface);
            std::vector<Eigen::Vector<T, 3>> int_params3 = bezier_curve_int_bezier_surface(curve3, *m_left_surface);
            std::vector<Eigen::Vector<T, 3>> int_params4 = bezier_curve_int_bezier_surface(curve4, *m_left_surface);

            //TODO:去重, 加细; 另一个曲面的四条边界
            std::vector<Eigen::Vector4<T>> int_params;
            for (auto param : int_params1)
            {
                int_params.push_back(Eigen::Vector4<T>(param[1], param[2], 0.0, param[0]));
            }
            for (auto param : int_params2)
            {
                int_params.push_back(Eigen::Vector4<T>(param[1], param[2], 1.0, param[0]));
            }
            for (auto param : int_params3)
            {
                int_params.push_back(Eigen::Vector4<T>(param[1], param[2], param[0], 0));
            }
            for (auto param : int_params4)
            {
                int_params.push_back(Eigen::Vector4<T>(param[1], param[2], param[0], 1.0));
            }

            for (auto& param : int_params)
            {
                Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> left_ders;
                m_left_surface->template derivative_on_surface<1>(param[0], param[1], left_ders);
                Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_ders;
                m_right_surface->template derivative_on_surface<1>(param[2], param[3], right_ders);
                std::vector<surf_surf_intersect_point<T, dim>> current_intpoint(1);// = new surf_surf_intersect_point<T, dim>();
                current_intpoint[0].m_uv = param;
                current_intpoint[0].m_point = left_ders(0, 0);

                Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
                left_normal.normalize();
                Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
                right_normal.normalize();
                
                bool transversal = true;
                double cosTheta = left_normal.dot(right_normal);
                double theta = std::acos(cosTheta);
                if (std::abs(cosTheta) > 1.0 - KNOTS_EPS<T>::value)
                {
                    transversal = false;
                }
                bool is_stop = false;
                bool is_arrived_singular;
                surfs_int_points_chat<T, dim> boundary_points;
                boundary_points.m_inter_points = current_intpoint;
                boundary_points.m_is_transversal = transversal;
                boundary_points.m_is_positive_direction = true;
                boundary_points.start_point_type = POINT_TYPE::BOUNDARY;
                if (trace_point(boundary_points, is_stop, is_arrived_singular) == ENUM_NURBS::NURBS_SUCCESS)
                {
                    // TODO:
                    assert(is_stop == false);
                    assert(is_arrived_singular == false);
                    
                    m_boundary_chats.push_back(boundary_points);
                    m_boundary_point_trace_dir.push_back('1');
                }
                boundary_points.m_is_positive_direction = false;
                boundary_points.m_inter_points = current_intpoint;
                if (trace_point(boundary_points, is_stop, is_arrived_singular) == ENUM_NURBS::NURBS_SUCCESS)
                {
                    // TODO:
                    assert(is_stop == false);
                    assert(is_arrived_singular == false);
                    
                    m_boundary_chats.push_back(boundary_points);
                    m_boundary_point_trace_dir.push_back('0');
                }
                //else
                //{
                //    // TODO: 
                //    assert(false);
                //}
            }
            
            
            while (m_boundary_chats.empty() == false || m_singular_chats.empty() == false)
            {
                bool is_stop = false;
                if (m_singular_chats.empty() == false)
                {
                    char direction = m_singular_point_trace_dir.back();
                    surfs_int_points_chat<T, dim> current_intpoints = m_singular_chats.back();
                    m_singular_chats.pop_back();
                    m_singular_point_trace_dir.pop_back();
                    //TODO:direction 不写死
                    if (trace_curve_by_point(current_intpoints, is_stop) != ENUM_NURBS::NURBS_SUCCESS)
                    {
                        return ENUM_NURBS::NURBS_ERROR;
                    }
                    if (is_stop == false)
                    {
                        m_singular_chats.push_back(current_intpoints);
                        m_singular_point_trace_dir.push_back(direction);
                    }
                    else
                    {
                        m_result.m_int_chats.push_back(current_intpoints);
                    }
                }
                else
                {
                    surfs_int_points_chat<T, dim> current_intpoints = m_boundary_chats.back();
                    m_boundary_chats.pop_back();
                    m_boundary_point_trace_dir.pop_back();
                    if (trace_curve_by_point(current_intpoints, is_stop) != ENUM_NURBS::NURBS_SUCCESS)
                    {
                        return ENUM_NURBS::NURBS_ERROR;
                    }
                    if (is_stop == false)
                    {
                        m_boundary_chats.push_back(current_intpoints);   
                    }
                    else
                    {
                        m_result.m_int_chats.push_back(current_intpoints);
                    }
                } 
            }            
            
            
            std::vector<std::pair<Box<T, 2>, Box<T, 2>>> int_boxes;
            std::vector<std::pair<Box<T, 2>, Box<T, 2>>> sub_int_boxes;
            Box<T, 2> uv_box, st_box;
            std::vector<Box<T, 4>> exact_int_boxes;
            std::vector<Box<T, 4>> exact_int_boxes2;
            m_left_surface->get_uv_box(uv_box);
            m_right_surface->get_uv_box(st_box);
            int_boxes.push_back(std::make_pair(uv_box, st_box));
            int num_count = 0;
            int num_count2 = 0;
            while (int_boxes.empty() == false)
            {
                num_count++;
                std::pair<Box<T, 2>, Box<T, 2>> domain_box = int_boxes.back();
                int_boxes.pop_back();
                std::vector<Box<T, 2>> uv_patchs = domain_box.first.split_at_middle();
                std::vector<Box<T, 2>> st_patchs = domain_box.second.split_at_middle();
                for (int i = 0; i < 4; ++i)
                {
                    surface_type sub_surface1;
                    m_left_surface->sub_divide(uv_patchs[i], sub_surface1);
                    Box<T, dim> sub_uv_box;
                    sub_surface1.get_box(sub_uv_box);
                    for (int j = 0; j < 4; ++j)
                    {
                        surface_type sub_surface2;
                        m_right_surface->sub_divide(st_patchs[j], sub_surface2);
                        Box<T, dim> sub_st_box;
                        sub_surface2.get_box(sub_st_box);
                        Box<T, dim> intersect_box;
                        if (sub_st_box.intersect(sub_uv_box, intersect_box) == true)
                        {
                            Box<T, 4> sub_int_domian;
                            sub_int_domian.Min.template block<2, 1>(0, 0) = uv_patchs[i].Min;
                            sub_int_domian.Min.template block<2, 1>(2, 0) = st_patchs[j].Min;

                            sub_int_domian.Max.template block<2, 1>(0, 0) = uv_patchs[i].Max;
                            sub_int_domian.Max.template block<2, 1>(2, 0) = st_patchs[j].Max;

                            std::array<Box<T, 4>, 1> param_ders;

                            std::array<Box<T, 3>, 1> space_ders;
                            bool may_arrived_singular = false;
                            ENUM_NURBS code = eval_preiamge_and_space_ders<1>(sub_int_domian, param_ders, space_ders, may_arrived_singular);
                            if (code != ENUM_NURBS::NURBS_SUCCESS)
                            {
                                // TODO: 目前排除不了相切的情况
                                if ((uv_patchs[i].Min - uv_patchs[i].Max).norm() < 0.001)
                                {
                                    sub_int_boxes.push_back(std::make_pair(uv_patchs[i], st_patchs[j]));
                                }

                                int_boxes.push_back(std::make_pair(uv_patchs[i], st_patchs[j]));
                                num_count2 += 1;
                                continue;
                            }

                            Box<T, 2> uv_tangent_box(param_ders[0].Min.template block<2, 1>(0, 0), param_ders[0].Max.template block<2, 1>(0, 0));
                            cone<T, 2> uv_tangent_cone = point_box(uv_tangent_box, Eigen::Vector2<T>(0, 0));
                            if (uv_tangent_cone.is_cotain_vector(Eigen::Vector2<double>(0, 1)))
                            {
                                if ((uv_patchs[i].Min - uv_patchs[i].Max).norm() < 0.001)
                                {
                                    bool add_int_box = true;

                                    Box<T, 2> intersect_box1, intersect_box2;
                                    for (Box<T, 4> &singular_box : m_singular_boxes)
                                    {
                                        Box<T, 2> uv_singular_box(singular_box.Min.template block<2, 1>(0, 0), singular_box.Max.template block<2, 1>(0, 0));
                                        Box<T, 2> st_singular_box(singular_box.Min.template block<2, 1>(2, 0), singular_box.Max.template block<2, 1>(2, 0));
                                        if (uv_patchs[i].intersect(uv_singular_box, intersect_box1) || st_patchs[j].intersect(st_singular_box, intersect_box2))
                                        {
                                            add_int_box = false;
                                            break;
                                        }
                                    }

                                    if (add_int_box)
                                    {
                                        sub_int_boxes.push_back(std::make_pair(uv_patchs[i], st_patchs[j]));
                                    }
                                    continue;
                                }
                                int_boxes.push_back(std::make_pair(uv_patchs[i], st_patchs[j]));
                            }
                        }
                    }
                }
            }
            //std::cout << num_count << std::endl;
            //std::vector<Box<T, 2>> save_uv_box, save_st_box;

            //std::vector<Eigen::Vector3d> save_point;
            //for (auto p : sub_int_boxes)
            //{
            //    save_uv_box.push_back(p.first);
            //    save_st_box.push_back(p.second);
            //    Eigen::Vector2d param = p.first.get_middle_point();
            //    Eigen::Vector3d pp;
            //    m_left_surface->point_on_surface(param[0], param[1], pp);
            //    save_point.push_back(pp);
            //}
            //save_box(save_uv_box, "uv_box.obj");
            //save_box(save_st_box, "st_box.obj");
            //save_obj(save_point, "point.obj");

            for (std::pair<Box<T, 2>, Box<T, 2>>& int_box : sub_int_boxes)
            {
                Eigen::Vector4<T> initial_param, intersect_param;
                initial_param.template block<2, 1>(0, 0) = int_box.first.get_middle_point();
                initial_param.template block<2, 1>(2, 0) = int_box.second.get_middle_point();
                if (intersect_point_iteration(m_product_box, initial_param, intersect_param) == ENUM_NURBS::NURBS_SUCCESS)
                {
                    bool flag = is_point_on_intcurve(intersect_param, true, 0);
                    if (flag == false)
                    {
                        surf_surf_intersect_point<T, dim> current_intpoint{};
                        current_intpoint.m_uv = intersect_param;
                        m_right_surface->point_on_surface(intersect_param[2], intersect_param[3], current_intpoint.m_point);
                        surfs_int_points_chat<T, dim> loop_points;
                        loop_points.m_inter_points.push_back(current_intpoint);
                        

                        // TODO:
                        loop_points.m_is_transversal = true;
                        loop_points.m_is_positive_direction = true;
                        loop_points.start_point_type = POINT_TYPE::LOOP;
                        m_loop_chats.push_back(loop_points);

                        while (m_loop_chats.empty() == false || m_singular_chats.empty() == false)
                        {
                            bool is_stop = false;
                            if (m_singular_chats.empty() == false)
                            {
                                char direction = m_singular_point_trace_dir.back();
                                surfs_int_points_chat<T, dim> current_intpoints = m_singular_chats.back();
                                m_singular_chats.pop_back();
                                m_singular_point_trace_dir.pop_back();
                                //TODO:direction 不写死
                                if (trace_curve_by_point(current_intpoints, is_stop) != ENUM_NURBS::NURBS_SUCCESS)
                                {
                                    return ENUM_NURBS::NURBS_ERROR;
                                }
                                if (is_stop == false)
                                {
                                    m_singular_chats.push_back(current_intpoints);
                                    m_singular_point_trace_dir.push_back(direction);
                                }
                                else
                                {
                                    m_result.m_int_chats.push_back(current_intpoints);
                                }
                            }
                            else
                            {
                                surfs_int_points_chat<T, dim> current_intpoints = m_loop_chats.back();
                                m_loop_chats.pop_back();
                                if (trace_curve_by_point(current_intpoints, is_stop) != ENUM_NURBS::NURBS_SUCCESS)
                                {
                                    return ENUM_NURBS::NURBS_ERROR;
                                }
                                if (is_stop == false)
                                {
                                    m_loop_chats.push_back(current_intpoints);
                                }
                                else
                                {
                                    m_result.m_int_chats.push_back(current_intpoints);
                                }
                            }
                        }
                    }
                }
            }            
            
            return ENUM_NURBS::NURBS_SUCCESS;  
        }
    };
}