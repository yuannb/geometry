#pragma once
#include "nurbs_surface.h"
#include "cone.h"
#include "declare.h"
#include "memory"
#include <unordered_set>
#include <deque>
#include "discret.h"
namespace tnurbs
{

    template<typename curve_type, typename T = typename curve_type::Type, int dim = curve_type::dimension>
    ENUM_NURBS find_nearst_point_on_curve_inner(const curve_type *cur, const Eigen::Vector<T, dim> &point,  T &u, const Box<T, 1> &u_box , Eigen::Vector<T, dim> &nearst_point, T eps = TDEFAULT_ERROR<T>::value)
    {
        T min = u_box.Min[0];
        T max = u_box.Max[0];

        cur->point_on_curve(u, nearst_point);
        T min_length = (point - nearst_point).squaredNorm();
        Eigen::Vector<Eigen::Vector<T, dim>, 3> ders_vec;
        T distance;
        T current_u = u;
        for (int loop_index = 0; loop_index < SURFACE_ITERATE_DEEP; ++loop_index)
        {
            cur->template derivative_on_curve<2>(current_u, ders_vec);
            Eigen::Vector<T, dim> dist_vec = ders_vec[0] - point;
            distance = dist_vec.squaredNorm();
            if (distance < eps * eps)
            {
                u = current_u;
                nearst_point = ders_vec[0];
                return ENUM_NURBS::NURBS_SUCCESS;
            }

            if (min_length > distance)
            {
                min_length = distance;
                u = current_u;
                nearst_point = ders_vec[0];
            }

            T cos_angle = ders_vec[1].dot(dist_vec);
            T tanget_vec_square_len = ders_vec[1].squaredNorm();
            if (cos_angle * cos_angle < tanget_vec_square_len * distance * eps * eps)
            {
                break;
            }
            
            T next_u = current_u - cos_angle / (ders_vec[2].dot(dist_vec) + tanget_vec_square_len);
            bool is_closed_flag = cur->is_closed();
            if (is_closed_flag)
            {
                if (next_u < min)
                    next_u = max - (min - next_u);
                else if (next_u > max)
                    next_u = min + (next_u - max);
            }
            else
            {
                if (next_u < min)
                    next_u = min;
                else if (next_u > max)
                    next_u = max;
            }
            if ((next_u - current_u) * (next_u - current_u) * tanget_vec_square_len < eps * eps)
                break;
            current_u = next_u;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename curve_type, typename T = typename curve_type::Type, int dim = curve_type::dimension>
    ENUM_NURBS find_nearst_point_on_curve_inner(const curve_type *cur, const curve_mesh_helper<curve_type> &mh,  const Eigen::Vector<T, dim> &point, T &min_dis, T &u, Eigen::Vector<T, dim> &nearst_point, T eps = TDEFAULT_ERROR<T>::value)
    {
        const curve_patch<curve_type> *current = mh.root;
        Box<T, 1> domain = cur->get_domain_box();
        //1.寻找曲线极值点; 先判断是否在包围盒内部
        while (current != nullptr)
        {
            bool exit_minimal_point = true;
            //先计算点与box的最近距离
            T dis = current->box.eval_minimal_distance(point);
            if (dis > min_dis)
            {
                exit_minimal_point = false;
            }
            
            if (exit_minimal_point == true)
            {
                bool flag = current->box.is_contain_point(point);
                if (flag == false)
                {
                    //先判断目前的patch是否可能存在极值点
                    //先计算point和包围盒成的锥
                    cone<T, dim> c = point_box(current->box, point);
                    T dot_product = c.m_dir.dot(current->u_cone.m_dir);
                    T angle = 0.0;
                    if (dot_product > 1.0)
                        angle = M_PI;
                    else if (dot_product < -1.0)
                        angle = 0.0;
                    else
                        angle = std::acos(dot_product);
                    //如果两个锥相交, 则merge
                    if (angle <= c.m_angle + current->u_cone.m_angle)
                    {
                        c.merge_cone(current->u_cone);
                        if (c.m_angle < M_PI_2 / 2.0 - TDEFAULT_ERROR<T>::value)
                        {
                            exit_minimal_point = false;
                        }
                    }
                    else
                    {
                        T min_angle = angle - c.m_angle - current->u_cone.m_angle;
                        T max_angle = angle + c.m_angle + current->u_cone.m_angle;
                        if ((min_angle < M_PI_2 + TDEFAULT_ERROR<T>::value && max_angle > M_PI_2 - TDEFAULT_ERROR<T>::value) || 
                                    (min_angle < 1.5 * M_PI + TDEFAULT_ERROR<T>::value && max_angle > 1.5 * M_PI - TDEFAULT_ERROR<T>::value))
                        {
                            // exit_minimal_point = true;
                        }
                        else
                        {
                            exit_minimal_point = false;
                        }
                    }
                }
                //再判断是否是叶子节点
                if (current->left != nullptr && exit_minimal_point == true)
                {
                    current = current->left;
                    continue;
                }
            }
            if (exit_minimal_point == true)
            {
                //迭代
                T current_u = (current->u_box.Min[0] + current->u_box.Max[0]) / 2.0;
                Eigen::Vector<T, dim> current_point;
                find_nearst_point_on_curve_inner<curve_type, typename curve_type::Type, curve_type::dimension>(cur, point, current_u, domain, current_point, eps);
                if (mh.root->u_box.is_contain_point(Eigen::Vector<T, 1>(u)) == false)
                {
                    continue;
                }
                T current_dis = (current_point - point).norm();
                if (min_dis > current_dis)
                {
                    nearst_point = current_point;
                    u = current_u;
                    min_dis = current_dis;
                    if (min_dis < eps)
                        break;
                }
            }

            //寻找下一个patch
            curve_patch<curve_type> *temp = current->root;
            if (temp == nullptr)
                break;
            if (temp->left == current)
                current = temp->right;
            else
            {
                while (temp)
                {
                    if (temp->root == nullptr)
                    {
                        temp = temp->root;
                    }
                    else if (temp->root->right == temp)
                    {
                        temp = temp->root;
                    }
                    else
                    {
                        temp = temp->root->right;
                        break;
                    }
                }
                current = temp;   
            }

        }
        

        //3.和端点比较
        Box<T, 1> ends_knots;
        cur->get_u_box(ends_knots);
        Eigen::Vector<T, dim> corner_point1;
        Eigen::Vector<T, dim> corner_point2;
        cur->point_on_curve(ends_knots.Min[0], corner_point1);
        cur->point_on_curve(ends_knots.Max[0], corner_point2);
        T d1 = (point - corner_point1).norm();
        T d2 = (point - corner_point2).norm();
        if (d1 < min_dis || d2 < min_dis)
        {
            if (d1 < d2)
            {
                u = ends_knots.Min[0];
                nearst_point = corner_point1;
                min_dis = d1;
            }
            else
            {
                u = ends_knots.Max[0];
                nearst_point = corner_point2;
                min_dis = d2;
            }
        }
        
       

        return ENUM_NURBS::NURBS_SUCCESS;
    }



    template<typename curve_type>
    class curve_nearest
    {
        using T = curve_type::Type;
        static constexpr int dim = curve_type::dimension;
    private:
        T m_eps = TDEFAULT_ERROR<T>::value;
        T m_dist_eps = 10000;
        T m_angle_eps = 0.1;
        T m_chord_eps = 10000;
        Box<T, dim> *m_box = nullptr;
        Box<T, 1> *m_u_box = nullptr;
        std::vector<curve_mesh_helper<curve_type>> m_mesh_help;
        std::vector<std::unique_ptr<curve_type>> m_curves;
        /* data */
    public:
        curve_nearest() = default;
        bool set_dis_eps(T dis_eps)
        {
            if (dis_eps < TDEFAULT_ERROR<T>::value)
                return false;
            m_dist_eps = dis_eps;
            return true;
        }
        bool set_angle_eps(T angle_eps)
        {
            if (angle_eps < TDEFAULT_ERROR<T>::value)
                return false;
            m_angle_eps = angle_eps;
            return true;
        }
        bool set_chord_eps(T chord_eps)
        {
            if (chord_eps < TDEFAULT_ERROR<T>::value)
                return false;
            m_chord_eps = chord_eps;
            return true;
        }
        bool set_u_box(const Box<T, 1> *u_box)
        {
            m_u_box = u_box;
            return true;
        }
        bool set_box(const Box<T, dim> *box)
        {
            m_box = box;
            return true;
        }
        bool set_eps(T eps)
        {
            if (eps < 0.0)
                return false;
            m_eps = eps;
            return true;
        }

        ENUM_NURBS init(const curve_type *cur)
        {
            std::vector<curve_type *> sperate_curves;

            //TODO : 删除nurbs的, 范化成curve
            cur->decompose_to_nurbs(sperate_curves);
            int curves_count = sperate_curves.size();
            m_curves.reserve(curves_count);
            m_mesh_help.resize(curves_count);
            for (int index = 0; index < curves_count; ++index)
            {
                m_curves.push_back(std::unique_ptr<curve_type>(sperate_curves[index]));
                disc_curve(sperate_curves[index], m_mesh_help[index], m_eps, m_dist_eps, m_angle_eps, m_chord_eps, m_box, m_u_box);
            }
            return ENUM_NURBS::NURBS_SUCCESS;  
        }

        ENUM_NURBS find_nearst_point_on_curve(const Eigen::Vector<T, dim> &point, T &u, Eigen::Vector<T, dim> &nearst_point, T &min_dis)
        {
            Box<T, 1> ends_knots;
            m_curves[0]->get_u_box(ends_knots);
            Eigen::Vector<T, dim> corner_point1;
            Eigen::Vector<T, dim> corner_point2;
            m_curves[0]->point_on_curve(ends_knots.Min[0], corner_point1);
            m_curves[0]->point_on_curve(ends_knots.Max[0], corner_point2);
            T d1 = (point - corner_point1).norm();
            T d2 = (point - corner_point2).norm();
            
            if (d1 < d2)
            {
                nearst_point = corner_point1;
                u = ends_knots.Min[0];
                min_dis = d1;
            }
            else
            {
                nearst_point = corner_point2;
                u = ends_knots.Max[0];
                min_dis = d2;
            }
            if (min_dis < m_eps)
                return ENUM_NURBS::NURBS_SUCCESS;
            int curves_count = m_curves.size();
            for (int index = 0; index < curves_count; ++index)
            {
                find_nearst_point_on_curve_inner<curve_type, curve_type::Type, curve_type::dimension>(m_curves[index].get(), m_mesh_help[index], point, min_dis, u, nearst_point, m_eps);
            }
            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ~curve_nearest() { };
    };
    


    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS find_nearst_point_on_surface_inner(const surface_type *surf, const Eigen::Vector<T, dim> &point,  T &u, T &v, const Box<T, 2> &uv_box , Eigen::Vector<T, dim> &nearst_point, T eps = TDEFAULT_ERROR<T>::value)
    {
        T current_u = u;
        T current_v = v;
        T distance;
        T u_min = uv_box.Min[0];
        T u_max = uv_box.Max[0];
        T v_min = uv_box.Min[1];
        T v_max = uv_box.Max[1];
        surf->point_on_surface(u, v, nearst_point);
        T min_length = (nearst_point - point).squaredNorm();
        Eigen::Matrix3<Eigen::Vector<T, dim>> ders_vec;
        bool is_u_closed_flag = surf->is_u_closed();
        bool is_v_closed_flag = surf->is_v_closed();

        for (int loop_index = 0; loop_index < SURFACE_ITERATE_DEEP; ++loop_index)
        {
            surf->template derivative_on_surface<2>(current_u, current_v, ders_vec);
            Eigen::Vector<T, dim> dist_vec = ders_vec(0, 0) - point;
            distance = dist_vec.squaredNorm();
            if (distance < eps * eps)
            {
                u = current_u;
                v = current_v;
                nearst_point = ders_vec(0, 0);
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            if (distance < min_length)
            {
                min_length = distance;
                u = current_u;
                v = current_v;
                nearst_point = ders_vec(0, 0);
            }

            T cos_angle_1 = ders_vec(1, 0).dot(dist_vec);
            T cos_angle_2 = ders_vec(0, 1).dot(dist_vec);
            T tanget_vec_square_len_1 = ders_vec(1, 0).squaredNorm();
            T tanget_vec_square_len_2 = ders_vec(0, 1).squaredNorm();
            bool flag1 = cos_angle_1 * cos_angle_1 < tanget_vec_square_len_1 * distance * ANGLE_ERROR * ANGLE_ERROR;
            bool flag2 = cos_angle_2 * cos_angle_2 < tanget_vec_square_len_2 * distance * ANGLE_ERROR * ANGLE_ERROR;
            if ((flag1 == true) && (flag2 == true))
            {
                break;
            }
            Eigen::Matrix<T, 2, 2> J;
            J(0, 0) = tanget_vec_square_len_1 + dist_vec.dot(ders_vec(2, 0));
            J(1, 0) = ders_vec(1, 0).dot(ders_vec(0, 1)) + dist_vec.dot(ders_vec(1, 1));
            J(0, 1) = J(1, 0);
            J(1, 1) = tanget_vec_square_len_2 + dist_vec.dot(ders_vec(0, 2));

            Eigen::Vector<T, 2> K;
            K[0] = -1 * dist_vec.dot(ders_vec(1, 0));
            K[1] = -1 * dist_vec.dot(ders_vec(0, 1));

            double det = J.determinant();
            if (std::abs(det) < eps)
            {
                //退化矩阵, break
                break;
            }
            Eigen::JacobiSVD<Eigen::MatrixX<T>,  Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(J);
            Eigen::Vector<T, 2> delta_param = matSvd.solve(K);
            
            //TODO ：根据曲线弧长和区间长度放缩delta_param
            T next_u = current_u + delta_param[0];
            T next_v = current_v + delta_param[1];
            
            if (is_u_closed_flag)
            {
                if (next_u < u_min)
                    next_u = u_max - (u_min - next_u);
                else if (next_u > u_max)
                    next_u = u_min + (next_u - u_max);
            }
            else
            {
                if (next_u < u_min)
                    next_u = u_min;
                else if (next_u > u_max)
                    next_u = u_max;
            }

            if (is_v_closed_flag)
            {
                if (next_v < v_min)
                    next_v = v_max - (v_min - next_v);
                else if (next_v > v_max)
                    next_v = v_min + (next_v - v_max);
            }
            else
            {
                if (next_v < v_min)
                    next_v = v_min;
                else if (next_v > v_max)
                    next_v = v_max;
            }
            if (((next_u - current_u) * ders_vec(1, 0)).squaredNorm() + ((next_v - current_v) * ders_vec(0, 1)).squaredNorm() < eps * eps)
            {
                break;
            }
            current_u = next_u;
            current_v = next_v;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS find_nearst_point_on_surface_inner(const surface_type *surf, const surface_mesh_helper<surface_type> &mh, const Eigen::Vector<T, dim> &point, T &min_dis, T &u, T &v, Eigen::Vector<T, dim> &nearst_point, T eps = TDEFAULT_ERROR<T>::value)
    {
        //先找曲面的极值点
        // int count = 0;
        surface_patch<surface_type> *current = mh.root;
        Box<T, 2> domain = surf->get_domain_box();

        //2.寻找曲面极值点; 先判断是否在包围盒内部
        while (current != nullptr)
        {
            // T u_min = current->uv_box.Min[0];
            // T u_max = current->uv_box.Max[0];
            // T v_min = current->uv_box.Min[1];
            // T v_max = current->uv_box.Max[1];

            bool exit_minimal_point = true;
            //先计算点与box的最近距离
            T dis = current->box.eval_minimal_distance(point);
            if (dis > min_dis)
            {
                exit_minimal_point = false;
            }
            
            if (exit_minimal_point == true)
            {
                bool flag = current->box.is_contain_point(point);
                if (flag == false)
                {
                    //先判断目前的patch是否可能存在极值点
                    //先计算point和包围盒成的锥
                    cone<T, dim> c = point_box(current->box, point);
                    //判断u向切向锥和c锥merge之后是否大于90度, 如果小于90度, 那么此patch就没有极值点
                    cone<T, dim> temp(c);
                    T dot_product = c.m_dir.dot(current->u_cone.m_dir);
                    T angle = 0.0;
                    if (dot_product > 1.0)
                        angle = M_PI;
                    else if (dot_product < -1.0)
                        angle = 0.0;
                    else
                        angle = std::acos(dot_product);
                    //如果两个锥相交, 则merge
                    if (angle <= c.m_angle + current->u_cone.m_angle)
                    {
                        temp.merge_cone(current->u_cone);
                        if (temp.m_angle < M_PI_2 / 2.0 - TDEFAULT_ERROR<T>::value)
                        {
                            exit_minimal_point = false;
                        }
                    }
                    else
                    {
                        T min_angle = angle - c.m_angle - current->u_cone.m_angle;
                        T max_angle = angle + c.m_angle + current->u_cone.m_angle;
                        if ((min_angle < M_PI_2 + TDEFAULT_ERROR<T>::value && max_angle > M_PI_2 - TDEFAULT_ERROR<T>::value) || 
                                    (min_angle < 1.5 * M_PI + TDEFAULT_ERROR<T>::value && max_angle > 1.5 * M_PI - TDEFAULT_ERROR<T>::value))
                        {
                            // exit_minimal_point = true;
                        }
                        else
                        {
                            exit_minimal_point = false;
                        }
                    }


                    if (exit_minimal_point == true)
                    {
                        //再判断v向
                        temp = c;
                        if (angle <= c.m_angle + current->v_cone.m_angle)
                        {
                            temp.merge_cone(current->v_cone);
                            if (temp.m_angle < M_PI_2 / 2.0 - TDEFAULT_ERROR<T>::value)
                            {
                                exit_minimal_point = false;
                            }
                        }
                        else
                        {
                            T min_angle = angle - c.m_angle - current->v_cone.m_angle;
                            T max_angle = angle + c.m_angle + current->v_cone.m_angle;
                            if ((min_angle < M_PI_2 + TDEFAULT_ERROR<T>::value && max_angle > M_PI_2 - TDEFAULT_ERROR<T>::value) ||
                                (min_angle < 1.5 * M_PI + TDEFAULT_ERROR<T>::value && max_angle > 1.5 * M_PI - TDEFAULT_ERROR<T>::value) )
                            {
                                // exit_minimal_point = true;
                            }
                            else
                            {
                                exit_minimal_point = false;
                            }
                        }
                    }
                }
                //再判断是否是叶子节点
                if (current->left != nullptr && exit_minimal_point == true)
                {
                    current = current->left;
                    continue;
                }
            }
            if (exit_minimal_point == true)
            {
                //迭代
                T current_u = (current->uv_box.Min[0] + current->uv_box.Max[0]) / 2.0;
                T current_v = (current->uv_box.Min[1] + current->uv_box.Max[1]) / 2.0;
                Eigen::Vector<T, dim> current_point;
                // num += 1;
                // std::cout << num << std::endl;
                // count += 1;
                find_nearst_point_on_surface_inner<surface_type, typename surface_type::Type, surface_type::dimension>(surf, point, current_u, current_v, domain, current_point, eps);
                if (mh.root->uv_box.is_contain_point(Eigen::Vector2<T>(u, v)) == false)
                {
                    continue;
                }
                T current_dis = (current_point - point).norm();
                if (min_dis > current_dis)
                {
                    nearst_point = current_point;
                    u = current_u;
                    v = current_v;
                    min_dis = current_dis;
                    if (min_dis < eps)
                        break;
                }
            }

            //寻找下一个patch
            surface_patch<surface_type> *temp = current->root;
            if (temp == nullptr)
                break;
            if (temp->left == current)
                current = temp->right;
            else
            {
                while (temp)
                {
                    if (temp->root == nullptr)
                    {
                        temp = temp->root;
                    }
                    else if (temp->root->right == temp)
                    {
                        temp = temp->root;
                    }
                    else
                    {
                        temp = temp->root->right;
                        break;
                    }
                }
                current = temp;   
            }

        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename surface_type>
    class surface_nearest
    {
        using T = surface_type::Type;
        using curve_type = surface_type::iso_curve_type;
        static constexpr int dim = surface_type::dimension;
    private:
        T m_eps = TDEFAULT_ERROR<T>::value;
        T m_dist_eps = 10000;
        T m_angle_eps = 0.1;
        T m_chord_eps = 10000;
        Box<T, dim> *m_box = nullptr;
        Box<T, 2> *m_uv_box = nullptr;

        Eigen::MatrixX<surface_mesh_helper<surface_type>> m_mesh_help;
        Eigen::MatrixX<std::unique_ptr<surface_type>> m_surfaces;

        std::vector<T> m_u_params, m_v_params;
        std::vector<std::unique_ptr<curve_type>> m_u_curves;
        std::vector<std::unique_ptr<curve_type>> m_v_curves;

        std::vector<curve_mesh_helper<curve_type>> m_u_curve_mesh_help;
        std::vector<curve_mesh_helper<curve_type>> m_v_curve_mesh_help;
        /* data */
    public:
        surface_nearest() = default;
        bool set_dis_eps(T dis_eps)
        {
            if (dis_eps < TDEFAULT_ERROR<T>::value)
                return false;
            m_dist_eps = dis_eps;
            return true;
        }
        bool set_angle_eps(T angle_eps)
        {
            if (angle_eps < TDEFAULT_ERROR<T>::value)
                return false;
            m_angle_eps = angle_eps;
            return true;
        }
        bool set_chord_eps(T chord_eps)
        {
            if (chord_eps < TDEFAULT_ERROR<T>::value)
                return false;
            m_chord_eps = chord_eps;
            return true;
        }
        bool set_uv_box(const Box<T, 2> *uv_box)
        {
            m_uv_box = uv_box;
            return true;
        }
        bool set_box(const Box<T, dim> *box)
        {
            m_box = box;
            return true;
        }
        bool set_eps(T eps)
        {
            if (eps < 0.0)
                return false;
            m_eps = eps;
            return true;
        }

        ENUM_NURBS init(const surface_type *surf)
        {
            Eigen::MatrixX<surface_type *> sperate_surfaces;

            //TODO : 删除nurbs的, 范化成curve
            surf->decompose_to_nurbs(sperate_surfaces);
            int v_count = sperate_surfaces.rows();
            int u_count = sperate_surfaces.cols();
            m_mesh_help.resize(v_count, u_count);
            m_surfaces.resize(v_count, u_count);
            for (int v_index = 0; v_index < v_count; ++v_index)
            {
                for (int u_index = 0; u_index < u_count; ++u_index)
                {
                    m_surfaces(v_index, u_index) = std::move(std::unique_ptr<surface_type>(sperate_surfaces(v_index, u_index)));
                    disc_surface(sperate_surfaces(v_index, u_index), m_mesh_help(v_index, u_index), m_eps, m_dist_eps, m_angle_eps, m_chord_eps, m_box, m_uv_box);
                }
            }

            //再计算曲线上的极值点
            //先将c0处的等参线取出
            std::vector<curve_type *> u_nurbs_curves, v_nurbs_curves;
            surf->get_c0_isoparameter_curve(u_nurbs_curves, m_u_params, v_nurbs_curves, m_v_params);
            int u_curves_count = u_nurbs_curves.size();
            int v_curves_count = v_nurbs_curves.size();
            m_u_curves.reserve(u_curves_count);
            m_v_curves.reserve(v_curves_count);
            m_u_curve_mesh_help.resize(u_curves_count);
            m_v_curve_mesh_help.resize(v_curves_count);
            Box<T, 1> *u_box = nullptr;
            Box<T, 1> temp_u_box;
            if (m_uv_box != nullptr)
            {
                temp_u_box = Box<T, 1>( Eigen::Vector<T, 1> (m_uv_box->Min[0]), Eigen::Vector<T, 1> { m_uv_box->Max[0] });
                u_box = &temp_u_box;
            }

             
            
            for (int index = 0; index < u_curves_count; ++index)
            {
                m_u_curves.push_back(std::unique_ptr<curve_type>(u_nurbs_curves[index]));
                disc_curve(u_nurbs_curves[index], m_u_curve_mesh_help[index], m_eps, m_dist_eps, m_angle_eps, m_chord_eps, m_box, u_box);
            }
            Box<T, 1> *v_box = nullptr;
            Box<T, 1> temp_v_box;
            if (m_uv_box != nullptr)
            {
                temp_v_box = Box<T, 1>(Eigen::Vector<T, 1>{ m_uv_box->Min[1] } , Eigen::Vector<T, 1> { m_uv_box->Max[1] } );
                v_box = &temp_v_box;
            }
            for (int index = 0; index < v_curves_count; ++index)
            {
                m_v_curves.push_back(std::unique_ptr<curve_type>(v_nurbs_curves[index]));
                disc_curve(v_nurbs_curves[index], m_v_curve_mesh_help[index], m_eps, m_dist_eps, m_angle_eps, m_chord_eps, m_box, v_box);
            }

            return ENUM_NURBS::NURBS_SUCCESS;  
        }

        ENUM_NURBS find_nearst_point_on_surface(const Eigen::Vector<T, dim> &point, T &u, T &v, Eigen::Vector<T, dim> &nearest_point, T &min_distance)
        {
            Box<T, 2> uv_box;
            m_surfaces(0, 0)->get_uv_box(uv_box);
            Eigen::Vector<T, dim> temp_point;
            std::vector<Eigen::Vector<T, dim>> corner_points;
            m_surfaces(0, 0)->point_on_surface(uv_box.Min[0], uv_box.Min[1], temp_point);
            corner_points.push_back(temp_point);
            m_surfaces(0, 0)->point_on_surface(uv_box.Min[0], uv_box.Max[1], temp_point);
            corner_points.push_back(temp_point);
            m_surfaces(0, 0)->point_on_surface(uv_box.Max[0], uv_box.Min[1], temp_point);
            corner_points.push_back(temp_point);
            m_surfaces(0, 0)->point_on_surface(uv_box.Max[0], uv_box.Max[1], temp_point);
            corner_points.push_back(temp_point);

            nearest_point = corner_points[0];
            u = uv_box.Min[0];
            v = uv_box.Min[1];
            min_distance = (corner_points[0] - point).norm();

            for (int j = 1; j < 4; ++j)
            {
                T distance = (corner_points[j] - point).norm();
                if (distance < min_distance)
                {
                    min_distance = distance;
                    nearest_point = corner_points[j];
                    u = j > 1 ? uv_box.Max[0] : uv_box.Min[0];
                    v = j % 2 == 0 ? uv_box.Min[1] : uv_box.Max[1];
                }
            }

            int v_count = m_surfaces.rows();
            int u_count = m_surfaces.cols();

            for (int v_index = 0; v_index < v_count; ++v_index)
            {
                for (int u_index = 0; u_index < u_count; ++u_index)
                {
                    surface_type* surf = nullptr;
                    // surface_mesh_helper<surface_type> temp_mh; 
                    //Eigen::Vector<T, dim> tmep_point; 
                    Eigen::Matrix<T, dim, 1> tmep_point;
                    T temp_min_dis, temp_u, temp_v;
                    Eigen::Vector<T, dim> temp_nearst_point;
                    find_nearst_point_on_surface_inner<surface_type, surface_type::Type, surface_type::dimension>(m_surfaces(v_index, u_index).get(), m_mesh_help(v_index, u_index), point, min_distance, u, v, nearest_point, m_eps);
                    //find_nearst_point_on_surface_inner<surface_type::Type, surface_type::dimension>(point);
                    if (min_distance < m_eps)
                        return ENUM_NURBS::NURBS_SUCCESS;
                }
            }

            int u_curves_count = m_u_curves.size();
            
            for (int index = 0; index < u_curves_count; ++index)
            {
                T temp_min_dis = min_distance;
                find_nearst_point_on_curve_inner<curve_type, curve_type::Type, curve_type::dimension>(m_u_curves[index].get(), m_u_curve_mesh_help[index], point, temp_min_dis, v, nearest_point, m_eps);
                if (temp_min_dis < min_distance)
                {
                    min_distance = temp_min_dis;
                    u = m_u_params[index];
                    if (min_distance < m_eps)
                        return ENUM_NURBS::NURBS_SUCCESS;
                }
            }
            int v_curves_count = m_v_curves.size();
            
            for (int index = 0; index < v_curves_count; ++index)
            {
                T temp_min_dis = min_distance;
                find_nearst_point_on_curve_inner<curve_type, curve_type::Type, curve_type::dimension>(m_v_curves[index].get(), m_v_curve_mesh_help[index], point, temp_min_dis, u, nearest_point, m_eps);
                if (temp_min_dis < min_distance)
                {
                    min_distance = temp_min_dis;
                    v = m_v_params[index];
                    if (min_distance < m_eps)
                        return ENUM_NURBS::NURBS_SUCCESS;
                }
            }



            return ENUM_NURBS::NURBS_SUCCESS;
        }

        ~surface_nearest() { };
    };
    
    
}