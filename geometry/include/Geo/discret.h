#pragma once
#include "nurbs_surface.h"
#include "cone.h"
#include "declare.h"
#include "memory"
#include <unordered_set>
#include <unordered_map>
#include <deque>
namespace tnurbs
{
    //计算nurbs的切向向锥
    template<typename T, int dim>
    ENUM_NURBS eval_nurbs_tangnet_cone(const Eigen::Matrix<T, dim, Eigen::Dynamic> &control_points,
        cone<T, dim> &tangent_cone)
    {
        int points_count = control_points.cols();
        Eigen::Vector<T, dim> dir = control_points.col(1) - control_points.col(0);
        dir.normalize();
        tangent_cone.m_angle = 0.0;
        tangent_cone.m_dir = dir;
        for (int u_index = 2; u_index < points_count; ++u_index)
        {
            dir = control_points.col(u_index) - control_points.col(u_index - 1);
            dir.normalize();
            tangent_cone.merge_vector(dir);
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    //计算nurbs surface的切向向锥
    template<typename T, int dim>
    ENUM_NURBS eval_nurbs_tangnet_cone(const std::vector<Eigen::Matrix<T, dim, Eigen::Dynamic>> &control_points,
        std::vector<cone<T, dim>> &tangent_cones)
    {
        int v_points_count = control_points.size();
        tangent_cones.reserve(v_points_count);
        for (int v_index = 0; v_index < v_points_count; ++v_index)
        {
            cone<T, dim> tangent_cone;
            eval_nurbs_tangnet_cone(control_points[v_index], tangent_cone);
            tangent_cones.push_back(tangent_cone);
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    //此类仅在curve_mesh_helper中使用, 内存也在其中管理
    template<typename curve_type>
    struct curve_patch
    {
        using T = typename curve_type::Type;
        static constexpr int dim = curve_type::dimension;

        std::array<int, 2> point_index;
        const curve_type *cur;
        Box<T, dim> box;
        Box<T, 1> u_box;
        curve_patch *left;
        curve_patch *right;
        curve_patch *root;
        cone<T, dim> u_cone;
        curve_patch() : cur(nullptr), left(nullptr), right(nullptr), root(nullptr) { };
        ~curve_patch () {  }
    };

    template<typename curve_type>
    struct curve_mesh_helper
    {
        using T = typename curve_type::Type;
        static constexpr int dim = curve_type::dimension;

        curve_mesh_helper() = default;
        ~curve_mesh_helper () 
        {    
            curve_patch<curve_type> *current = root;
            while (current)
            {
                if (current->left == nullptr && current->right == nullptr)
                {
                    curve_patch<curve_type> *temp = current;
                    current = current->root;
                    if (current == nullptr)
                    {
                        delete temp;
                        break;
                    }

                    if (temp == current->left)
                        current->left = nullptr;
                    else
                        current->right = nullptr;
                    delete temp;
                }
                else if (current->left != nullptr)
                    current = current->left;
                else
                    current = current->right;
                // if (current->left)
                // {
                //     current = current->left;
                //     continue;
                // }
                // curve_patch<curve_type> *temp = current->root;
                
                // //寻找下一个节点
                // if (temp == nullptr)
                // {
                //     delete current;
                //     break;
                // }
                // if (temp->left == current)
                // {
                //     temp->left = nullptr;
                // }
                // else
                // {
                //     temp->right = nullptr;
                // }
                // if (temp->right != nullptr)
                //     current = temp->right;
                // else
                // {
                //     //找到下一个节点
                //     while (temp)
                //     {
                //         if (temp->root == nullptr)
                //         {
                //             current = temp;
                //             temp = temp->root;
                //             delete current;
                //         }
                //         else if (temp->root->right == temp)
                //         {
                //             current = temp;
                //             temp = temp->root;
                //             delete current;
                //         }
                //         else
                //         {
                //             current = temp;
                //             temp = temp->root->right;
                //             delete current;
                //             break;
                //         }
                //     }
                //     current = temp;   
                // }
            }
            
        }
        std::vector<Eigen::Vector<Eigen::Vector<T, dim>, 2>> ders;
        // std::vector<std::array<int, 2>> point_indexs; //TODO : 删除
        curve_patch<curve_type> *root;

    };


    template<typename curve_type, typename T = typename curve_type::Type, int dim = curve_type::dimension>
    bool check_patch(curve_patch<curve_type> &pat, std::vector<Eigen::Vector<Eigen::Vector<T, dim>, 2>> &ders, T dis_eps, T angel_eps, T chord_eps)
    {
        Box<T, 1> &u_box = pat.u_box;
        T mid = (u_box.Min[0] + u_box.Max[0]) / 2.0;

        Eigen::Vector<Eigen::Vector<T, dim>, 2> p0_ders = ders[pat.point_index[0]];

        Eigen::Vector<Eigen::Vector<T, dim>, 2> p1_ders;
        pat.cur->template derivative_on_curve<1>(mid, p1_ders);

        Eigen::Vector<Eigen::Vector<T, dim>, 2> p2_ders = ders[pat.point_index[1]];

        T u_dis = (p0_ders[0] - p1_ders[0]).norm() + (p2_ders[0] - p1_ders[0]).norm();

        
        Eigen::Vector<T, dim> u_dir = p1_ders(1);
        u_dir.normalize();

        Eigen::Vector<T, dim> temp = p1_ders[0] - p0_ders[0];
        temp -= (temp.dot(u_dir) * u_dir);
        T chord_dis = temp.norm();
        temp = p1_ders[0] - p2_ders[0];
        temp -= (temp.dot(u_dir) * u_dir);
        chord_dis = std::max(chord_dis, temp.norm());

        

        curve_type cur;
        pat.cur->sub_divide(u_box, cur);
        cone<T, dim> u_cone;
        eval_nurbs_tangnet_cone<T, dim>(cur.get_nonhomo_control_points(), u_cone);
        pat.u_cone = u_cone;

        if (u_cone.m_angle > angel_eps && u_box.Max[0] - u_box.Min[0] > TDEFAULT_ERROR<T>::value)
        {
            ders.push_back(p1_ders);
            return false;
        }

        if (chord_dis > chord_eps || u_dis > dis_eps)
        {
            ders.push_back(p1_ders);
            return false;
        }
            
        return true;
    }

    template<typename curve_type, typename T = typename curve_type::Type, int dim = curve_type::dimension>
    ENUM_NURBS disc_curve(const curve_type *cur, curve_mesh_helper<curve_type> &mh, T eps = TDEFAULT_ERROR<T>::value, T dist_eps = TDEFAULT_ERROR<T>::value, 
        T angle_eps = TDEFAULT_ERROR<T>::value, T chord_eps = TDEFAULT_ERROR<T>::value, const Box<T, dim> *box = nullptr, const Box<T, 1> *u_box = nullptr)
    {
        curve_patch<curve_type> *root = new curve_patch<curve_type>();
        mh.ders.clear();
        mh.root = root;

        root->root = nullptr;
        root->left = nullptr;
        root->right = nullptr;
        root->cur = cur;
        if (u_box != nullptr)
        {
            //不检测box的合法性
            root->u_box = *u_box;
        }
        else
        {
            cur->get_u_box(root->u_box);
        }
        if (box != nullptr)
        {
            //不检测box的合法性
            root->box = *box;
        }
        else
        {
            curve_type temp_cur(*cur);
            temp_cur.sub_divide(root->u_box);
            temp_cur.get_box(root->box);
        }


        //计算2个角点
        Eigen::Vector2<Eigen::Vector<T, dim>> left, right;
        cur->template derivative_on_curve<1>(root->u_box.Min[0], left);
        cur->template derivative_on_curve<1>(root->u_box.Max[0], right);

        mh.ders.push_back(left);
        mh.root->point_index[0] = 0;
        mh.ders.push_back(right);
        mh.root->point_index[1] = 1;

        curve_patch<curve_type> *current = root;

        do
        {
            bool flag = check_patch<curve_type, typename curve_type::Type, curve_type::dimension>(*current, mh.ders, dist_eps, angle_eps, chord_eps);
            
            if (flag == true)
            {
                // mh.point_indexs.push_back(current->point_index); //TODO ： 待删除
                curve_patch<curve_type> *temp = current->root;
                if (temp == nullptr)
                    break;
                if (temp->left == current)
                    current = temp->right;
                else
                {
                    //找到下一个节点
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
            else
            {
                int m_ders_size = mh.ders.size() - 1;
                curve_patch<curve_type> *left = new curve_patch<curve_type>();
                curve_patch<curve_type> *right = new curve_patch<curve_type>();
                left->u_box.Min = current->u_box.Min;
                left->u_box.Max = (current->u_box.Max + current->u_box.Min) / 2.0;
                right->u_box.Min = left->u_box.Max;
                right->u_box.Max = current->u_box.Max;
                left->point_index = std::array<int, 2> { current->point_index[0], m_ders_size };
                right->point_index = std::array<int, 2> { m_ders_size, current->point_index[1] };
                
                current->left = left;
                current->right = right;
                right->root = current;
                left->root = current;
                right->cur = cur;
                left->cur = cur;            
                //计算box
                curve_type temp_curve(*cur);
                temp_curve.sub_divide(left->u_box);
                temp_curve.get_box(left->box);
                
                curve_type temp_curve2(*cur);
                temp_curve2.sub_divide(right->u_box);
                temp_curve2.get_box(right->box);

                current = left;
            }
        } while (current != nullptr);
        
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    //此类仅在surface_mesh_helper中使用, 内存也在其中管理
    template<typename surface_type>
    struct surface_patch
    {
        using T = typename surface_type::Type;
        static constexpr int dim = surface_type::dimension;

        std::array<int, 4> point_index;
        const surface_type *surface;
        Box<T, dim> box;
        Box<T, 2> uv_box;
        surface_patch *left;
        surface_patch *right;
        surface_patch *root;
        cone<T, dim> u_cone;
        cone<T, dim> v_cone;

        //离散时使用, indexs记录分割边的所以点
        std::vector<int> indexs;
        ENUM_DIRECTION split_direction;

        surface_patch() : surface(nullptr), left(nullptr), right(nullptr), root(nullptr) { };
        ~surface_patch () {  }
    };


    template<typename surface_type>
    struct surface_mesh_helper
    {
        using T = typename surface_type::Type;
        static constexpr int dim = surface_type::dimension;

        surface_mesh_helper() = default;
        ~surface_mesh_helper () 
        {
            surface_patch<surface_type> *current = root;
            while (current)
            {
                if (current->left == nullptr && current->right == nullptr)
                {
                    surface_patch<surface_type> *temp = current;
                    current = current->root;
                    if (current == nullptr)
                    {
                        delete current;
                        current = nullptr;
                        break;
                    }

                    if (temp == current->left)
                        current->left = nullptr;
                    else
                        current->right = nullptr;
                    delete temp;
                }
                else if (current->left != nullptr)
                    current = current->left;
                else
                    current = current->right;
            }
            
        }

        std::vector<Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2>> ders;
        // std::vector<std::array<int, 4>> point_indexs;
        surface_patch<surface_type> *root;

    };


    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    bool check_patch(surface_patch<surface_type> &pat, std::vector<Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2>> &ders, 
        ENUM_DIRECTION &direction, std::vector<Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2>> &new_ders,  T dis_eps, T angel_eps, T chord_eps)
    {
        new_ders.clear();
        Box<T, 2> &uv_box = pat.uv_box;
        Eigen::Vector<T, 2> ul_mid((uv_box.Min[0] + uv_box.Max[0]) / 2.0, uv_box.Min[1]);
        Eigen::Vector<T, 2> uh_mid(ul_mid[0], uv_box.Max[1]);
        Eigen::Vector<T, 2> vl_mid(uv_box.Min[0], (uv_box.Min[1] + uv_box.Max[1]) / 2.0);
        Eigen::Vector<T, 2> vh_mid(uv_box.Max[0], vl_mid[1]);
        Eigen::Vector<T, 2> uv_mid(ul_mid[0], vl_mid[1]);


        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> p0_ders = ders[pat.point_index[0]];

        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> p1_ders;
        pat.surface->template derivative_on_surface<1>(ul_mid[0], ul_mid[1], p1_ders);

        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> p2_ders = ders[pat.point_index[1]];

        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> p3_ders;
        pat.surface->template derivative_on_surface<1>(vh_mid[0], vh_mid[1], p3_ders);

        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> p4_ders = ders[pat.point_index[2]];

        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> p5_ders;
        pat.surface->template derivative_on_surface<1>(uh_mid[0], uh_mid[1], p5_ders);

        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> p6_ders = ders[pat.point_index[3]];

        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> p7_ders;
        pat.surface->template derivative_on_surface<1>(vl_mid[0], vl_mid[1], p7_ders);

        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> p8_ders;
        pat.surface->template derivative_on_surface<1>(uv_mid[0], uv_mid[1], p8_ders);

        std::vector<Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2>> current_points{p0_ders, p1_ders, p2_ders, p3_ders, p4_ders, p5_ders, p6_ders, p7_ders, p8_ders};
        T d01 = (p0_ders(0, 0) - p1_ders(0, 0)).norm() + (p2_ders(0, 0) - p1_ders(0, 0)).norm();
        T d12 = (p3_ders(0, 0) - p2_ders(0, 0)).norm() + (p4_ders(0, 0) - p3_ders(0, 0)).norm();
        T d23 = (p5_ders(0, 0) - p4_ders(0, 0)).norm() + (p6_ders(0, 0) - p5_ders(0, 0)).norm();
        T d30 = (p7_ders(0, 0) - p6_ders(0, 0)).norm() + (p0_ders(0, 0) - p7_ders(0, 0)).norm();

        T d = (p0_ders(0, 0) - p8_ders(0, 0)).norm();
        for (int index = 1; index < 8; ++index)
        {
            T current_d = (current_points[index](0, 0) - p8_ders(0, 0)).norm();
            if (d < current_d)
                d = current_d;
        }

        T chord_dis = 0;
        Eigen::Vector<T, dim> u_dir = p8_ders(1, 0);
        u_dir.normalize();
        Eigen::Vector<T, dim> v_dir = p8_ders(0, 1);
        //将其与u_dir正交
        v_dir -= v_dir.dot(u_dir) * u_dir;
        v_dir.normalize();

        for (int index = 0; index < 8; ++index)
        {
            Eigen::Vector<T, dim> temp = current_points[index](0, 0) - p8_ders(0,0);
            T u_coeff = temp.dot(u_dir);
            T v_coeff = temp.dot(v_dir);
            temp -= (u_coeff * u_dir + v_coeff * v_dir);
            T current_dis = std::abs(temp.norm());
            if (current_dis > chord_dis)
                chord_dis = current_dis;
        }

        

        surface_type surf;
        pat.surface->sub_divide(uv_box, surf);
        std::vector<cone<T, dim>> u_cones, v_cones;
        eval_nurbs_tangnet_cone<T, dim>(surf.get_nonhomo_control_points(), u_cones);
        cone<T, dim> u_cone(u_cones[0]);
        int u_cones_count = u_cones.size();
        for (int index = 1; index < u_cones_count; ++index)
        {
            u_cone.merge_cone(u_cones[index]);
        }
        pat.u_cone = u_cone;
        
        surf.reverse_uv();
        eval_nurbs_tangnet_cone<T, dim>(surf.get_nonhomo_control_points(), v_cones);
        cone<T, dim> v_cone(v_cones[0]);
        int v_cones_count = v_cones.size();
        for (int index = 1; index < v_cones_count; ++index)
        {
            v_cone.merge_cone(v_cones[index]);
        }
        pat.v_cone = v_cone;

        T u_dis = d01 + d23;
        T v_dis = d12 + d30;

        bool flag1 = u_dis > 100 * v_dis;
        bool flag2 = v_dis > 100 * u_dis;

        if ((u_cone.m_angle > angel_eps && u_cone.m_angle > v_cone.m_angle) && uv_box.Max[0] - uv_box.Min[0] > TDEFAULT_ERROR<T>::value && flag2 == false)
        {
            direction =  ENUM_DIRECTION::U_DIRECTION;
            new_ders.push_back(p1_ders);
            new_ders.push_back(p5_ders);
            return false;
        }
        else if (v_cone.m_angle > angel_eps && uv_box.Max[1] - uv_box.Min[1] > TDEFAULT_ERROR<T>::value && flag1 == false)
        {
            direction = ENUM_DIRECTION::V_DIRECTION;
            new_ders.push_back(p3_ders);
            new_ders.push_back(p7_ders);
            return false;
        }

        if (chord_dis > chord_eps)
        {
            if (u_dis >= v_dis)
            {
                direction = ENUM_DIRECTION::U_DIRECTION;
                new_ders.push_back(p1_ders);
                new_ders.push_back(p5_ders);  
                return false;
            }
            else
            {
                direction = ENUM_DIRECTION::V_DIRECTION;
                new_ders.push_back(p3_ders);
                new_ders.push_back(p7_ders);
                return false;
            }
        }
        if (d01 > dis_eps || d12 > dis_eps || d23 > dis_eps || d30 > dis_eps || d > dis_eps)
        {
            if (u_dis >= v_dis)
            {
                direction = ENUM_DIRECTION::U_DIRECTION;
                new_ders.push_back(p1_ders);
                new_ders.push_back(p5_ders);
                return false;
            }
            else
            {
                direction = ENUM_DIRECTION::V_DIRECTION;
                new_ders.push_back(p3_ders);
                new_ders.push_back(p7_ders);
                return false;
            }
        }
            
        return true;
    }


    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS disc_surface(const surface_type *surf, surface_mesh_helper<surface_type> &mh, T eps = TDEFAULT_ERROR<T>::value, T dist_eps = TDEFAULT_ERROR<T>::value, 
        T angle_eps = TDEFAULT_ERROR<T>::value, T chord_eps = TDEFAULT_ERROR<T>::value, const Box<T, dim> *box = nullptr, const Box<T, 2> *uv_box = nullptr)
    {
        surface_patch<surface_type> *root = new surface_patch<surface_type>();
        mh.ders.clear();
        mh.root = root;

        root->root = nullptr;
        root->left = nullptr;
        root->right = nullptr;
        root->surface = surf;
        if (uv_box != nullptr)
        {
            //不检测box的合法性
            root->uv_box = *uv_box;
        }
        else
        {
            surf->get_uv_box(root->uv_box);
        }
        if (box != nullptr)
        {
            //不检测box的合法性
            root->box = *box;
        }
        else
        {
            surface_type temp_surf(*surf);
            temp_surf.sub_divide(root->uv_box);
            temp_surf.get_box(root->box);
        }


        //计算4个角点
        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> left_bottom, left_top, right_bottom, right_top;
        surf->template derivative_on_surface<1>(root->uv_box.Min[0], root->uv_box.Min[1], left_bottom);
        surf->template derivative_on_surface<1>(root->uv_box.Max[0], root->uv_box.Min[1], right_bottom);
        surf->template derivative_on_surface<1>(root->uv_box.Max[0], root->uv_box.Max[1], right_top);
        surf->template derivative_on_surface<1>(root->uv_box.Min[0], root->uv_box.Max[1], left_top);

        mh.ders.push_back(left_bottom);
        mh.root->point_index[0] = 0;
        mh.ders.push_back(right_bottom);
        mh.root->point_index[1] = 1;
        mh.ders.push_back(right_top);
        mh.root->point_index[2] = 2;
        mh.ders.push_back(left_top);
        mh.root->point_index[3] = 3;

        surface_patch<surface_type> *curent = root;
        // int num = 0;
        // int countt = 0;
        do
        {
            // std::cout << countt << std::endl;
            // ++countt;
            std::vector<Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2>> new_ders;
            ENUM_DIRECTION dir;
            bool flag = check_patch<surface_type, typename surface_type::Type, surface_type::dimension>(*curent, mh.ders, dir, new_ders, dist_eps, angle_eps, chord_eps);
            if (flag == true)
            {
                // num += 1;
                // mh.point_indexs.push_back(curent->point_index);
                surface_patch<surface_type> *temp = curent->root;
                if (temp == nullptr)
                    break;
                if (temp->left == curent)
                    curent = temp->right;
                else
                {
                    //找到下一个节点
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
                    curent = temp;   
                }
            }
            else
            {
                mh.ders.insert(mh.ders.end(), new_ders.begin(), new_ders.end());
                int m_ders_size = mh.ders.size() - 1;
                surface_patch<surface_type> *left = new surface_patch<surface_type>();
                surface_patch<surface_type> *right = new surface_patch<surface_type>();
                if (dir == ENUM_DIRECTION::U_DIRECTION)
                {
                    left->uv_box.Min = curent->uv_box.Min;
                    left->uv_box.Max[1] = curent->uv_box.Max[1];
                    left->uv_box.Max[0] = (curent->uv_box.Min[0] + curent->uv_box.Max[0]) / 2.0;
                    left->point_index = std::array<int, 4>{ curent->point_index[0], m_ders_size - 1, m_ders_size, curent->point_index[3] };

                    right->uv_box.Max = curent->uv_box.Max;
                    right->uv_box.Min[1] = curent->uv_box.Min[1];
                    right->uv_box.Min[0] = left->uv_box.Max[0];
                    right->point_index = std::array<int, 4>{ m_ders_size - 1, curent->point_index[1], curent->point_index[2], m_ders_size };
                    curent->split_direction = ENUM_DIRECTION::U_DIRECTION;
                }
                else
                {
                    left->uv_box.Min = curent->uv_box.Min;
                    left->uv_box.Max[0] = curent->uv_box.Max[0];
                    left->uv_box.Max[1] = (curent->uv_box.Min[1] + curent->uv_box.Max[1]) / 2.0;
                    left->point_index = std::array<int, 4>{ curent->point_index[0], curent->point_index[1], m_ders_size - 1, m_ders_size };
                    
                    right->uv_box.Max = curent->uv_box.Max;
                    right->uv_box.Min[0] = curent->uv_box.Min[0];
                    right->uv_box.Min[1] = left->uv_box.Max[1];
                    right->point_index = std::array<int, 4>{ m_ders_size, m_ders_size - 1, curent->point_index[2], curent->point_index[3] };
                    curent->split_direction = ENUM_DIRECTION::V_DIRECTION;
                }
                curent->left = left;
                curent->right = right;
                right->root = curent;
                left->root = curent;
                right->surface = surf;
                left->surface = surf;            
                //计算box
                surface_type temp_surf(*surf);
                temp_surf.sub_divide(left->uv_box);
                temp_surf.get_box(left->box);
                
                surface_type temp_surf2(*surf);
                temp_surf2.sub_divide(right->uv_box);
                temp_surf2.get_box(right->box);

                curent = left;
            }
        } while (curent != nullptr);
        // std::cout << "disc " << num << std::endl;
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    /// @brief mesh 类
    /// @tparam T double, float ...
    /// @tparam n (0 : 表示点; 1 : 表示曲线; 2 : 表示曲面; 其余的暂未实现)
    /// @tparam dim  维数
    template<int n, typename T = double, int dim = 3>
    class mesh
    {
    public:
        // \phi(x_1, x_2, ..., x_n)为参数方程; m_ders[0] 表示点, m_ders[i] = \partial{\phi} / \partial{x_i}
        std::vector<Eigen::Matrix<Eigen::Vector<T, dim>, n, n>> m_ders;
        // 表示点的连接关系
        std::vector<std::vector<int>> m_indexs;
    public:
        mesh() = default;
        ~mesh() { };
    };

    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS remove_multiple_point(surface_mesh_helper<surface_type> &mh)
    {
        struct help
        {
            int m_index;
            T m_eps;
            Eigen::Vector<T, dim> m_point;
            help(const Eigen::Vector<T, dim> &point, int index, T eps = TDEFAULT_ERROR<T>::value) { m_index = index; m_point = point; m_eps = eps; }
            ~help() { }
            const bool operator==(const help &right) const
            {
                if ((m_point - right.m_point).squaredNorm() < m_eps)
                    return true;
                return false;
            }
        };
        class hash_help {
        public:
            size_t operator()(const help& hero)const {
                std::string temp;
               
                for (int index = 0; index < dim; ++index)
                {
                    temp += std::to_string(hero.m_point[index]);
                }
                return std::hash<std::string>()(temp);
            }
        };
        std::unordered_set<help, hash_help> new_ders;
        std::vector<Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2>> result_vec;
        std::unordered_map<int, int> index_map;
        int difference = 0;
        int mh_ders_size = mh.ders.size();
        result_vec.reserve(mh_ders_size);
        for (int index = 0; index < mh_ders_size; ++index)
        {
            help temp(mh.ders[index](0, 0), index - difference);
            auto it = new_ders.find(temp);
            if (it == new_ders.end())
            {
                new_ders.insert(temp);
                result_vec.push_back(mh.ders[index]);
                index_map.insert(std::make_pair(index, index - difference));
            }
            else
            {
                difference += 1;
                index_map.insert(std::make_pair(index, it->m_index));
            }
        }

        surface_patch<surface_type> *current = mh.root;
        while (current)
        {
            if (current->left)
                current = current->left;
            else
                break;
        }
        int count = 0;
        while (current)
        {
            count += 1;
            for (int i = 0; i < 4; ++i)
                current->point_index[i] = index_map[current->point_index[i]];
            surface_patch<surface_type> *temp = current->root;
            if (temp == nullptr)
                break;
            if (temp->right == current)
            {
                current = temp;
            }
            else
            {
                current = temp->right;
                while (current->left != nullptr)
                {
                    current = current->left;
                }
            }
                
        }
        mh.ders = result_vec;
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS eval_boundary_index(const surface_patch<surface_type> *patch, std::vector<int> &indexs)
    {

        const surface_patch<surface_type> *patch_root = patch->root;
        if (patch_root == nullptr)
        {
            //错误
            return ENUM_NURBS::NURBS_ERROR;
        }

        std::deque<const surface_patch<surface_type> *> patchs;
        indexs.clear();
        if (patch_root->split_direction == ENUM_DIRECTION::U_DIRECTION)
        {
            if (patch_root->left == patch)
            {
                indexs.push_back(patch->point_index[1]);
                const surface_patch<surface_type> *current = patch;
                
                while (current->left && current->split_direction == ENUM_DIRECTION::U_DIRECTION)
                {
                    current = current->right;
                }
                if (current->right != nullptr)
                    patchs.push_back(current->right);
                // else
                //     patchs.push_back(current);
                current = current->left ? current->left : current;
                
                while (true)
                {
                    if (current->left == nullptr)
                    {
                        indexs.push_back(current->point_index[2]);
                        if (patchs.empty() == false)
                        {
                            current = patchs.back();
                            patchs.pop_back();
                        }      
                        else
                            break;
                        // current = patchs.back();
                        // patchs.pop_back();
                    }
                    else
                    {
                        if (current->split_direction == ENUM_DIRECTION::U_DIRECTION)
                        {
                            current = current->right;
                        }
                        else
                        {
                            patchs.push_back(current->right);
                            current = current->left;
                        }
                    }
                }
                
            }
            else
            {
                indexs.push_back(patch->point_index[0]);
                const surface_patch<surface_type> *current = patch;
                
                while (current->left && current->split_direction == ENUM_DIRECTION::U_DIRECTION)
                {
                    current = current->left;
                }

                if (current->right != nullptr)
                    patchs.push_back(current->right);
                // else
                //     patchs.push_back(current);
                current = current->left ? current->left : current;
                
                // patchs.push_back(current->right);
                // current = current->left;
                
                while (true)
                {
                    if (current->left == nullptr)
                    {
                        indexs.push_back(current->point_index[3]);
                        
                        if (patchs.empty() == false)
                        {
                            current = patchs.back();
                            patchs.pop_back();
                        }      
                        else
                            break;
                    }
                    else
                    {
                        if (current->split_direction == ENUM_DIRECTION::U_DIRECTION)
                        {
                            current = current->left;
                        }
                        else
                        {
                            patchs.push_back(current->right);
                            current = current->left;
                        }
                    }
                }
           
            }
        }
        else
        {
            if (patch_root->left == patch)
            {
                indexs.push_back(patch->point_index[3]);
                const surface_patch<surface_type> *current = patch;
                
                while (current->left && current->split_direction == ENUM_DIRECTION::V_DIRECTION)
                {
                    current = current->right;
                }

                if (current->right != nullptr)
                    patchs.push_back(current->right);
                // else
                //     patchs.push_back(current);
                current = current->left ? current->left : current;

                // patchs.push_back(current->right);
                // current = current->left;
                
                while (true)
                {
                    if (current->left == nullptr)
                    {
                        indexs.push_back(current->point_index[2]);
                        if (patchs.empty() == false)
                        {
                            current = patchs.back();
                            patchs.pop_back();
                        }      
                        else
                            break;
                    }
                    else
                    {
                        if (current->split_direction == ENUM_DIRECTION::V_DIRECTION)
                        {
                            current = current->right;
                        }
                        else
                        {
                            patchs.push_back(current->right);
                            current = current->left;
                        }
                    }
                }
                
            }
            else
            {
                indexs.push_back(patch->point_index[0]);
                const surface_patch<surface_type> *current = patch;
                
                while (current->left && current->split_direction == ENUM_DIRECTION::V_DIRECTION)
                {
                    current = current->left;
                }
                
                if (current->right != nullptr)
                    patchs.push_back(current->right);
                // else
                //     patchs.push_back(current);
                current = current->left ? current->left : current;
                // patchs.push_back(current->right);
                // current = current->left;
                
                while (true)
                {
                    if (current->left == nullptr)
                    {
                        indexs.push_back(current->point_index[1]);
                        if (patchs.empty() == false)
                        {
                            current = patchs.back();
                            patchs.pop_back();
                        }      
                        else
                            break;
                    }
                    else
                    {
                        if (current->split_direction == ENUM_DIRECTION::V_DIRECTION)
                        {
                            current = current->left;
                        }
                        else
                        {
                            patchs.push_back(current->right);
                            current = current->left;
                        }
                    }
                }
           
            }
        }

        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS eval_mesh_from_mh(surface_mesh_helper<surface_type> &mh, mesh<2, T, dim> &surface_mesh)
    {
        // surface_mesh.m_ders = mh.ders;
        int ders_count = mh.ders.size();
        surface_mesh.m_ders.reserve(ders_count);
        for (auto &der : mh.ders)
            surface_mesh.m_ders.push_back(der);
        surface_patch<surface_type> *current = mh.root;
        if (mh.root == nullptr)
        {
            std::vector<int> indexs(mh.root->point_index.begin(), mh.root->point_index.end());
            surface_mesh.m_indexs.push_back(indexs);
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        while (current)
        {
            if (current->left)
                current = current->left;
            else
                break;
        }

        int count_loop = 0;
        while (current)
        {
            count_loop += 1;

            std::vector<int> indexs;
            surface_patch<surface_type> *temp = current->root, *temp2 = current;
            //edge u_bottom
            while (temp)
            {
                if (temp->split_direction == ENUM_DIRECTION::V_DIRECTION && temp->right == temp2)
                {
                    break;
                }
                else
                {
                    temp2 = temp;
                    temp = temp->root;
                }
            }
            if (temp == nullptr)
            {
                indexs.push_back(current->point_index[0]);
            }
            else
            {
                std::vector<int> temp_indexs;
                eval_boundary_index(temp->left, temp_indexs);
                auto begin_it = std::find(temp_indexs.begin(), temp_indexs.end(), current->point_index[0]);
                auto end_it = std::find(temp_indexs.begin(), temp_indexs.end(), current->point_index[1]);
                
                if (begin_it == temp_indexs.end() || end_it == temp_indexs.end())
                    indexs.push_back(current->point_index[0]);
                else
                    indexs.insert(indexs.end(), begin_it, end_it);
            }
            
            //edge v_right
            temp = current->root, temp2 = current;
            while (temp)
            {
                if (temp->split_direction == ENUM_DIRECTION::U_DIRECTION && temp->left == temp2)
                {
                    break;
                }
                else
                {
                    temp2 = temp;
                    temp = temp->root;
                }
            }
            if (temp == nullptr)
            {
                indexs.push_back(current->point_index[1]);
            }
            else
            {
                std::vector<int> temp_indexs;
                eval_boundary_index(temp->right, temp_indexs);
                auto begin_it = std::find(temp_indexs.begin(), temp_indexs.end(), current->point_index[1]);
                auto end_it = std::find(temp_indexs.begin(), temp_indexs.end(), current->point_index[2]);
                
                if (begin_it == temp_indexs.end() || end_it == temp_indexs.end())
                    indexs.push_back(current->point_index[1]);
                else
                    indexs.insert(indexs.end(), begin_it, end_it);
            }
            //edge u_up
            temp = current->root, temp2 = current;
            while (temp)
            {
                if (temp->split_direction == ENUM_DIRECTION::V_DIRECTION && temp->left == temp2)
                {
                    break;
                }
                else
                {
                    temp2 = temp;
                    temp = temp->root;
                }
            }
            if (temp == nullptr)
            {
                indexs.push_back(current->point_index[2]);
            }
            else
            {
                std::vector<int> temp_indexs;
                eval_boundary_index(temp->right, temp_indexs);
                std::reverse(temp_indexs.begin(), temp_indexs.end());
                auto begin_it = std::find(temp_indexs.begin(), temp_indexs.end(), current->point_index[2]);
                auto end_it = std::find(temp_indexs.begin(), temp_indexs.end(), current->point_index[3]);
                
                if (begin_it == temp_indexs.end() || end_it == temp_indexs.end())
                    indexs.push_back(current->point_index[2]);
                else if (begin_it > end_it)
                {
                    //todo
                    continue;
                }
                else
                    indexs.insert(indexs.end(), begin_it, end_it);
            }
            //edge v_left
            temp = current->root, temp2 = current;
            while (temp)
            {
                if (temp->split_direction == ENUM_DIRECTION::U_DIRECTION && temp->right == temp2)
                {
                    break;
                }
                else
                {
                    temp2 = temp;
                    temp = temp->root;
                }
            }
            if (temp == nullptr)
            {
                indexs.push_back(current->point_index[3]);
            }
            else
            {
                std::vector<int> temp_indexs;
                eval_boundary_index(temp->left, temp_indexs);
                std::reverse(temp_indexs.begin(), temp_indexs.end());
                auto begin_it = std::find(temp_indexs.begin(), temp_indexs.end(), current->point_index[3]);
                auto end_it = std::find(temp_indexs.begin(), temp_indexs.end(), current->point_index[0]);
                
                if (begin_it == temp_indexs.end() || end_it == temp_indexs.end())
                    indexs.push_back(current->point_index[3]);
                else
                    indexs.insert(indexs.end(), begin_it, end_it);
            }
            surface_mesh.m_indexs.push_back(indexs);

            //寻找下一个叶子节点
            if (current->root == nullptr)
                break;
            if (current == current->root->left)
            {
                current = current->root->right;
                while (current->left)
                {
                    current = current->left;
                }
            }
            else
            {
                surface_patch<surface_type> *temp = current->root;
                while (temp != nullptr && current == temp->right)
                {
                    current = temp;
                    temp = temp->root;
                }
                if (temp == nullptr)
                    break;
                else
                    current = temp->right;
                while (current->left)
                {
                    current = current->left;
                }
                
                
            }

           

        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename surface_type, typename T = typename surface_type::Type, int dim = surface_type::dimension>
    ENUM_NURBS mesh_help_to_mesh(surface_mesh_helper<surface_type> &mh, mesh<2, T, dim> &surface_mesh)
    {
        remove_multiple_point(mh);
        eval_mesh_from_mh(mh, surface_mesh);
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename curve_type, typename T = typename curve_type::Type, int dim = curve_type::dimension>
    ENUM_NURBS mesh_help_to_mesh(curve_mesh_helper<curve_type> &mh, mesh<1, T, dim> &curve_mesh)
    {
        curve_mesh.m_ders.clear();
        for (auto it = mh.ders.begin(); it != mh.ders.end(); ++it)
        {
            Eigen::Matrix<Eigen::Vector<T, dim>, 1, 1> vec;
            vec(0, 0) = (*it)[0];
            curve_mesh.m_ders.push_back(vec);
        }
        curve_mesh.m_indexs.clear();

        curve_patch<curve_type> *current = mh.root;
        while (current)
        {
            if (current->left == nullptr && current->right == nullptr)
            {        
                curve_patch<curve_type> *temp = current;
                current = current->root;
                std::vector<int> indexs(2);
                indexs[0] = temp->point_index[0];
                indexs[1] = temp->point_index[1];
                curve_mesh.m_indexs.push_back(indexs);
                if (current == nullptr)
                {
                    // delete temp;
                    break;
                }
                

                if (temp == current->left)
                    current = current->right;
                else
                {
                    while (current)
                    {
                        if (current->root && current->root->right == current)
                            current = current->root;
                        else if (current->root == nullptr)
                        {
                            current = nullptr;
                            break;
                        }
                        else
                        {
                            current = current->root->right;
                            break;
                        }
                    }
                    //if (current == nullptr)
                    //    break;
                    //else
                    //{
                    //    current = current->root->right;
                    //}
                }
            }
            else if (current->left != nullptr)
                current = current->left;
            else
                current = current->right;
        }
        return ENUM_NURBS::NURBS_SUCCESS;    
    }



}

