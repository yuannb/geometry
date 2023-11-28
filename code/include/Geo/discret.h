#pragma once
#include "nurbs_surface.h"
#include "cone.h"
#include "declare.h"
#include "memory"
#include <unordered_set>
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
    ENUM_NURBS eval_nurbs_tangnet_cone(const Eigen::VectorX<Eigen::Matrix<T, dim, Eigen::Dynamic>> &control_points,
        std::vector<cone<T, dim>> &tangent_cones)
    {
        int v_points_count = control_points.rows();
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
        std::vector<std::array<int, 2>> point_indexs; //TODO : 删除
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
            bool flag = check_patch(*current, mh.ders, dist_eps, angle_eps, chord_eps);
            
            if (flag == true)
            {
                mh.point_indexs.push_back(current->point_index); //TODO ： 待删除
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


    template<typename T, int dim, bool is_rational>
    ENUM_NURBS find_nearst_point_on_curve_inner(const nurbs_curve<T, dim, is_rational, -1, -1> &cur, const Eigen::Vector<T, dim> &point,  T &u, const Box<T, 1> &u_box , Eigen::Vector<T, dim> &nearst_point, T eps = TDEFAULT_ERROR<T>::value)
    {
        T min = u_box.Min[0];
        T max = u_box.Max[0];

        cur.point_on_curve(u, nearst_point);
        T min_length = (point - nearst_point).squaredNorm();
        Eigen::Vector<Eigen::Vector<T, dim>, 3> ders_vec;
        T distance;
        T current_u = u;
        for (int loop_index = 0; loop_index < SURFACE_ITERATE_DEEP; ++loop_index)
        {
            cur.derivative_on_curve<2>(current_u, ders_vec);
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
            bool is_closed_flag = cur.is_closed();
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


    template<typename T, int dim, bool is_rational>
    ENUM_NURBS find_nearst_point_on_curve_inner(const nurbs_curve<T, dim, is_rational, -1, -1> &cur, const std::vector<Eigen::Vector<T, dim>> &points, const std::vector<T> &min_distances,  std::vector<T> &us, std::vector<Eigen::Vector<T, dim>> &nearst_points)
    {
        //1.离散
        using curve_type = nurbs_curve<T, dim, is_rational, -1, -1>;
        curve_mesh_helper<curve_type> mh;
        clock_t start_time = clock();
        disc_curve(&cur, mh, TDEFAULT_ERROR<double>::value, 0.1, 0.1, 1.0);
        clock_t end_time=clock();
        std::cout << "find_nearst_point_on_surface: disc_surface " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;
        Box<T, 1> ends_knots;
        cur.get_u_box(ends_knots);
        Eigen::Vector<T, dim> corner_point1;
        Eigen::Vector<T, dim> corner_point2;
        cur.point_on_curve(ends_knots.Min[0], corner_point1);
        cur.point_on_curve(ends_knots.Max[0], corner_point2);
        int points_count = points.size();
        
        us.reserve(points_count);
        // nearst_points.reserve(points_count);
        clock_t start_time2 = clock();
        for (int index = 0; index < points_count; ++index)
        {
            // std::cout << "index " << index << std::endl;
            curve_patch<curve_type> *current = mh.root;
            Eigen::Vector<T, dim> point = points[index];
            T u = us[index];
            Eigen::Vector<T, dim> nearst_point = nearst_points[index];
            T min_dis = min_distances[index];
            
            //2.寻找曲线极值点; 先判断是否在包围盒内部
            while (current != nullptr)
            {
                // count += 1;
                // std::cout << count << std::endl;
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
                    // num += 1;
                    // std::cout << num << std::endl;
                    find_nearst_point_on_curve_inner(cur, point, current_u, current->u_box, current_point);
                    T current_dis = (current_point - point).norm();
                    if (min_dis > current_dis)
                    {
                        nearst_point = current_point;
                        u = current_u;
                        min_dis = current_dis;
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
            
            T d1 = (points[index] - corner_point1).norm();
            T d2 = (points[index] - corner_point2).norm();
            if (d1 < min_dis || d2 < min_dis)
            {
                if (d1 < d2)
                {
                    us[index] = ends_knots.Min[0];
                    nearst_points[index] = corner_point1;
                }
                else
                {
                    us[index] = ends_knots.Max[0];
                    nearst_points[index] = corner_point2;
                }
            }
            else
            {
                us[index] = u;
                nearst_points[index] = nearst_point;
            }
            
        }
        clock_t end_time2 = clock();
        std::cout << "find_nearst_point_on_surface: else " <<(double)(end_time2 - start_time2) / CLOCKS_PER_SEC << "s" << std::endl;
       

        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename T, int dim, bool is_rational>
    ENUM_NURBS find_nearst_point_on_curve(const nurbs_curve<T, dim, is_rational, -1, -1> &cur, const std::vector<Eigen::Vector<T, dim>> &points, std::vector<T> &us, std::vector<Eigen::Vector<T, dim>> &nearst_points)
    {

        //先将曲面沿着重节点=阶数的节点切开, 防止曲线切向在此处波动过大, 以至于离散效果很差
        std::vector<nurbs_curve<T, dim, is_rational, -1, -1>*> sperate_curves;
        cur.decompose_to_nurbs(sperate_curves);
        //将sperate_curves放入unique_ptr管理内存
        int curves_count = sperate_curves.size();
        std::vector<std::unique_ptr<nurbs_curve<T, dim, is_rational, -1, -1>>> manger_curves;
        for (int index = 0; index < curves_count; ++index)
        {
            manger_curves.push_back(std::unique_ptr<nurbs_curve<T, dim, is_rational, -1, -1>>(sperate_curves[index]));
        }

        // using curve_type = nurbs_curve<T, dim, is_rational, -1, -1>;

        Box<T, 1> ends_knots;
        cur.get_u_box(ends_knots);
        Eigen::Vector<T, dim> corner_point1;
        Eigen::Vector<T, dim> corner_point2;
        cur.point_on_curve(ends_knots.Min[0], corner_point1);
        cur.point_on_curve(ends_knots.Max[0], corner_point2);
        int points_count = points.size();
        
        us.reserve(points_count);
        nearst_points.resize(points_count);
        // clock_t start_time2 = clock();
        std::vector<T> min_distance(points_count);
        for (int index = 0; index < points_count; ++index)
        {
            Eigen::Vector<T, dim> point = points[index];
            T min_dis = (point - corner_point1).norm();
            T temp = (point - corner_point2).norm();
            if (min_dis < temp)
            {
                us.push_back(ends_knots.Min[0]);
                nearst_points[index] = corner_point1;
            }
            else
            {
                us.push_back(ends_knots.Max[0]);
                min_dis = temp;
                nearst_points[index] = corner_point2;
            }
            min_distance[index] = min_dis;
        }
        

        for (int index = 0; index < curves_count; ++index)
        {
            find_nearst_point_on_curve_inner(*sperate_curves[index], points, min_distance, us, nearst_points);
        }
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
        std::vector<std::array<int, 4>> point_indexs;
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
        int num = 0;
        // int countt = 0;
        do
        {
            // std::cout << countt << std::endl;
            // ++countt;
            std::vector<Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2>> new_ders;
            ENUM_DIRECTION dir;
            bool flag = check_patch(*curent, mh.ders, dir, new_ders, dist_eps, angle_eps, chord_eps);
            if (flag == true)
            {
                num += 1;
                mh.point_indexs.push_back(curent->point_index);
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
        std::cout << "disc " << num << std::endl;
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename T, int dim, bool is_rational>
    ENUM_NURBS find_nearst_point_on_surface_inner(const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &surf, const Eigen::Vector<T, dim> &point,  T &u, T &v, const Box<T, 2> &uv_box , Eigen::Vector<T, dim> &nearst_point, T eps = TDEFAULT_ERROR<T>::value)
    {
        T current_u = u;
        T current_v = v;
        T distance;
        T u_min = uv_box.Min[0];
        T u_max = uv_box.Max[0];
        T v_min = uv_box.Min[1];
        T v_max = uv_box.Max[1];
        surf.point_on_surface(u, v, nearst_point);
        T min_length = (nearst_point - point).squaredNorm();
        Eigen::Matrix3<Eigen::Vector<T, dim>> ders_vec;
        bool is_u_closed_flag = surf.is_u_closed();
        bool is_v_closed_flag = surf.is_v_closed();

        for (int loop_index = 0; loop_index < SURFACE_ITERATE_DEEP; ++loop_index)
        {
            surf.derivative_on_surface<2>(current_u, current_v, ders_vec);
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
            if (((next_u - current_u) * ders_vec(1, 0) + (next_v - next_v) * ders_vec(0, 1)).squaredNorm() < eps * eps)
            {
                break;
            }
            current_u = next_u;
            current_v = next_v;
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    template<typename T, int dim, bool is_rational>
    ENUM_NURBS find_nearst_point_on_surface_inner(const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &surf, const std::vector<Eigen::Vector<T, dim>> &points, std::vector<T> &min_distances, std::vector<T> &us, std::vector<T> &vs, std::vector<Eigen::Vector<T, dim>> &nearst_points)
    {
        //先找曲面的极值点
        int count = 0;
        //1.离散
        using surface_type = nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>;
        surface_mesh_helper<nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>> mh;
        clock_t start_time = clock();
        disc_surface(&surf, mh, TDEFAULT_ERROR<double>::value, 0.5, 0.1, 1.0);
        clock_t end_time=clock();
        std::cout << "find_nearst_point_on_surface: disc_surface " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;
        // Box<T, 2> ends_knots;
        // surf.get_uv_box(ends_knots);
        // Eigen::Vector<T, dim> corner_point;
        int points_count = points.size();
        
        // us.reserve(points_count);
        // vs.reserve(points_count);
        // nearst_points.reserve(points_count);
        clock_t start_time2 = clock();
        for (int index = 0; index < points_count; ++index)
        {
            // std::cout << "index " << index << std::endl;
            surface_patch<surface_type> *current = mh.root;
            Eigen::Vector<T, dim> point = points[index];
            T u = us[index], v = vs[index];
            Eigen::Vector<T, dim> nearst_point = nearst_points[index];
            T min_dis = min_distances[index];
            // surf.point_on_surface(ends_knots.Min[0], ends_knots.Min[1], corner_point);
            // T min_dis = (point - corner_point).norm();
            
            // surf.point_on_surface(ends_knots.Min[0], ends_knots.Max[1], corner_point);
            // T temp = (point - corner_point).norm();
            // min_dis = std::min(min_dis, temp);
            
            // surf.point_on_surface(ends_knots.Max[0], ends_knots.Min[1], corner_point);
            // temp = (point - corner_point).norm();
            // min_dis = std::min(min_dis, temp);
            
            // surf.point_on_surface(ends_knots.Max[0], ends_knots.Max[1], corner_point);
            // temp = (point - corner_point).norm();
            // min_dis = std::min(min_dis, temp);

            // int num = 0;
            

            //2.寻找曲面极值点; 先判断是否在包围盒内部
            while (current != nullptr)
            {
                // T u_min = current->uv_box.Min[0];
                // T u_max = current->uv_box.Max[0];
                // T v_min = current->uv_box.Min[1];
                // T v_max = current->uv_box.Max[1];
                // count += 1;
                // std::cout << "count " << count << std::endl;
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
                    count += 1;
                    find_nearst_point_on_surface_inner(surf, point, current_u, current_v, current->uv_box, current_point);
                    T current_dis = (current_point - point).norm();
                    if (min_dis > current_dis)
                    {
                        nearst_point = current_point;
                        u = current_u;
                        v = current_v;
                        min_dis = current_dis;
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
            
            // std::cout << "num " << num << std::endl;
    
            // us.push_back(u);
            // vs.push_back(v);
            // nearst_points.push_back(nearst_point);
            us[index] = u;
            vs[index] = v;
            nearst_points[index] = nearst_point;
            min_distances[index] = min_dis;
        }
        clock_t end_time2 = clock();
        std::cout << "find_nearst_point_on_surface: else " <<(double)(end_time2 - start_time2) / CLOCKS_PER_SEC << "s" << std::endl;
        std::cout << "iterator_count " << count << std::endl;

        return ENUM_NURBS::NURBS_SUCCESS;
    }

    template<typename T, int dim, bool is_rational>
    ENUM_NURBS find_nearst_point_on_surface(const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &surf, const std::vector<Eigen::Vector<T, dim>> &points, std::vector<T> &us, std::vector<T> &vs, std::vector<Eigen::Vector<T, dim>> &nearst_points)
    {
         //先将曲面沿着重节点=阶数的节点切开, 防止曲线切向在此处波动过大, 以至于离散效果很差
        Eigen::MatrixX<nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>*> sperate_surfaces;
        surf.decompose_to_nurbs(sperate_surfaces);
        
        //将sperate_curves放入unique_ptr管理内存
        int v_count = sperate_surfaces.rows();
        int u_count = sperate_surfaces.cols();


        std::vector<std::unique_ptr<nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>>> manger_surface;
        for (int v_index = 0; v_index < v_count; ++v_index)
        {
            for (int u_index = 0; u_index < u_count; ++u_index)
            {
                manger_surface.push_back(std::unique_ptr<nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>>(sperate_surfaces(v_index, u_index)));
            }
        }

        Box<T, 2> uv_box;
        surf.get_uv_box(uv_box);
        Eigen::Vector<T, dim> temp_point;
        std::vector<Eigen::Vector<T, dim>> corner_points;
        surf.point_on_surface(uv_box.Min[0], uv_box.Min[1], temp_point);
        corner_points.push_back(temp_point);
        surf.point_on_surface(uv_box.Min[0], uv_box.Max[1], temp_point);
        corner_points.push_back(temp_point);
        surf.point_on_surface(uv_box.Max[0], uv_box.Min[1], temp_point);
        corner_points.push_back(temp_point);
        surf.point_on_surface(uv_box.Max[0], uv_box.Max[1], temp_point);
        corner_points.push_back(temp_point);

        int points_count = points.size();
        
        us.reserve(points_count);
        vs.reserve(points_count);
        nearst_points.reserve(points_count);
        
        std::vector<T> min_distance;
        min_distance.reserve(points_count);
        for (int index = 0; index < points_count; ++index)
        {
            nearst_points.push_back(corner_points[0]);
            us.push_back(uv_box.Min[0]);
            vs.push_back(uv_box.Min[1]);
            min_distance.push_back((corner_points[0] - points[index]).norm());
            for (int j = 1; j < 4; ++j)
            {
                T distance = (corner_points[j] - points[index]).norm();
                if (distance < min_distance[index])
                {
                    min_distance[index] = distance;
                    nearst_points[index] = corner_points[j];
                    us[index] = j > 1 ? uv_box.Max[0] : uv_box.Min[0];
                    vs[index] = j % 2 == 0 ? uv_box.Min[1] : uv_box.Max[1];
                }
            }
        }

        //先计算曲面上的极值点
        for (int v_index = 0; v_index < v_count; ++v_index)
        {
            for (int u_index = 0; u_index < u_count; ++u_index)
            {
                find_nearst_point_on_surface_inner(*sperate_surfaces(v_index, u_index), points, min_distance, us, vs, nearst_points);
            }
        }

        //再计算曲线上的极值点
        //先将c0处的等参线取出
        std::vector<T> u_params, v_params;
        std::vector<nurbs_curve<T, dim, is_rational, -1, -1> *> u_nurbs_curves, v_nurbs_curves;
        surf.get_c0_isoparameter_curve(u_nurbs_curves, u_params, v_nurbs_curves, v_params);
        int u_curves_count = u_nurbs_curves.size();
        int v_curves_count = v_nurbs_curves.size();

        std::vector<std::unique_ptr<nurbs_curve<T, dim, is_rational, -1, -1>>> manger_curves1, manger_curves2;
        for (int index = 0; index < u_curves_count; ++index)
        {
            manger_curves1.push_back(std::unique_ptr<nurbs_curve<T, dim, is_rational, -1, -1>>(u_nurbs_curves[index]));
        }
        for (int index = 0; index < v_curves_count; ++index)
        {
            manger_curves2.push_back(std::unique_ptr<nurbs_curve<T, dim, is_rational, -1, -1>>(v_nurbs_curves[index]));
        }


        for (int index = 0; index < u_curves_count; ++index)
        {
            std::vector<T> params;
            std::vector<Eigen::Vector<T, dim>> curve_nearst_points;
            find_nearst_point_on_curve(*u_nurbs_curves[index], points, params, curve_nearst_points);
            for (int j = 0; j < points_count; ++j)
            {
                T distance = (points[j] - nearst_points[j]).norm();
                if (distance < min_distance[j])
                {
                    us[j] = u_params[index];
                    vs[j] = params[j];
                    min_distance[j] = distance;
                    nearst_points[j] = curve_nearst_points[j];
                }
            }
        }

        
        for (int index = 0; index < v_curves_count; ++index)
        {
            std::vector<T> params;
            std::vector<Eigen::Vector<T, dim>> curve_nearst_points;
            find_nearst_point_on_curve(*v_nurbs_curves[index], points, params, curve_nearst_points);
            for (int j = 0; j < points_count; ++j)
            {
                T distance = (points[j] - nearst_points[j]).norm();
                if (distance < min_distance[j])
                {
                    us[j] = params[j];
                    vs[j] = params[index];
                    min_distance[j] = distance;
                    nearst_points[j] = curve_nearst_points[j];
                }
            }
        }

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
            int a1 = current->point_index[0];
            int a2 = current->point_index[1];
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
            a1 = current->point_index[1];
            a2 = current->point_index[2];
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
            a1 = current->point_index[2];
            a2 = current->point_index[3];
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
            a1 = current->point_index[3];
            a2 = current->point_index[0];
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



}

