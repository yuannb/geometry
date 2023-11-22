#pragma once
#include "nurbs_surface.h"
#include "cone.h"
#include "declare.h"
#include "memory"

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
                if (current->left)
                {
                    current = current->left;
                    continue;
                }
                curve_patch<curve_type> *temp = current->root;
                
                //寻找下一个节点
                if (temp == nullptr)
                {
                    delete current;
                    break;
                }
                if (temp->left == current)
                {
                    temp->left = nullptr;
                }
                else
                {
                    temp->right = nullptr;
                }
                if (temp->right != nullptr)
                    current = temp->right;
                else
                {
                    //找到下一个节点
                    while (temp)
                    {
                        if (temp->root == nullptr)
                        {
                            current = temp;
                            temp = temp->root;
                            delete current;
                        }
                        else if (temp->root->right == temp)
                        {
                            current = temp;
                            temp = temp->root;
                            delete current;
                        }
                        else
                        {
                            current = temp;
                            temp = temp->root->right;
                            delete current;
                            break;
                        }
                    }
                    current = temp;   
                }
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
            
            
            T a = ders_vec[2].dot(dist_vec);
            T b = ders_vec[2].dot(dist_vec) + tanget_vec_square_len;
            
            
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
            std::cout << "index " << index << std::endl;
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
                if (dis > min_dis && min_dis < 0.0)
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
                        // //判断u向切向锥和c锥merge之后是否大于90度, 如果小于90度, 那么此patch就没有极值点
                        // cone<T, dim> temp(c);
                        // temp.merge_cone(current->u_cone);
                        // if (temp.m_angle < M_PI_2 - TDEFAULT_ERROR<T>::value)
                        // {
                        //     exit_minimal_point = false;
                        // }
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
                if (current->left)
                {
                    current = current->left;
                    continue;
                }
                surface_patch<surface_type> *temp = current->root;
                
                //寻找下一个节点
                if (temp == nullptr)
                {
                    delete current;
                    break;
                }
                if (temp->left == current)
                {
                    temp->left = nullptr;
                }
                else
                {
                    temp->right = nullptr;
                }
                if (temp->right != nullptr)
                    current = temp->right;
                else
                {
                    //找到下一个节点
                    while (temp)
                    {
                        if (temp->root == nullptr)
                        {
                            current = temp;
                            temp = temp->root;
                            delete current;
                        }
                        else if (temp->root->right == temp)
                        {
                            current = temp;
                            temp = temp->root;
                            delete current;
                        }
                        else
                        {
                            current = temp;
                            temp = temp->root->right;
                            delete current;
                            break;
                        }
                    }
                    current = temp;   
                }
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
        int countt = 0;
        do
        {
            std::cout << countt << std::endl;
            ++countt;
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
    ENUM_NURBS find_nearst_point_on_surface(const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &surf, const Eigen::Vector<T, dim> &point,  T &u, T &v, const Box<T, 2> &uv_box , Eigen::Vector<T, dim> &nearst_point, T eps = TDEFAULT_ERROR<T>::value)
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
    ENUM_NURBS find_nearst_point_on_surface(const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational> &surf, const std::vector<Eigen::Vector<T, dim>> &points, std::vector<T> &us, std::vector<T> &vs, std::vector<Eigen::Vector<T, dim>> &nearst_points)
    {

        //TODO: 先将曲面沿着重节点=阶数的节点切开, 防止曲面切向在此处波动过大, 以至于离散效果很差

        //先找曲面的极值点
        int count = 0;
        //1.离散
        using surface_type = nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>;
        surface_mesh_helper<nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>> mh;
        clock_t start_time = clock();
        disc_surface(&surf, mh, TDEFAULT_ERROR<double>::value, 0.5, 0.1, 1.0);
        clock_t end_time=clock();
        std::cout << "find_nearst_point_on_surface: disc_surface " <<(double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << std::endl;
        Box<T, 2> ends_knots;
        surf.get_uv_box(ends_knots);
        Eigen::Vector<T, dim> corner_point;
        int points_count = points.size();
        
        us.reserve(points_count);
        vs.reserve(points_count);
        nearst_points.reserve(points_count);
        clock_t start_time2 = clock();
        for (int index = 0; index < points_count; ++index)
        {
            // std::cout << "index " << index << std::endl;
            surface_patch<surface_type> *current = mh.root;
            Eigen::Vector<T, dim> point = points[index];
            T u, v;
            Eigen::Vector<T, dim> nearst_point;
            surf.point_on_surface(ends_knots.Min[0], ends_knots.Min[1], corner_point);
            T min_dis = (point - corner_point).norm();
            
            surf.point_on_surface(ends_knots.Min[0], ends_knots.Max[1], corner_point);
            T temp = (point - corner_point).norm();
            min_dis = std::min(min_dis, temp);
            
            surf.point_on_surface(ends_knots.Max[0], ends_knots.Min[1], corner_point);
            temp = (point - corner_point).norm();
            min_dis = std::min(min_dis, temp);
            
            surf.point_on_surface(ends_knots.Max[0], ends_knots.Max[1], corner_point);
            temp = (point - corner_point).norm();
            min_dis = std::min(min_dis, temp);

            // int num = 0;

            //2.寻找曲面极值点; 先判断是否在包围盒内部
            while (current != nullptr)
            {
                T u_min = current->uv_box.Min[0];
                T u_max = current->uv_box.Max[0];
                T v_min = current->uv_box.Min[1];
                T v_max = current->uv_box.Max[1];
                // count += 1;
                // std::cout << "count " << count << std::endl;
                bool exit_minimal_point = true;
                //先计算点与box的最近距离
                T dis = current->box.eval_minimal_distance(point);
                if (dis > min_dis && min_dis < 0.0)
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
                        // temp.merge_cone(current->v_cone);
                        // if (temp.m_angle < M_2_PI - TDEFAULT_ERROR<T>::value)
                        //     exit_minimal_point = false;
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
                    find_nearst_point_on_surface(surf, point, current_u, current_v, current->uv_box, current_point);
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
            //TODO: 3.找曲线的极值点
    
            us.push_back(u);
            vs.push_back(v);
            nearst_points.push_back(nearst_point);
        }
        clock_t end_time2 = clock();
        std::cout << "find_nearst_point_on_surface: else " <<(double)(end_time2 - start_time2) / CLOCKS_PER_SEC << "s" << std::endl;
        std::cout << "iterator_count " << count << std::endl;

        return ENUM_NURBS::NURBS_SUCCESS;
    }

}

