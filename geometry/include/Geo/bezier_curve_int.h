#pragma once
#include "bezier_curve.h"
#include "bezier_surface.h"
#include "nurbs_curve.h"
#include "nurbs_surface.h"
#include <array>
#include <deque>
#include <unordered_set>
#include <Eigen/Eigenvalues>
#include "convex_hell.h"

// #define DCSI

#ifdef DCSI
#include <chrono>
using namespace std::chrono;
#endif // DSSI
namespace tnurbs
{

    template<int variaty_count>
    struct Index
    {
        std::array<int, variaty_count> m_index;
        std::array<int, variaty_count> m_bounds;
        Index() 
        { 
            std::fill(m_index.begin(), m_index.end(), 0); 
            std::fill(m_bounds.begin(), m_bounds.end(), 0);
        }
        
        bool add_one()
        {
            m_index[variaty_count - 1] += 1;
            for (int i = variaty_count - 1; i > 0; --i)
            {
                if (m_index[i] != m_bounds[i] + 1)
                    break;
                else
                {
                    m_index[i] = 0;
                    m_index[i - 1] += 1;
                }
            }
            if (m_index[0] == m_bounds[0] + 1)
                return false;
            return true;
        }

        //保持第index个数不变，最小位加1
        bool add_one(int index)
        {
            static_assert(variaty_count > 1, "variaty_count <= 1");
            switch (index)
            {
            case 0:
                m_index[variaty_count - 1] += 1;
                for (int i = variaty_count - 1; i > 1; --i)
                {
                    if (m_index[i] < m_bounds[i] + 1)
                        break;
                    else
                    {
                        m_index[i] = 0;
                        m_index[i - 1] += 1;
                    }
                }
                if (m_index[1] >= m_bounds[1] + 1)
                    return false;
                break;
            case variaty_count - 1:
                m_index[variaty_count - 2] += 1;
                for (int i = variaty_count - 2; i > 0; --i)
                {
                    if (m_index[i] < m_bounds[i] + 1)
                        break;
                    else
                    {
                        m_index[i] = 0;
                        m_index[i - 1] += 1;
                    }
                }
                if (m_index[0] >= m_bounds[0] + 1)
                    return false;
                break;
            
            default:
                m_index[variaty_count - 1] += 1;
                for (int i = variaty_count - 1; i > 0; --i)
                {
                    if (i == index)
                    {
                        continue;
                    }
                    if (m_index[i] < m_bounds[i] + 1)
                        break;
                    else
                    {
                        m_index[i] = 0;
                        if (i - 1 == index)
                            m_index[i - 2] += 1;
                        else
                            m_index[i - 1] += 1;
                    }
                }
                if (m_index[0] >= m_bounds[0] + 1)
                    return false;
                break;
            }
            return true;
        }


        bool set_bounds(int index, int bound) 
        { 
            if (index > variaty_count - 1 || index < 0) return false; 
            m_bounds[index] = bound; return true; 
        }
        bool set_bounds(const std::array<int, (unsigned)variaty_count> &degrees)
        {
            m_bounds = degrees;
            return true;
        }

        bool reset()
        {
            std::fill(m_index.begin(), m_index.end(), 0);
            return true;
        }
    };

    template<typename T, int dim>
    struct help
    {
        T m_eps;
        Eigen::Vector<T, dim> m_point;
        help(const Eigen::Vector<T, dim> &point, T eps = TDEFAULT_ERROR<T>::value) { m_point = point; m_eps = eps; }
        ~help() { }
        const bool operator==(const help &right) const
        {
            if ((m_point - right.m_point).squaredNorm() < m_eps * m_eps)
                return true;
            return false;
        }
    };

    template<typename T, int dim>
    class hash_help 
    {
    public:
        size_t operator()(const help<T, dim>& hero)const {
            std::string temp;
            
            for (int index = 0; index < dim; ++index)
            {
                temp += (std::to_string(hero.m_point[index]) + ",");
            }
            return std::hash<std::string>()(temp);
        }
    };
        // 

    //基于bernstien基解多元多次多项式方程(此方法适合求解为有限多个的情况, 其余情况目前不支持)
    template<typename T, int equation_count, int variaty_count>
    class smspe
    {
    private:
        using curve_type = nurbs_curve<T, equation_count, false, -1, -1>;
        Eigen::Matrix<T, equation_count, Eigen::Dynamic> m_coeff;
        mutable Index<variaty_count> m_index;
        std::vector<int> m_skip;
        curve_compute<curve_type> m_curve_compute;


		std::vector<Eigen::Vector2<T>> project_point;
		std::vector<int> indices;
		std::vector<Eigen::Vector2<T>> convex_hell;

    public:
        smspe() = default;
        ~smspe() { }

        bool init(const Eigen::Matrix<T, equation_count, Eigen::Dynamic>& coeff, const std::array<int, variaty_count>& degrees)
        {
            int max_degree = 0;
            int coeff_count = 1;
            for (int degree : degrees)
            {
                coeff_count *= (degree + 1);
                max_degree = std::max(max_degree, degree);
            }
            if (coeff_count != coeff.cols())
                return false;

            m_index.reset();
            m_index.m_bounds = degrees;
            m_coeff = coeff;
            m_curve_compute.init(max_degree + 1, 2 * max_degree + 2, max_degree);
            return true;
        }

        bool compute_ipp(std::vector<Eigen::Vector<T, variaty_count>> &int_params)
        {
            if (sign_change(m_coeff) == false)
                return true;
            //init interval 
            Eigen::Vector<T, variaty_count> min, max;
            min.setConstant(0.0);
            max.setConstant(1.0);
            Box<T, variaty_count> init_box(min, max);
            std::deque<Box<T, variaty_count>> boxes { init_box };

            //int loop_count = 0;
            while (boxes.empty() == false)
            {
                m_skip.clear();
                //std::cout << "loop count: " << loop_count << std::endl;
                //++loop_count;
                Box<T, variaty_count> current_box = boxes.back();
                boxes.pop_back();
                if ((current_box.Max - current_box.Min).norm() < TDEFAULT_ERROR<T>::value)
                {
                    //逻辑需要处理
                    Eigen::Vector<T, variaty_count> param =  (current_box.Max + current_box.Min) / 2.0;
                    int_params.push_back(param);
                    continue;
                }
                //将bezier细分
                Eigen::Matrix<T, equation_count, Eigen::Dynamic> control_points = m_coeff;
                for (int index = 0; index < variaty_count; ++index)
                {
                    int bezier_control_points_count = m_index.m_bounds[index] + 1;
                    //std::cout << "index: " << index << std::endl;
                    Eigen::Matrix<T, equation_count, Eigen::Dynamic> bezier_control_points(equation_count, bezier_control_points_count);
                    Box<T, 1> sub_box(Eigen::Vector<T, 1>(current_box.Min[index]), Eigen::Vector<T, 1>(current_box.Max[index]));
                    if (sub_box.Max[0] - sub_box.Min[0] < KNOTS_EPS<T>::value)
                    {
                        sub_box.Max[0] = std::min(1.0, sub_box.Min[0] + KNOTS_EPS<T>::value * 5);
                        sub_box.Min[0] = std::max(0.0, sub_box.Min[0] - KNOTS_EPS<T>::value * 5);
                        
                        //m_skip.push_back(index);
                        //int step = 1;
                        //for (int j = index + 1; j < variaty_count; ++j)
                        //    step *= (m_index.m_bounds[j] + 1);
                        //int first_index = m_index.m_index[variaty_count - 1];

                        //int temp = m_index.m_bounds[variaty_count - 1] + 1;
                        //for (int j = variaty_count - 2; j >= 0; --j)
                        //{
                        //    first_index += m_index.m_index[j] * temp;
                        //    temp *= (m_index.m_bounds[j] + 1);
                        //}
                        //for (int i = 0; i < bezier_control_points_count; ++i)
                        //{
                        //    control_points.col(first_index + i * step).setConstant(0.0);
                        //}
                        //continue;
                    }
                    //int x = 0;
                    do
                    {
                        //std::cout << "x: " << x << std::endl;
                        //++x;
                        int step = 1;
                        for (int j = index + 1; j < variaty_count; ++j)
                            step *= (m_index.m_bounds[j] + 1);
                        int first_index = m_index.m_index[variaty_count - 1];

                        int temp = m_index.m_bounds[variaty_count - 1] + 1;
                        for (int j = variaty_count - 2; j >= 0 ; --j)
                        {
                            first_index += m_index.m_index[j] * temp;
                            temp *= (m_index.m_bounds[j] + 1);
                        }
                            
                        //构造bezier曲线
                        for (int i = 0; i < bezier_control_points_count; ++i)
                        {
                            bezier_control_points.col(i) = control_points.col(first_index + i * step);
                        }
                        nurbs_curve<T, equation_count, false, -1, -1> nurbs(bezier_control_points);
                        // std::cout << bezier_control_points << std::endl;
                        // std::cout << control_points << std::endl;
                        //Box<T, 1> sub_box(Eigen::Vector<T, 1>(current_box.Min[index]), Eigen::Vector<T, 1>(current_box.Max[index]));

                        //if (sub_box.Max[0] - sub_box.Min[0] < KNOTS_EPS<T>::value)
                        //{
                        //    Eigen::Vector<T, equation_count> point;
                        //    nurbs.point_on_curve((sub_box.Max[0] + sub_box.Min[0]) / 2.0, point);
                        //    for (int i = 0; i < bezier_control_points_count; ++i)
                        //    {
                        //        control_points.col(first_index + i * step) = point;
                        //    }
                        //}
                        //else
                        //{
                        //    nurbs.sub_divide(sub_box);
                        //    bezier_control_points = nurbs.get_control_points();
                        //    // std::cout << bezier_control_points << std::endl;
                        //    for (int i = 0; i < bezier_control_points_count; ++i)
                        //    {
                        //        control_points.col(first_index + i * step) = bezier_control_points.col(i);
                        //    }
                        //}
                        nurbs.sub_divide(sub_box);
                        bezier_control_points = nurbs.get_control_points();
                        // std::cout << bezier_control_points << std::endl;
                        for (int i = 0; i < bezier_control_points_count; ++i)
                        {
                            control_points.col(first_index + i * step) = bezier_control_points.col(i);
                        }
                        
                    } while (m_index.add_one(index));
                    m_index.reset();
                }

                Box<T, variaty_count> project_box; 
                if (interval_projected_polyhedron(control_points, project_box) == false)
                {
                    //boxes.pop_back();
                    continue;
                }
                //std::cout << control_points << std::endl;
                //control_points.setConstant(0.0);
                //cut box
                std::vector<Box<T, variaty_count>> split_boxes;
                split_box_by_ratio(project_box, current_box, split_boxes);
                for (const auto &b : split_boxes)
                {
                    Box<T, variaty_count> cbox = current_box;
                    cut_box(cbox, b);
                    boxes.push_back(cbox);
                }
            }
            
            //TODO: 去重
            std::unordered_set<help<T, variaty_count>, hash_help<T, variaty_count>> remove_mult;
            for (const Eigen::Vector<T, variaty_count> &param : int_params)
            {
                help<T, variaty_count> h(param);
                remove_mult.insert(h);
            }
            int_params.clear();
            for (auto it = remove_mult.begin(); it != remove_mult.end(); ++it)
            {
                int_params.push_back(it->m_point);
            }

            return true;
        }
    
        bool compute_ipp(std::vector<Box<T, variaty_count>>& int_params, T box_size = TDEFAULT_ERROR<T>::value)
        {
            if (sign_change(m_coeff) == false)
                return true;
            //init interval 
            Eigen::Vector<T, variaty_count> min, max;
            min.setConstant(0.0);
            max.setConstant(1.0);
            Box<T, variaty_count> init_box(min, max);
            std::deque<Box<T, variaty_count>> boxes{ init_box };

            Eigen::Matrix<T, equation_count, Eigen::Dynamic> control_points(equation_count, m_coeff.cols());// = m_coeff;
            //int loop_count = 0;
            while (boxes.empty() == false)
            {
                m_skip.clear();
                //std::cout << "loop count: " << loop_count << std::endl;
                //++loop_count;
                Box<T, variaty_count> current_box = boxes.back();
                boxes.pop_back();
                if ((current_box.Max - current_box.Min).norm() < box_size)
                {
                    //逻辑需要处理
                    int_params.push_back(current_box);
                    continue;
                }
                //将bezier细分
                control_points = m_coeff;
                for (int index = 0; index < variaty_count; ++index)
                {
                    int bezier_control_points_count = m_index.m_bounds[index] + 1;
                    Eigen::Matrix<T, equation_count, Eigen::Dynamic> bezier_control_points(equation_count, bezier_control_points_count);
                    Box<T, 1> sub_box(Eigen::Vector<T, 1>(current_box.Min[index]), Eigen::Vector<T, 1>(current_box.Max[index]));
                    if (sub_box.Max[0] - sub_box.Min[0] < KNOTS_EPS<T>::value)
                    {
                        // sub_box.Max[0] = std::min(1.0, sub_box.Min[0] + KNOTS_EPS<T>::value * 5);
                        // sub_box.Min[0] = std::max(0.0, sub_box.Min[0] - KNOTS_EPS<T>::value * 5);
                        int x = 0;
                    }
                    do
                    {
                        int step = 1;
                        for (int j = index + 1; j < variaty_count; ++j)
                            step *= (m_index.m_bounds[j] + 1);
                        int first_index = m_index.m_index[variaty_count - 1];

                        int temp = m_index.m_bounds[variaty_count - 1] + 1;
                        for (int j = variaty_count - 2; j >= 0; --j)
                        {
                            first_index += m_index.m_index[j] * temp;
                            temp *= (m_index.m_bounds[j] + 1);
                        }

                        //构造bezier曲线
                        for (int i = 0; i < bezier_control_points_count; ++i)
                        {
                            bezier_control_points.col(i) = control_points.col(first_index + i * step);
                        }

                        // nurbs_curve<T, equation_count, false, -1, -1> nurbs(bezier_control_points);
                        // nurbs.sub_divide(sub_box);
						
                        m_curve_compute.insert_bezier_knots(bezier_control_points, sub_box.Min[0], sub_box.Max[0]);
						// std::cout << "ipp: " << time1.count() << "um" << std::endl;
                        const auto& new_control_points = m_curve_compute.auxiliary_points;
                        // const auto& new_control_points = nurbs.get_control_points_ref();
                        for (int i = 0; i < bezier_control_points_count; ++i)
                        {
                            // Eigen::Vector<T, equation_count> diff = new_control_points2.col(i) - new_control_points.col(i); std::cout << "norm: " << diff.norm() << std::endl;
                            // if (diff.norm() > 1e-12)
                            // {
                            //    std::cout << new_control_points2 << std::endl;
                            //    std::cout << "11111" << std::endl;
                            //    std::cout << new_control_points << std::endl;
                            //     std::cout << "EEEEEEEEE" << std::endl;
                            // }
                            control_points.col(first_index + i * step) = new_control_points.col(i);
                        }

                    } while (m_index.add_one(index));
                    m_index.reset();
                }

			    // std::cout << "insert_bezier_knots: " << times.count() << "um" << std::endl;
                Box<T, variaty_count> project_box;
                if (interval_projected_polyhedron(control_points, project_box) == false)
                {
                    continue;
                }
                std::vector<Box<T, variaty_count>> split_boxes;
                split_box_by_ratio(project_box, current_box, split_boxes);
                for (const auto& b : split_boxes)
                {
                    Box<T, variaty_count> cbox = current_box;
                    cut_box(cbox, b);
                    boxes.push_back(cbox);
                }
            }

            return true;
        }


    private:
        Eigen::Matrix<T, equation_count, equation_count> fi_fj_mat(const Eigen::Matrix<T, equation_count, Eigen::Dynamic> &control_points)
        {
            Eigen::Matrix<T, equation_count, equation_count> result;
            for (int row = 0; row < equation_count; ++row)
            {
                result(row, row) = 0.0;
                for (int col = row + 1; col < equation_count; ++col)
                {
                    Eigen::VectorX<T> points = control_points.row(row) - control_points.row(col);
                    T dot = points.squaredNorm();
                    result(row, col) = dot;
                    result(col, row) = dot;
                }
            }
            return result;
        }
        
        bool sign_change(const Eigen::Matrix<T, equation_count, Eigen::Dynamic> &control_points) const
        {
            int count = control_points.cols();
            bool flag = false;

            for (int index = 0; index < equation_count; ++index)
            {
                for (int i = 1; i < count; ++i)
                {
                    if (control_points(index, i) * control_points(index, i - 1) <= 0)
                    {
                        flag = true;
                        break;
                    }
                }
                if (flag == false && (std::abs(control_points(index, 0)) < TDEFAULT_ERROR<T>::value || std::abs(control_points(index, count - 1)) < TDEFAULT_ERROR<T>::value))
                    flag = true;
                if (flag == false)
                    return false;
                flag = false;
            }
            return true;
        }
      
        bool interval_projected_polyhedron(const Eigen::Matrix<T, equation_count, Eigen::Dynamic> &control_points, Box<T, variaty_count> &reuslt)
        {
            bool temp = true;


            int points_count = control_points.cols();
            reuslt.Min.setConstant(0.0);
            reuslt.Max.setConstant(1.0);

			// std::vector<Eigen::Vector2<T>> project_point;
			project_point.reserve(points_count);
			// std::vector<int> indices;
			indices.reserve(points_count);
            // std::vector<Eigen::Vector2<T>> convex_hell;
            convex_hell.reserve(points_count);
            
            m_index.reset();
            Eigen::Matrix<T, variaty_count + 1, Eigen::Dynamic> current(variaty_count + 1, points_count);
            for (int f_index = 0; f_index < equation_count; ++f_index)
            {
                for (int i = 0; i < points_count; ++i)
                {
                    for (int j = 0; j < variaty_count; ++j)
                    {
                        current(j, i) = static_cast<T>(m_index.m_index[j]) / static_cast<T>(m_index.m_bounds[j]);
                    }
                    current(variaty_count, i) = control_points(f_index, i);
                    m_index.add_one();
                }
                m_index.reset();
                auto skip_it = m_skip.begin();
                for (int index = 0; index < variaty_count; ++index)
                {
                    if (skip_it != m_skip.end() && index == *skip_it)
                    {
                        continue;
                    }
                    //project(index, 2) to plane
					int current_point_count = 0;
                    for (int i = 0; i < points_count; ++i)
                    {
                        Eigen::Vector2<T> p { current(index, i) ,current(variaty_count, i) };
                        // if (i > 0)
                        // {
                        //     if ((p - project_point.back()).norm() < 1e-8)
                        //     {
                        //         continue;
                        //     }
                        // }
                        project_point.push_back(p);
                        indices.push_back(current_point_count);
                        current_point_count += 1;
                    }
                    // std::vector<Eigen::Vector2<T>> convex_hell2;
                    // graham_scan(project_point, convex_hell2);
                    graham_scan(project_point, indices, convex_hell);


                    // std::array<T, 2> interval2; 
                    // std::array<T, 2> interval3; 
                    // bool f1 = convex_hell_int_xaxis(convex_hell2, interval2);
                    // bool f2 = convex_hell_int_xaxis2(convex_hell, interval3);
                    // if (f1 != f2)
                    // {
                    //     std::cout << "error" << std::endl;
					// 	project_point.clear();
					// 	indices.clear();
					// 	convex_hell.clear();
                    // }
                    // if (f1 == true && f2 == true)
                    // {
					// 	if (interval3[0] > interval2[0] + 1e-8 || interval3[1] < interval2[1] - 1e-8)
					// 	{
					// 		std::cout << "error" << std::endl;
					// 		project_point.clear();
					// 		indices.clear();
					// 		convex_hell.clear();
					// 	}
                    // }

                    // std::array<T, 2> interval = interval3; 
                    std::array<T, 2> interval; 
                    // if (f2 == true)
                    // bool flag = convex_hell_int_xaxis(convex_hell2, interval);
                    bool flag = convex_hell_int_xaxis2(convex_hell, interval);
                    if (flag == true)
                    {
                        reuslt.Min[index] = std::max(reuslt.Min[index], interval[0]);
                        reuslt.Max[index] = std::min(reuslt.Max[index], interval[1]);
                        if (reuslt.Max[index] < reuslt.Min[index])
                        {
                            temp = false;
							project_point.clear();
							indices.clear();
							convex_hell.clear();
                            goto end_loop;
                            // return false;
                        }
                    }
                    else
                    {
                        temp = false;
						project_point.clear();
						indices.clear();
						convex_hell.clear();
                        goto end_loop;
                        // return false;
                    }
                    project_point.clear(); 
                    indices.clear();
                    convex_hell.clear();
                    int i = 0;
                }
            }
end_loop:

            return temp;
            // return true;
        }
    
        bool split_box_by_ratio(const Box<T, variaty_count> &box, const Box<T, variaty_count> &origin_box, std::vector<Box<T, variaty_count>> &split_boxes) const
        {
            split_boxes.push_back(box);
            for (int index = 0; index < variaty_count; ++index)
            {
                bool flag = box.Max[index] - box.Min[index] < 0.8;
                //if (flag == false && origin_box.Max[index] - origin_box.Min[index] > TDEFAULT_ERROR<T>::value * 10)
                if (flag == false && origin_box.Max[index] - origin_box.Min[index] > TDEFAULT_ERROR<T>::value)
                {
                    std::vector<Box<T, variaty_count>> current_split_boxes;
                    for (Box<T, variaty_count> &box : split_boxes)
                    {
                        std::array<Box<T, variaty_count>, 2> bs = box.split_at_middle(index);
                        current_split_boxes.insert(current_split_boxes.end(), bs.begin(), bs.end());
                    }
                    split_boxes = current_split_boxes;
                }
            }
            return true;
        }

        bool cut_box(Box<T, variaty_count> &box, const Box<T, variaty_count> &sub_box)
        {
            for (int index = 0; index < variaty_count; ++index)
            {
                T len = box.Max[index] - box.Min[index];
                box.Max[index] = box.Min[index] + sub_box.Max[index] * len;
                box.Min[index] += sub_box.Min[index] * len;
            }
            return true;
        }
        
        bool convex_hell_int_xaxis(const std::vector<Eigen::Vector2<T>> &convex_hell, std::array<T, 2> &intersect_interval) const
        {
            Interval<T> result;
            int count = 2;
            size_t points_count = convex_hell.size();
            size_t point_index = 0;
            //std::array<T, 2> extreams(0, 1);
            while (count > 0 && point_index < points_count)
            {
                int next_point_index = (point_index + 1) % points_count;
                if (std::abs(convex_hell[point_index][1] - 0.0) < KNOTS_EPS<T>::value)
                {
                    intersect_interval[count - 1] = convex_hell[point_index][0];
                    --count;
                }
                else if (std::abs(convex_hell[next_point_index][1] - 0.0) > KNOTS_EPS<T>::value && convex_hell[point_index][1] * convex_hell[next_point_index][1] < 0.0)
                {
                    T t = convex_hell[point_index][1] / (convex_hell[point_index][1] - convex_hell[next_point_index][1]);
                    T num = t * (convex_hell[next_point_index][0] - convex_hell[point_index][0]) + convex_hell[point_index][0];
                    //if (count == 2)
                    //{
                    //    num = std::min(std::max(convex_hell[point_index][0], convex_hell[next_point_index][0]), num + KNOTS_EPS<T>::value);
                    //}
                    //else
                    //{
                    //    num = std::max(std::min(convex_hell[point_index][0], convex_hell[next_point_index][0]), num - KNOTS_EPS<T>::value);
                    //}
                    
                    intersect_interval[count - 1] = num;
                    --count;
                }
                ++point_index;
            }
            if (count == 2)
                return false;
            else if (count == 1)
            {
                intersect_interval[0] = intersect_interval[1];
                return true;
            }

            if (intersect_interval[0] > intersect_interval[1])
            {
                std::swap(intersect_interval[0], intersect_interval[1]);
            }
            //intersect_interval[0] = std::max(0.0, intersect_interval[0] - KNOTS_EPS<T>::value);
            //intersect_interval[1] = std::min(1.0, intersect_interval[1] + KNOTS_EPS<T>::value);
                
            return true;
        }    
        
        void line_int_x_axis(const Eigen::Vector2<T>& start_point, const Eigen::Vector2<T>& end_point, Interval<T>& int_param, T tol = KNOTS_EPS<T>::value * 100) const
        {
            int is_coinciendce = 0;
            if (std::abs(start_point[1]) < tol)
            {
                int_param.m_interval[0] = std::min(start_point[0], int_param.m_interval[0]);
                int_param.m_interval[1] = std::max(start_point[0], int_param.m_interval[1]);
                is_coinciendce += 1;
            }
            if (std::abs(end_point[1]) < tol)
            {
                int_param.m_interval[0] = std::min(end_point[0], int_param.m_interval[0]);
                int_param.m_interval[1] = std::max(end_point[0], int_param.m_interval[1]);
                is_coinciendce += 1;
            }
            if (start_point[1] * end_point[1] > 0 || is_coinciendce == 2)
            {
                return;
            }

			T t = start_point[1] / (start_point[1] - end_point[1]);
			T num = t * (end_point[0] - start_point[0]) + start_point[0];
			int_param.m_interval[0] = std::min(num - KNOTS_EPS<T>::value, int_param.m_interval[0]);
			int_param.m_interval[1] = std::max(num + KNOTS_EPS<T>::value, int_param.m_interval[1]);
			int_param.m_interval[0] = std::max(0.0, int_param.m_interval[0]);
			int_param.m_interval[1] = std::min(1.0, int_param.m_interval[1]);
            return;
        }
        
        
        bool convex_hell_int_xaxis2(const std::vector<Eigen::Vector2<T>> &convex_hell, std::array<T, 2> &intersect_interval) const
        {
            Interval<T> result(1000000, -1000000);
            size_t points_count = convex_hell.size();
            size_t point_index = 0;
            while (point_index < points_count)
            {
                int next_point_index = (point_index + 1) % points_count;
                line_int_x_axis(convex_hell[point_index], convex_hell[next_point_index], result);
				++point_index;
            }
            if (result.m_interval[1] >= result.m_interval[0])
            {
                intersect_interval[0] = result.m_interval[0];
                intersect_interval[1] = result.m_interval[1];
                return true;
            }
            return false;
        }    
    };

    //TODO: is_rationl = true
    template<typename T, int dim, bool is_rational = false>
    std::vector<Eigen::Vector<T, 2>> beziers_int(const bezier_curve<T, dim, is_rational, -1> &left, const bezier_curve<T, dim, is_rational, -1> &right)
    {
        int left_degree = left.get_degree();
        int right_degree = right.get_degree();
        std::array<int, 2> degrees { left_degree, right_degree };
        Eigen::Matrix<T, dim, Eigen::Dynamic> left_control_points;
        left.get_control_points(left_control_points);
        Eigen::Matrix<T, dim, Eigen::Dynamic> right_control_points;
        right.get_control_points(right_control_points);
        Eigen::Matrix<T, dim, Eigen::Dynamic> coeff(dim, (left_degree + 1) * (right_degree + 1));
        for (int i = 0; i <= left_degree; ++i)
        {
            for (int j = 0; j <= right_degree; ++j)
            {
                int col = i * (right_degree + 1) + j;
                coeff.col(col) = left_control_points.col(i) - right_control_points.col(j);
            }
        }
        smspe<T, dim, 2> solver;
        solver.init(coeff, degrees);
        std::vector<Eigen::Vector<T, 2>> int_params;
        solver.compute_ipp(int_params);
        return int_params;
    }

    enum class curve_curve_int_type 
    {
        UNKNOWN = -1,
        NORMAL = 0,
        TANGENT = 1
    };

    template<typename T, int dim>
    struct curve_curve_int_point
    {
        Eigen::Vector2<T> m_int_param;
        Eigen::Vector<T, dim> m_int_point;
        curve_curve_int_type int_type{ curve_curve_int_type::UNKNOWN};
    };

    // 需要保证两条曲线的合法性，此类不管理曲线的内存(还未测试)
    template<typename left_curve_type, typename right_curve_type>
    class curve_curve_int
    {
        static_assert(std::is_same<typename left_curve_type::Type, typename right_curve_type::Type>::value, "left curve type is not same with right curve type");
        static_assert(left_curve_type::dimension == right_curve_type::dimension, "left curve dimension not equal right curve dimensin");
        
        using T = typename left_curve_type::Type;
        static constexpr int dim = left_curve_type::dimension;
        left_curve_type* m_left_curve;
        right_curve_type* m_right_curve;
        
        struct remove_help
        {
            T m_eps;
            curve_curve_int_point<T, dim> m_int_point;
            remove_help(const curve_curve_int_point<T, dim>& int_point, T eps = TDEFAULT_ERROR<T>::value)
            {
                m_eps = eps; 
                m_int_point = int_point; 
            }
            ~remove_help() { }
            const bool operator==(const remove_help &right) const
            {
                if ((m_int_point.m_int_point - right.m_int_point.m_int_point).squaredNorm() < m_eps * m_eps)
                    return true;
                return false;
            }
        };

        class remove_hash_help 
        {
        public:
            size_t operator()(const remove_help& hero)const {
                std::string temp;

                for (int index = 0; index < dim; ++index)
                {
                    temp += (std::to_string(hero.m_int_point.m_int_point[index]) + ",");
                }
                return std::hash<std::string>()(temp);
            }
        };

    public:
        std::vector<curve_curve_int_point<T, dim>> m_int_points;
        T m_tolerance;
        ENUM_NURBS m_error_code;


		curve_curve_int() = default;
		ENUM_NURBS reset_curves(left_curve_type* left_curve, right_curve_type* right_curve)
		{
			m_left_curve = left_curve;
			m_right_curve = right_curve;
			m_int_points.clear();
			return ENUM_NURBS::NURBS_SUCCESS;
		}
		ENUM_NURBS init(left_curve_type* left_curve, right_curve_type* right_curve, T tol = PRECISION<T>::value)
		{
			m_left_curve = left_curve;
			m_right_curve = right_curve;
			m_int_points.clear();
            m_tolerance = tol;
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		//迭代加细交点(周期性曲面待处理)(TODO: 重构)
		ENUM_NURBS intersect_point_iteration(const Box<T, 3>& domian, Eigen::Vector<T, 2> current_param, Eigen::Vector<T, 2>& intersect_param)
		{
			// TODO: closed
			Eigen::Vector<Eigen::Vector<T, dim>, 2> left_point_ders;
			Eigen::Vector<Eigen::Vector<T, dim>, 2> right_point_ders;
			Eigen::Vector<T, dim> vec;
			T min_distance = 100000;
			intersect_param = current_param;

			//迭代次数需要修改
			for (int loop_index = 0; loop_index < 2 * SURFACE_ITERATE_DEEP; ++loop_index)
			{
				m_left_curve->template derivative_on_curve<1>(current_param[0], left_point_ders);
				m_right_curve->template derivative_on_curve<1>(current_param[1], right_point_ders);
				vec = right_point_ders(0, 0) - left_point_ders[0];
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

				Eigen::Matrix<T, dim, 2> mat;
				mat.col(0) = left_point_ders[1];
				mat.col(1) = -1.0 * right_point_ders[1];

				Eigen::JacobiSVD<Eigen::Matrix<T, dim, 3>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
				Eigen::Vector<T, 2> delta = matSvd.solve(vec);
				if (matSvd.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;

				Eigen::Vector<T, 2> next_param = current_param + delta;
				for (int index = 0; index < 2; ++index)
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


		//迭代加细交点(周期性曲面待处理)(TODO: 重构)
		ENUM_NURBS intersect_point_singular_iteration(const Box<T, 3>& domian, Eigen::Vector<T, 2> current_param, Eigen::Vector<T, 2>& intersect_param)
		{
			// TODO: closed
			Eigen::Vector<Eigen::Vector<T, dim>, 3> left_point_ders;
			Eigen::Vector<Eigen::Vector<T, dim>, 3> right_point_ders;
			//Eigen::Vector<T, dim> left_point, right_point;
			Eigen::Vector<T, dim + 1> vec;
			T min_distance = 100000;
			intersect_param = current_param;

			//迭代次数需要修改
			for (int loop_index = 0; loop_index < 2 * SURFACE_ITERATE_DEEP; ++loop_index)
			{
				m_left_curve->template derivative_on_curve<2>(current_param[0], left_point_ders);
				m_right_curve->template derivative_on_curve<2>(current_param[1], right_point_ders);
				vec.template block<dim, 1>(0, 0) = right_point_ders[0] - left_point_ders[0];
                Eigen::Vector<T, dim> cross_vec = right_point_ders[1].cross(left_point_ders[1]);
				T tangent_len1 = left_point_ders[1].norm();
				T tangent_len2 = right_point_ders[1].norm();
				T coeff = -1.0 / (tangent_len1 * tangent_len2);
				vec[dim] = cross_vec.squaredNorm() * coeff;

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

				Eigen::Matrix<T, dim + 1, 2> mat;
				mat.template block<dim, 1>(0, 0) = left_point_ders[1];
                mat.template block<dim, 1>(0, 1) = -1.0 * right_point_ders[1];
				mat(dim, 0) = 2 * cross_vec.dot(right_point_ders[1].cross(left_point_ders[2])) * coeff;
				mat(dim, 2) = 2 * cross_vec.dot(right_point_ders[2].cross(left_point_ders[1])) * coeff;

				Eigen::JacobiSVD<Eigen::Matrix<T, dim + 1, 2>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
				Eigen::Vector<T, 2> delta = matSvd.solve(vec);
				if (matSvd.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;

				Eigen::Vector<T, 2> next_param = current_param + delta;
				for (int index = 0; index < 2; ++index)
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

        std::vector<curve_curve_int_point<T, dim>> bezier_curve_int_bezier_curve(const left_curve_type& left, const right_curve_type& right)
		{
			int left_degree = left.get_degree();
			int right_degree = left.get_degree();
			std::array<int, 2> degrees{ left_degree, right_degree };
			auto left_control_points = left.get_control_points();
			auto right_control_points = right.get_control_points();

			int count = (left_degree + 1) * (right_degree + 1);
			Eigen::Matrix<T, dim, Eigen::Dynamic> coeff(dim, count);
			Index<2> index;
			index.reset();
			index.set_bounds(degrees);
			for (int i = 0; i < count; ++i)
			{
				if constexpr (left_curve_type::is_ratio == false && right_curve_type::is_ratio == false)
				{
					coeff.col(i) = left_control_points.col(index.m_index[0]) - right_control_points.col(index.m_index[1]);
				}
				else if constexpr (left_curve_type::is_ratio == true && right_curve_type::is_ratio == true)
				{
					coeff.col(i) = left_control_points.template block<dim, 1>(0, index.m_index[0]) * right_control_points(dim, index.m_index[1]) -
						right_control_points.template block<dim, 1>(0, index.m_index[1]) * left_control_points(dim, index.m_index[0]);
				}
				else if constexpr (left_curve_type::is_ratio == true && right_curve_type::is_ratio == false)
				{
					coeff.col(i) = left_control_points.template block<dim, 1>(0, index.m_index[0]) -
						right_control_points.col(index.m_index[1]) * left_control_points(dim, index.m_index[0]);
				}
				else
				{
					coeff.col(i) = left_control_points.col(index.m_index[0]) * right_control_points(dim, index.m_index[1]) -
						right_control_points.template block<dim, 1>(0, index.m_index[1]);
				}
				index.add_one();
			}

			smspe<T, dim, 2> solver;
			solver.init(coeff, degrees);
			std::vector<Box<T, 3>> int_boxes;
			solver.compute_ipp(int_boxes);
            std::vector<curve_curve_int_point<T, dim>> int_points;
			Box<T, 2> domian;
            std::array<T, 2> curve_ends_knots;
            left.get_ends_knots(curve_ends_knots);
			domian.Min[0] = curve_ends_knots[0];
			domian.Max[0] = curve_ends_knots[1];
            
            right.get_ends_knots(curve_ends_knots);
			domian.Min[1] = curve_ends_knots[0];
			domian.Max[1] = curve_ends_knots[1];
            Eigen::Vector<T, 2> result_param;
			for (Box<T, 2>&int_box : int_boxes)
			{
				bool flag = false;
				int_box.Min[0] = domian.Min[0] + (domian.Max[0] - domian.Min[0]) * int_box.Min[0];
				int_box.Max[0] = domian.Min[0] + (domian.Max[0] - domian.Min[0]) * int_box.Max[0];
				int_box.Min[1] = domian.Min[1] + (domian.Max[1] - domian.Min[1]) * int_box.Min[1];
				int_box.Max[1] = domian.Min[1] + (domian.Max[1] - domian.Min[1]) * int_box.Max[1];

                Eigen::Vector<T, 2> initial_param = int_box.get_middle_point();
                Eigen::Vector<T, 2> param;
				if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration(domian, initial_param, param))
				{
					if (int_box.is_contain_point(param, TDEFAULT_ERROR<T>::value))
					{
						flag = true;
						result_param = param;
						initial_param = param;
					}
				}
                bool is_tangent = false;
				if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_singular_iteration(domian, initial_param, param))
				{
					if (int_box.is_contain_point(param, TDEFAULT_ERROR<T>::value))
					{
						flag = true;
                        is_tangent = true;
						result_param = param;
					}
				}
				if (flag == true)
				{
                    int_points.push_back(curve_curve_int_point<T, dim>{});
                    curve_curve_int_point<T, dim>& end = int_points.back();
                    end.int_type = is_tangent == true ? curve_curve_int_type::TANGENT : curve_curve_int_type::NORMAL;
                    end.m_int_param = result_param;
                    m_left_curve->point_on_curve(result_param[0], end.m_int_point);
				}
			}

			//TODO: 去重
			// std::unordered_set<help<T, 3>, hash_help<T, 3>> remove_mult;
			// for (const Eigen::Vector<T, 3>&param : iter_int_params)
			// {
			// 	help<T, 3> h(param);
			// 	remove_mult.insert(h);
			// }
			// iter_int_params.clear();
			// for (auto it = remove_mult.begin(); it != remove_mult.end(); ++it)
			// {
			// 	iter_int_params.push_back(it->m_point);
			// }

			return int_points;
		}

        void get_curves_int_pair(std::vector<std::pair<left_curve_type, right_curve_type>>& int_curves_pairs)
        {
            // TOOD: 分层碰撞; 内存管理

            std::vector<left_curve_type *> left_sub_curves;
            m_left_curve->decompose_to_bezier(left_sub_curves);

            std::vector<left_curve_type *> right_sub_curves;
            m_right_curve->decompose_to_bezier(right_sub_curves);
            
            size_t left_curves_count = left_sub_curves.size();
            size_t right_curves_count = right_sub_curves.size();
            std::vector<Box<T, dim>> left_curve_boxes;
            std::vector<Box<T, dim>> right_curve_boxes;
            for (size_t curve_index = 0; curve_index < left_curves_count; ++curve_index)
            {
                left_curve_boxes.push_back(Box<T, dim>());
                left_sub_curves[curve_index]->get_box(left_curve_boxes.back());
            }
            for (size_t curve_index = 0; curve_index < right_curves_count; ++curve_index)
            {
                right_curve_boxes.push_back(Box<T, dim>());
                right_sub_curves[curve_index]->get_box(right_curve_boxes.back());
            }
            Box<T, dim> int_box;
            for (size_t left_curve_index = 0; left_curve_index < left_curves_count; ++left_curve_index)
            {
                for (size_t right_curve_index = 0; right_curve_index < right_curves_count; ++right_curve_index)
                {
					if (left_curve_boxes[left_curve_index].intersect(right_curve_boxes[right_curve_index], int_box) == true)
					{
						int_curves_pairs.push_back(std::make_pair(left_curve_type(*left_sub_curves[left_curve_index]), right_curve_type(*right_sub_curves[right_curve_index])));
					}
                }
            }
            for (size_t curve_index = 0; curve_index < left_curves_count; ++curve_index)
            {
                delete left_sub_curves[curve_index];
            }
            for (size_t curve_index = 0; curve_index < right_curves_count; ++curve_index)
            {
                delete right_sub_curves[curve_index];
            }
            return;
        }
        void run_intersect()
        {
            std::vector<std::pair<left_curve_type, right_curve_type>> curves_pairs;
            get_curves_int_pair(curves_pairs);
           
			std::unordered_set<remove_help, remove_hash_help> remove_mult;
            for (auto& curve_pair : curves_pairs)
            {
                std::vector<curve_curve_int_point<T, dim>> int_points = bezier_curve_int_bezier_curve(curves_pairs.first, curves_pairs.second); 
                for (curve_curve_int_point<T, dim>& int_point : int_points)
                {
                    remove_help h(int_point);
                    remove_mult.insert(h);
                }
            }
            // m_int_points.insert(m_int_points.end(), remove_mult.begin(), remove_mult.end()); 
            for (auto it = remove_mult.begin(); it != remove_mult.end(); ++it)
            {
                m_int_points.push_back(it->m_int_point);
            }
            return;
        }
    };

    //TODO: is_rationl = true
    template<typename T, int dim, bool is_rational = false>
    std::vector<Eigen::Vector<T, 3>> bezier_curve_int_bezier_surface(const bezier_curve<T, dim, is_rational, -1>& left, const bezier_surface<T, dim, -1, -1, is_rational>& right)
    {
        int left_degree = left.get_degree();
        Eigen::MatrixX<Eigen::Vector<T, is_rational ? dim + 1 : dim>> right_control_points;
        right.get_control_points(right_control_points);
        int v_right_degree = right_control_points.cols() - 1;
        int u_right_degree = right_control_points.rows() - 1;
        std::array<int, 3> degrees{ left_degree, u_right_degree, v_right_degree };
        Eigen::Matrix<T, dim, Eigen::Dynamic> left_control_points;
        left.get_control_points(left_control_points);

        int count = (left_degree + 1) * (u_right_degree + 1) * (v_right_degree + 1);
        Eigen::Matrix<T, dim, Eigen::Dynamic> coeff(dim, count);
        Index<3> index;
        index.reset();
        index.set_bounds(degrees);
        for (int i = 0; i < count; ++i)
        {
            //std::cout << "i: " << i << std::endl;
            coeff.col(i) = left_control_points.col(index.m_index[0]) - right_control_points(index.m_index[1], index.m_index[2]);
            index.add_one();
        }
        
        smspe<T, dim, 3> solver;
        solver.init(coeff, degrees);
        std::vector<Eigen::Vector<T, 3>> int_params;
        solver.compute_ipp(int_params);
        return int_params;
    }


    //迭代加细交点(周期性曲面待处理)(TODO: 重构)
    template<typename T, int dim, bool is_rational = false>
    ENUM_NURBS intersect_point_iteration(const nurbs_curve<T, dim, is_rational, -1, -1>& left, const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>& right,
                                        const Box<T, 3>& domian, Eigen::Vector<T, 3> current_param, Eigen::Vector<T, 3>& intersect_param)
    {
        // TODO: closed
        Eigen::Vector<Eigen::Vector<T, dim>, 2> left_point_ders;
        Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_point_ders;
        //Eigen::Vector<T, dim> left_point, right_point;
        Eigen::Vector<T, dim> vec;
        T min_distance = 100000;
        intersect_param = current_param;

        //迭代次数需要修改
        for (int loop_index = 0; loop_index < 2 * SURFACE_ITERATE_DEEP; ++loop_index)
        {
            left.template derivative_on_curve<1>(current_param[0], left_point_ders);
            right.template derivative_on_surface<1>(current_param[1], current_param[2], right_point_ders);
            vec = right_point_ders(0, 0) - left_point_ders[0];
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

            Eigen::Matrix<T, dim, 3> mat;
            mat.col(0) = left_point_ders[1];
            mat.col(1) = -1.0 * right_point_ders(1, 0);
            mat.col(2) = -1.0 * right_point_ders(0, 1);

            Eigen::JacobiSVD<Eigen::Matrix<T, dim, 3>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
            Eigen::Vector<T, 3> delta = matSvd.solve(vec);
            if (matSvd.info() != Eigen::Success)
                return ENUM_NURBS::NURBS_ERROR;

            Eigen::Vector<T, 3> next_param = current_param + delta;
            for (int index = 0; index < 3; ++index)
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


    //迭代加细交点(周期性曲面待处理)(TODO: 重构)
    template<typename T, int dim, bool is_rational = false>
    ENUM_NURBS intersect_point_singular_iteration(const nurbs_curve<T, dim, is_rational, -1, -1>& left, const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>& right,
        const Box<T, 3>& domian, Eigen::Vector<T, 3> current_param, Eigen::Vector<T, 3>& intersect_param)
    {
        // TODO: closed
        Eigen::Vector<Eigen::Vector<T, dim>, 3> left_point_ders;
        Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_point_ders;
        //Eigen::Vector<T, dim> left_point, right_point;
        Eigen::Vector<T, dim + 1> vec;
        T min_distance = 100000;
        intersect_param = current_param;

        //迭代次数需要修改
        for (int loop_index = 0; loop_index < 2 * SURFACE_ITERATE_DEEP; ++loop_index)
        {
            left.template derivative_on_curve<2>(current_param[0], left_point_ders);
            right.template derivative_on_surface<2>(current_param[1], current_param[2], right_point_ders);
            vec.template block<dim, 1>(0, 0) = right_point_ders(0, 0) - left_point_ders[0];
            Eigen::Vector<T, dim> surface_normal = right_point_ders(1, 0).cross(right_point_ders(0, 1));
            T normal_len = surface_normal.norm();
            T tangent_len = left_point_ders[1].norm();
            T coeff = -1.0 / (normal_len * tangent_len);
            vec[dim] = -surface_normal.dot(left_point_ders[1]) * coeff;

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

            Eigen::Matrix<T, dim + 1, 3> mat;
            mat.template block<dim, 1>(0, 0) = left_point_ders[1];
            mat.template block<dim, 1>(0, 1) = -1.0 * right_point_ders(1, 0);
            mat.template block<dim, 1>(0, 2) = -1.0 * right_point_ders(0, 1);
            mat(dim, 0) = right_point_ders(1, 0).cross(right_point_ders(0, 1)).dot(left_point_ders[2]) * coeff;
            mat(dim, 1) = (right_point_ders(2, 0).cross(right_point_ders(0, 1)) + right_point_ders(1, 0).cross(right_point_ders(1, 1))).dot(left_point_ders[1]) * coeff;
            mat(dim, 2) = (right_point_ders(1, 0).cross(right_point_ders(0, 2)) + right_point_ders(1, 1).cross(right_point_ders(0, 1))).dot(left_point_ders[1]) * coeff;

            Eigen::JacobiSVD<Eigen::Matrix<T, dim + 1, 3>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
            Eigen::Vector<T, 3> delta = matSvd.solve(vec);
            if (matSvd.info() != Eigen::Success)
                return ENUM_NURBS::NURBS_ERROR;

            Eigen::Vector<T, 3> next_param = current_param + delta;
            for (int index = 0; index < 3; ++index)
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



    template<typename T, int dim, bool is_rational = false>
    std::vector<Eigen::Vector<T, 3>> bezier_curve_int_bezier_surface(const nurbs_curve<T, dim, is_rational, -1, -1>& left, const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>& right)
    {
        int left_degree = left.get_degree();
        auto right_control_points = right.get_control_points();
        int u_right_degree = right_control_points[0].cols() - 1;
        int v_right_degree = right_control_points.size() - 1;
        std::array<int, 3> degrees{ left_degree, u_right_degree, v_right_degree };
        auto left_control_points = left.get_control_points();
         

        int count = (left_degree + 1) * (u_right_degree + 1) * (v_right_degree + 1);
        Eigen::Matrix<T, dim, Eigen::Dynamic> coeff(dim, count);
        Index<3> index;
        index.reset();
        index.set_bounds(degrees);
        for (int i = 0; i < count; ++i)
        {
            if constexpr (is_rational == false)
            {
                coeff.col(i) = left_control_points.col(index.m_index[0]) - right_control_points[index.m_index[2]].col(index.m_index[1]);
            }
            else
            {
                coeff.col(i) = left_control_points.template block<dim, 1>(0, index.m_index[0]) * right_control_points[index.m_index[2]](dim, index.m_index[1]) -
                    right_control_points[index.m_index[2]].template block<dim, 1>(0, index.m_index[1]) * left_control_points(dim, index.m_index[0]);
            }
            index.add_one();
        }

        smspe<T, dim, 3> solver;
        solver.init(coeff, degrees);
        std::vector<Box<T, 3>> int_boxes;
        solver.compute_ipp(int_boxes);
        std::vector<Eigen::Vector<T, 3>> iter_int_params;
        Box<T, 3> domian;
        domian.Min.setConstant(0.0);
        domian.Max.setConstant(1.0);
        Eigen::Vector<T, 3> result_param;
        for (Box<T, 3>& int_box : int_boxes)
        {
            bool flag = false;
            Eigen::Vector<T, 3> initial_param = int_box.get_middle_point();
            Eigen::Vector<T, 3> param;
            if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration<T, dim, is_rational>(left, right, domian, initial_param, param))
            {
                if (int_box.is_contain_point(param, TDEFAULT_ERROR<T>::value))
                {
                    flag = true;
                    result_param = param;
                    initial_param = param;
                }
            }
            if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_singular_iteration<T, dim, is_rational>(left, right, domian, initial_param, param))
            {
                if (int_box.is_contain_point(param, TDEFAULT_ERROR<T>::value))
                {
                    flag = true;
                    result_param = param;
                }
            }
            if (flag == true)
            {
                iter_int_params.push_back(result_param);
            }
        }

        //TODO: 去重
        std::unordered_set<help<T, 3>, hash_help<T, 3>> remove_mult;
        for (const Eigen::Vector<T, 3>& param : iter_int_params)
        {
            help<T, 3> h(param);
            remove_mult.insert(h);
        }
        iter_int_params.clear();
        for (auto it = remove_mult.begin(); it != remove_mult.end(); ++it)
        {
            iter_int_params.push_back(it->m_point);
        }

        return iter_int_params;
    }
    enum class curve_surface_int_type 
    {
        UNKNOWN = -1,
        NORMAL = 0,
        TANGENT = 1
    };


    template<typename T, int dim>
    struct curve_surface_int_point
    {
        Eigen::Vector3<T> m_int_param;
        Eigen::Vector<T, dim> m_int_point;
        curve_surface_int_type int_type{curve_surface_int_type::UNKNOWN};
    };

    template<typename curve_type, typename surface_type>
    class curve_surface_int
    {
        using T = typename curve_type::Type;
        static constexpr int dim = curve_type::dimension;

        static_assert(std::is_same<typename curve_type::Type, typename surface_type::Type>::value, "curve type is not same with surface type");
        static_assert(curve_type::dimension == surface_type::dimension, "curve dimension not equal surface dimensin");

        curve_type* m_curve;
        surface_type* m_surface;
        struct remove_help
        {
            T m_eps;
            curve_surface_int_point<T, dim> m_int_point;
            remove_help(const curve_surface_int_point<T, dim>& int_point, T eps = TDEFAULT_ERROR<T>::value)
            {
                m_eps = eps; 
                m_int_point = int_point; 
            }
            ~remove_help() { }
            const bool operator==(const remove_help &right) const
            {
                if ((m_int_point.m_int_point - right.m_int_point.m_int_point).squaredNorm() < m_eps * m_eps)
                    return true;
                return false;
            }
        };

        class remove_hash_help 
        {
        public:
            size_t operator()(const remove_help& hero)const {
                std::string temp;

                for (int index = 0; index < dim; ++index)
                {
                    temp += (std::to_string(hero.m_int_point.m_int_point[index]).substr(0, 3) + ",");
                }
                return std::hash<std::string>()(temp);
            }
        };
        // 


        
    public: 
        std::vector<curve_surface_int_point<T, dim>> m_int_points;
        T m_tolerance;
        ENUM_NURBS m_error_code;
    public:
        curve_surface_int() = default;
        ENUM_NURBS reset_curve(curve_type* curve)
        {
            m_curve = curve;
            m_int_points.clear();
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        ENUM_NURBS init(const curve_type* curve, const surface_type* surface, T tol = PRECISION<T>::value)
        {
            m_curve = new curve_type(*curve);
            m_surface = new surface_type(*surface);
            m_tolerance = tol;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
		//迭代加细交点(周期性曲面待处理)(TODO: 重构)
		ENUM_NURBS intersect_point_iteration(const Box<T, 3>& domian, Eigen::Vector<T, 3> current_param, Eigen::Vector<T, 3>& intersect_param)
		{
			// TODO: closed
			Eigen::Vector<Eigen::Vector<T, dim>, 2> left_point_ders;
			Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_point_ders;
			Eigen::Vector<T, dim> vec;
			T min_distance = 100000;
			intersect_param = current_param;

			//迭代次数需要修改
			for (int loop_index = 0; loop_index < 2 * SURFACE_ITERATE_DEEP; ++loop_index)
			{
				m_curve->template derivative_on_curve<1>(current_param[0], left_point_ders);
				m_surface->template derivative_on_surface<1>(current_param[1], current_param[2], right_point_ders);
				vec = right_point_ders(0, 0) - left_point_ders[0];
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

				Eigen::Matrix<T, dim, 3> mat;
				mat.col(0) = left_point_ders[1];
				mat.col(1) = -1.0 * right_point_ders(1, 0);
				mat.col(2) = -1.0 * right_point_ders(0, 1);

				Eigen::JacobiSVD<Eigen::Matrix<T, dim, 3>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
				Eigen::Vector<T, 3> delta = matSvd.solve(vec);
				if (matSvd.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;

				Eigen::Vector<T, 3> next_param = current_param + delta;
				for (int index = 0; index < 3; ++index)
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


		//迭代加细交点(周期性曲面待处理)(TODO: 重构)
		ENUM_NURBS intersect_point_singular_iteration(const Box<T, 3>& domian, Eigen::Vector<T, 3> current_param, Eigen::Vector<T, 3>& intersect_param)
		{
			// TODO: closed
			Eigen::Vector<Eigen::Vector<T, dim>, 3> left_point_ders;
			Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_point_ders;
			//Eigen::Vector<T, dim> left_point, right_point;
			Eigen::Vector<T, dim + 1> vec;
			T min_distance = 100000;
			intersect_param = current_param;
            int decrease_slow = 0;

			//迭代次数需要修改
			for (int loop_index = 0; loop_index < 2 * SURFACE_ITERATE_DEEP; ++loop_index)
			{
				m_curve->template derivative_on_curve<2>(current_param[0], left_point_ders);
				m_surface->template derivative_on_surface<2>(current_param[1], current_param[2], right_point_ders);
				vec.template block<dim, 1>(0, 0) = right_point_ders(0, 0) - left_point_ders[0];
				Eigen::Vector<T, dim> surface_normal = right_point_ders(1, 0).cross(right_point_ders(0, 1));
				T normal_len = surface_normal.norm();
				T tangent_len = left_point_ders[1].norm();
				T coeff = -1.0 / (normal_len * tangent_len);
				vec[dim] = -surface_normal.dot(left_point_ders[1]) * coeff;

				T distance = vec.squaredNorm();
				if (distance < min_distance)
				{
                    if (distance > 0.9 * min_distance)
                    {
                        decrease_slow += 1;
                    }
                    else
                    {
                        decrease_slow = 0;
                    }
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

                    if (decrease_slow > 2)
                    {
                        return ENUM_NURBS::NURBS_ERROR;
                    }

				}

				Eigen::Matrix<T, dim + 1, 3> mat;
				mat.template block<dim, 1>(0, 0) = left_point_ders[1];
				mat.template block<dim, 1>(0, 1) = -1.0 * right_point_ders(1, 0);
				mat.template block<dim, 1>(0, 2) = -1.0 * right_point_ders(0, 1);
				mat(dim, 0) = right_point_ders(1, 0).cross(right_point_ders(0, 1)).dot(left_point_ders[2]) * coeff;
				mat(dim, 1) = (right_point_ders(2, 0).cross(right_point_ders(0, 1)) + right_point_ders(1, 0).cross(right_point_ders(1, 1))).dot(left_point_ders[1]) * coeff;
				mat(dim, 2) = (right_point_ders(1, 0).cross(right_point_ders(0, 2)) + right_point_ders(1, 1).cross(right_point_ders(0, 1))).dot(left_point_ders[1]) * coeff;

				Eigen::JacobiSVD<Eigen::Matrix<T, dim + 1, 3>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
				Eigen::Vector<T, 3> delta = matSvd.solve(vec);
				if (matSvd.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;

				Eigen::Vector<T, 3> next_param = current_param + delta;
				for (int index = 0; index < 3; ++index)
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

        std::vector<curve_surface_int_point<T, dim>> bezier_curve_int_bezier_surface(const curve_type& left, const surface_type& right)
		{
			int left_degree = left.get_degree();
			auto right_control_points = right.get_control_points();
			int u_right_degree = right_control_points[0].cols() - 1;
			int v_right_degree = right_control_points.size() - 1;
			std::array<int, 3> degrees{ left_degree, u_right_degree, v_right_degree };
			auto left_control_points = left.get_control_points();

			int count = (left_degree + 1) * (u_right_degree + 1) * (v_right_degree + 1);
			Eigen::Matrix<T, dim, Eigen::Dynamic> coeff(dim, count);
			Index<3> index;
			index.reset();
			index.set_bounds(degrees);
			for (int i = 0; i < count; ++i)
			{
				if constexpr (curve_type::is_ratio == false && surface_type::is_ratio == false)
				{
					coeff.col(i) = left_control_points.col(index.m_index[0]) - right_control_points[index.m_index[2]].col(index.m_index[1]);
				}
				else if constexpr (curve_type::is_ratio == true && surface_type::is_ratio == true)
				{
					coeff.col(i) = left_control_points.template block<dim, 1>(0, index.m_index[0]) * right_control_points[index.m_index[2]](dim, index.m_index[1]) -
						right_control_points[index.m_index[2]].template block<dim, 1>(0, index.m_index[1]) * left_control_points(dim, index.m_index[0]);
				}
				else if constexpr (curve_type::is_ratio == true && surface_type::is_ratio == false)
				{
					coeff.col(i) = left_control_points.template block<dim, 1>(0, index.m_index[0]) -
						right_control_points[index.m_index[2]].col(index.m_index[1]) * left_control_points(dim, index.m_index[0]);
				}
				else
				{
					coeff.col(i) = left_control_points.col(index.m_index[0]) * right_control_points[index.m_index[2]](dim, index.m_index[1]) -
						right_control_points[index.m_index[2]].template block<dim, 1>(0, index.m_index[1]);
				}
				index.add_one();
			}
#ifdef DCSI
            auto s1 = steady_clock::now();
#endif // DSSI
			smspe<T, dim, 3> solver;
			solver.init(coeff, degrees);
			std::vector<Box<T, 3>> int_boxes;
			solver.compute_ipp(int_boxes);
#ifdef DCSI
            auto e1 = steady_clock::now();
            auto time1 = duration_cast<microseconds>(e1 - s1);
            std::cout << "ipp: " << time1.count() << "um" << std::endl;
#endif // DSSI
            std::vector<curve_surface_int_point<T, dim>> int_points;
			Box<T, 3> domian;
            std::array<T, 2> curve_ends_knots;
            left.get_ends_knots(curve_ends_knots);
			domian.Min[0] = curve_ends_knots[0];
			domian.Max[0] = curve_ends_knots[1];
			std::array<T, 2> surface_u_ends_knots, surface_v_ends_knots;
            right.get_uv_knots_end(surface_u_ends_knots, surface_v_ends_knots);
            domian.Min[1] = surface_u_ends_knots[0];
            domian.Max[1] = surface_u_ends_knots[1];
            domian.Min[2] = surface_v_ends_knots[0];
            domian.Max[2] = surface_v_ends_knots[1];
            Eigen::Vector<T, 3> result_param;
#ifdef DCSI
            auto s2 = steady_clock::now();
#endif // DSSI
			for (Box<T, 3>&int_box : int_boxes)
			{
				bool flag = false;
				int_box.Min[0] = domian.Min[0] + (domian.Max[0] - domian.Min[0]) * int_box.Min[0];
				int_box.Max[0] = domian.Min[0] + (domian.Max[0] - domian.Min[0]) * int_box.Max[0];
				int_box.Min[1] = domian.Min[1] + (domian.Max[1] - domian.Min[1]) * int_box.Min[1];
				int_box.Max[1] = domian.Min[1] + (domian.Max[1] - domian.Min[1]) * int_box.Max[1];
				int_box.Min[2] = domian.Min[2] + (domian.Max[2] - domian.Min[2]) * int_box.Min[2];
				int_box.Max[2] = domian.Min[2] + (domian.Max[2] - domian.Min[2]) * int_box.Max[2];

                Eigen::Vector<T, 3> initial_param = int_box.get_middle_point();
                Eigen::Vector<T, 3> param;
				if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration(domian, initial_param, param))
				{
					if (int_box.is_contain_point(param, TDEFAULT_ERROR<T>::value))
					{
						flag = true;
						result_param = param;
						initial_param = param;
					}
				}
                bool is_tangent = false;

				Eigen::Vector<Eigen::Vector<T, dim>, 2> left_point_ders;
				Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_point_ders;
                m_curve->template derivative_on_curve<1>(initial_param[0], left_point_ders);
                m_surface->template derivative_on_surface<1>(initial_param[1], initial_param[2], right_point_ders);
                left_point_ders[1].normalize();
                right_point_ders(1, 0).normalize();
                right_point_ders(0, 1).normalize();
                Eigen::Vector<T, dim> normal = right_point_ders(1, 0).cross(right_point_ders(0, 1));
                if (std::abs(normal.dot(left_point_ders[1])) < 1e-3)
                {
					if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_singular_iteration(domian, initial_param, param))
					{
						if (int_box.is_contain_point(param, TDEFAULT_ERROR<T>::value))
						{
							flag = true;
							is_tangent = true;
							result_param = param;
						}
					}
                }
 
                
				if (flag == true)
				{
                    int_points.push_back(curve_surface_int_point<T, dim>{});
                    curve_surface_int_point<T, dim>& end = int_points.back();
                    end.int_type = is_tangent == true ? curve_surface_int_type::TANGENT : curve_surface_int_type::NORMAL;
                    end.m_int_param = result_param;
                    m_curve->point_on_curve(result_param[0], end.m_int_point);
				}
			}

#ifdef DCSI
            auto e2 = steady_clock::now();
            auto time2 = duration_cast<microseconds>(e2 - s2);
            std::cout << "iteration: " << time2.count() << "um" << std::endl;
#endif // DSSI
			//TODO: 去重
			// std::unordered_set<help<T, 3>, hash_help<T, 3>> remove_mult;
			// for (const Eigen::Vector<T, 3>&param : iter_int_params)
			// {
			// 	help<T, 3> h(param);
			// 	remove_mult.insert(h);
			// }
			// iter_int_params.clear();
			// for (auto it = remove_mult.begin(); it != remove_mult.end(); ++it)
			// {
			// 	iter_int_params.push_back(it->m_point);
			// }

			return int_points;
		}

        void get_curve_surface_int_pair(std::vector<std::pair<curve_type, surface_type>>& int_curve_surface_pairs)
        {
            // TOOD: 分层碰撞; 内存管理

            std::vector<curve_type *> sub_curves;
            m_curve->decompose_to_bezier(sub_curves);

            Eigen::MatrixX<surface_type *> sub_surfaces;
            m_surface->decompose_to_bezier(sub_surfaces);
            size_t curves_count = sub_curves.size();
            Eigen::Index rows = sub_surfaces.rows();
            Eigen::Index cols = sub_surfaces.cols();
            std::vector<Box<T, dim>> curve_boxes;
            Eigen::MatrixX<Box<T, dim>> surface_boxes;
            surface_boxes.resize(rows, cols);
            for (size_t curve_index = 0; curve_index < curves_count; ++curve_index)
            {
                curve_boxes.push_back(Box<T, dim>());
                sub_curves[curve_index]->get_box(curve_boxes.back());
            }
            for (Eigen::Index col_index = 0; col_index < cols; ++col_index)
            {
                for (Eigen::Index row_index = 0; row_index < rows; ++row_index)
                {
                    sub_surfaces(row_index, col_index)->get_box(surface_boxes(row_index, col_index));
                }
            }
            Box<T, dim> int_box;
            for (size_t curve_index = 0; curve_index < curves_count; ++curve_index)
            {
                for (Eigen::Index col_index = 0; col_index < cols; ++col_index)
                {
                    for (Eigen::Index row_index = 0; row_index < rows; ++row_index)
                    {
                        if (curve_boxes[curve_index].intersect(surface_boxes(row_index, col_index), int_box) == true)
                        {
                            int_curve_surface_pairs.push_back(std::make_pair(curve_type(*sub_curves[curve_index]), surface_type(*sub_surfaces(row_index, col_index))));
                        }
                    }
                }
            }
            for (size_t curve_index = 0; curve_index < curves_count; ++curve_index)
            {
                delete sub_curves[curve_index];
            }
            for (Eigen::Index col_index = 0; col_index < cols; ++col_index)
            {
                for (Eigen::Index row_index = 0; row_index < rows; ++row_index)
                {
                    delete sub_surfaces(row_index, col_index);
                }
            }
            return;
        }
            
        void run_intersect()
        {
            std::vector<std::pair<curve_type, surface_type>> curve_surface_pairs;
            get_curve_surface_int_pair(curve_surface_pairs);
           
			std::unordered_set<remove_help, remove_hash_help> remove_mult;
            for (auto& curve_surface_pair : curve_surface_pairs)
            {
                std::vector<curve_surface_int_point<T, dim>> int_points = bezier_curve_int_bezier_surface(curve_surface_pair.first, curve_surface_pair.second); 
                for (curve_surface_int_point<T, dim>& int_point : int_points)
                {
                    remove_help h(int_point);
                    remove_mult.insert(h);
                }
            }
            // m_int_points.insert(m_int_points.end(), remove_mult.begin(), remove_mult.end()); 
            for (auto it = remove_mult.begin(); it != remove_mult.end(); ++it)
            {
                m_int_points.push_back(it->m_int_point);
            }
            return;
        }
    };

    template<typename T, int dim, bool is_rational = false>
    std::vector<Box<T, 4>> bezier_surface_int_bezier_surface(const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>& left, const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>& right, T box_size = TDEFAULT_ERROR<T>::value)
    {
        static_assert(is_rational == false, "is_rational != false");
        auto left_control_points = left.get_control_points();
        int u_left_degree = left_control_points[0].cols() - 1;
        int v_left_degree = left_control_points.rows() - 1;

        auto right_control_points = right.get_control_points();
        int u_right_degree = right_control_points[0].cols() - 1;
        int v_right_degree = right_control_points.rows() - 1;
        std::array<int, 4> degrees{ u_left_degree, v_left_degree, u_right_degree, v_right_degree };

        int count = (v_left_degree + 1) * (u_left_degree + 1) * (u_right_degree + 1) * (v_right_degree + 1);
        Eigen::Matrix<T, dim, Eigen::Dynamic> coeff(dim, count);
        Index<4> index;
        index.reset();
        index.set_bounds(degrees);
        for (int i = 0; i < count; ++i)
        {
            coeff.col(i) = left_control_points[index.m_index[1]].col(index.m_index[0]) - right_control_points[index.m_index[3]].col(index.m_index[2]);
            index.add_one();
        }

        smspe<T, dim, 4> solver;
        solver.init(coeff, degrees);
        std::vector<Box<T, 4>> int_boxes;
        solver.compute_ipp(int_boxes, box_size);
        return int_boxes;
    }




}
