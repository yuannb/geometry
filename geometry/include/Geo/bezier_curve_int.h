#pragma once
#include "bezier_curve.h"
#include "bezier_surface.h"
#include "nurbs_curve.h"
#include "nurbs_surface.h"
#include <array>
#include <deque>
#include <Eigen/Eigenvalues>
#include "convex_hell.h"


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
        Eigen::Matrix<T, equation_count, Eigen::Dynamic> m_coeff;
        mutable Index<variaty_count> m_index;
    public:
        smspe() = default;
        ~smspe() { }

        bool init(const Eigen::Matrix<T, equation_count, Eigen::Dynamic> &coeff, const std::array<int, variaty_count> &degrees)
        {
            int coeff_count = 1;
            for (int degree : degrees)
                coeff_count *= (degree + 1);
            if (coeff_count != coeff.cols())
                return false;
            
            m_index.reset();
            m_index.m_bounds = degrees;
            m_coeff = coeff;
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
                        Box<T, 1> sub_box(Eigen::Vector<T, 1>(current_box.Min[index]), Eigen::Vector<T, 1>(current_box.Max[index]));

                        if (sub_box.Max[0] - sub_box.Min[0] < KNOTS_EPS<T>::value)
                        {
                            Eigen::Vector<T, equation_count> point;
                            nurbs.point_on_curve((sub_box.Max[0] + sub_box.Min[0]) / 2.0, point);
                            for (int i = 0; i < bezier_control_points_count; ++i)
                            {
                                control_points.col(first_index + i * step) = point;
                            }
                        }
                        else
                        {
                            nurbs.sub_divide(sub_box);
                            bezier_control_points = nurbs.get_control_points();
                            // std::cout << bezier_control_points << std::endl;
                            for (int i = 0; i < bezier_control_points_count; ++i)
                            {
                                control_points.col(first_index + i * step) = bezier_control_points.col(i);
                            }
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
      
        bool interval_projected_polyhedron(const Eigen::Matrix<T, equation_count, Eigen::Dynamic> &control_points, Box<T, variaty_count> &reuslt) const
        {
            int points_count = control_points.cols();
            reuslt.Min.setConstant(0.0);
            reuslt.Max.setConstant(1.0);

            m_index.reset();
            for (int f_index = 0; f_index < equation_count; ++f_index)
            {
                Eigen::Matrix<T, variaty_count + 1, Eigen::Dynamic> current(variaty_count + 1, points_count);
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
                for (int index = 0; index < variaty_count; ++index)
                {
                    //project(index, 2) to plane
                    std::vector<Eigen::Vector2<T>> project_point;
                    project_point.reserve(points_count);
                    for (int i = 0; i < points_count; ++i)
                    {
                        Eigen::Vector2<T> p { current(index, i) ,current(variaty_count, i) };
                        project_point.push_back(p);
                    }
                    std::vector<Eigen::Vector2<T>> convex_hell = graham_scan(project_point);
                    std::array<T, 2> interval; 
                    if (convex_hell_int_xaxis(convex_hell, interval) == true)
                    {
                        reuslt.Min[index] = std::max(reuslt.Min[index], interval[0]);
                        reuslt.Max[index] = std::min(reuslt.Max[index], interval[1]);
                    }
                    else
                    {
                        return false;
                    }
                    
                }
            }
            return true;
        }
    
        bool split_box_by_ratio(const Box<T, variaty_count> &box, const Box<T, variaty_count> &origin_box, std::vector<Box<T, variaty_count>> &split_boxes) const
        {
            split_boxes.push_back(box);
            for (int index = 0; index < variaty_count; ++index)
            {
                bool flag = box.Max[index] - box.Min[index] < 0.8;
                if (flag == false && origin_box.Max[index] - origin_box.Min[index] > TDEFAULT_ERROR<T>::value * 10)
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
            std::array<T, 2> nums;
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
                std::swap(intersect_interval[0], intersect_interval[1]);
            return true;
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

    template<typename T, int dim, bool is_rational = false>
    std::vector<Eigen::Vector<T, 3>> bezier_curve_int_bezier_surface(const nurbs_curve<T, dim, is_rational, -1, -1>& left, const nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>& right)
    {
        int left_degree = left.get_degree();
        auto right_control_points = right.get_control_points();
        int u_right_degree = right_control_points[0].cols() - 1;
        int v_right_degree = right_control_points.rows() - 1;
        std::array<int, 3> degrees{ left_degree, u_right_degree, v_right_degree };
        auto left_control_points = left.get_control_points();
        
        int count = (left_degree + 1) * (u_right_degree + 1) * (v_right_degree + 1);
        Eigen::Matrix<T, dim, Eigen::Dynamic> coeff(dim, count);
        Index<3> index;
        index.reset();
        index.set_bounds(degrees);
        for (int i = 0; i < count; ++i)
        {
            coeff.col(i) = left_control_points.col(index.m_index[0]) - right_control_points[index.m_index[2]].col(index.m_index[1]);
            index.add_one();
        }

        smspe<T, dim, 3> solver;
        solver.init(coeff, degrees);
        std::vector<Eigen::Vector<T, 3>> int_params;
        solver.compute_ipp(int_params);
        return int_params;
    }


}
