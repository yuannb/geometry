#pragma once

namespace tnurbs
{
    //目前仅在计算切向锥的时候使用;
    template<typename T, int dim>
    struct cone
    {
        //半角
        T m_angle;
        Eigen::Vector<T, dim> m_dir;
        cone(const cone &c) : m_angle(c.m_angle), m_dir(c.m_dir) {  }
        cone(const T &angle) { m_angle = angle; m_dir.setConstant(0.0); m_dir[0] = 1.0;}
        cone(const T &angle, const Eigen::Vector<T, dim> &dir) : m_angle(angle), m_dir(dir) { }
        cone() { m_angle = 0.0; m_dir.setConstant(0.0); m_dir[0] = 1.0; }
        ~cone() { };
        ENUM_NURBS merge_vector(const Eigen::Vector<T, dim> &v)
        {
            if (v.norm() < TDEFAULT_ERROR<T>::value)
            {
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            Eigen::Vector<T, dim> unit_v = v;
            unit_v.normalize();

            T cos_angle = m_dir.dot(unit_v);
            T angle;
            if (std::abs(cos_angle) > 1.0 - TDEFAULT_ERROR<T>::value)
            {
                if (cos_angle > 1.0 /*- TDEFAULT_ERROR<T>::value*/ && cos_angle < 1.0 + TDEFAULT_ERROR<T>::value)
                    angle = 0.0;
                else if (cos_angle < -1.0 /*+ TDEFAULT_ERROR<T>::value*/ && cos_angle > -1.0 - TDEFAULT_ERROR<T>::value)
                    angle = M_PI;
                else
                    return ENUM_NURBS::NURBS_ERROR;
            }
            else
                angle = std::acos(cos_angle);
            if (angle < m_angle)
                return ENUM_NURBS::NURBS_SUCCESS;
            //else
            Eigen::Vector<T, dim> u_normal = unit_v - m_dir.dot(unit_v) * m_dir;
            if (u_normal.norm() < TDEFAULT_ERROR<T>::value)
            {
                int max_index = 0;
                T max_element = 0.0;
                for (int index = 0; index < dim; ++index)
                {
                    if (std::abs(m_dir[index]) > max_element)
                    {
                        max_element = std::abs(m_dir[index]);
                        max_index = index;
                    }                        
                }
                u_normal.setConstant(0.0);
                int next_index = max_index % dim;
                u_normal[max_index] = -1.0 * m_dir[next_index];
                u_normal[next_index] = m_dir[max_index];
            }
            u_normal.normalize();
            m_dir = (std::cos(m_angle) + std::cos(angle)) * m_dir + (std::sin(angle) - std::sin(m_angle)) * u_normal;
            m_dir.normalize();
            
            m_angle += angle;
            if (m_angle > M_PI)
                m_dir *= (-1.0);
            m_angle /= 2.0;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        ENUM_NURBS merge_cone(const cone &c)
        {
            Eigen::Vector<T, dim> unit_v = c.m_dir;

            T cos_angle = m_dir.dot(unit_v);
            T angle;
            if (std::abs(cos_angle) > 1.0 - TDEFAULT_ERROR<T>::value)
            {
                if (cos_angle > 1.0 - TDEFAULT_ERROR<T>::value && cos_angle < 1.0 + TDEFAULT_ERROR<T>::value)
                    angle = 0.0;
                else if (cos_angle < -1.0 + TDEFAULT_ERROR<T>::value && cos_angle > -1.0 - TDEFAULT_ERROR<T>::value)
                    angle = M_PI;
                else
                    return ENUM_NURBS::NURBS_ERROR;
            }
            else
                angle = std::acos(cos_angle);
            if (angle + c.m_angle < m_angle)
                return ENUM_NURBS::NURBS_SUCCESS;
            //else
            Eigen::Vector<T, dim> u_normal = unit_v - m_dir.dot(unit_v) * m_dir;
            if (u_normal.norm() < TDEFAULT_ERROR<T>::value)
            {
                int max_index = 0;
                T max_element = 0.0;
                for (int index = 0; index < dim; ++index)
                {
                    if (std::abs(m_dir[index]) > max_element)
                    {
                        max_element = std::abs(m_dir[index]);
                        max_index = index;
                    }                        
                }
                u_normal.setConstant(0.0);
                int next_index = max_index % dim;
                u_normal[max_index] = -1.0 * m_dir[next_index];
                u_normal[next_index] = m_dir[max_index];
            }
            u_normal.normalize();
            // Eigen::Vector<T, dim> v1 = std::cos(m_angle) * m_dir - std::sin(m_angle) * u_normal;
            // v1.normalize();
            // Eigen::Vector<T, dim> v2 = std::cos(angle)
            m_dir = (std::cos(m_angle) + std::cos(angle + c.m_angle)) * m_dir + (std::sin(angle + c.m_angle) - std::sin(m_angle)) * u_normal;
            m_dir.normalize();
            
            m_angle += (angle + c.m_angle);
            if (m_angle > M_PI)
                m_dir *= (-1.0);
            m_angle /= 2.0;
            return ENUM_NURBS::NURBS_SUCCESS;
        }
        
        bool is_cotain_vector(const Eigen::Vector<T, dim>& dir) const
        {
            Eigen::Vector<T, dim> normal_dir = dir.normalized();
            double dot_vecs = normal_dir.dot(m_dir);
            double cos_theta = std::cos(m_angle);
            if (dot_vecs >= cos_theta - PRECISION<T>::value)
            {
                return true;
            }
            return false;
        }


        cone & operator=(const cone &c)
        {
            if (this != &c)
            {
                m_angle = c.m_angle;
                m_dir = c.m_dir;
            }
            return *this;
        }
    
    };


    //以point为奇异点, 包含box的锥;
    template<typename T, int dim>
    cone<T, dim> point_box(const Box<T, dim> &box, const Eigen::Vector<T, dim> &point)
    {
        struct Index
        {
            Eigen::Vector<int, dim> index;
            Index() { index.setConstant(0); }
            ENUM_NURBS add_one()
            {
                index[dim - 1] += 1;
                for (int i = dim - 1; i > 0; --i)
                {
                    if (index[i] != 2)
                        break;
                    else
                    {
                        index[i] -= 2;
                        index[i - 1] += 1;
                    }
                }
                if (index[0] == 2)
                    return ENUM_NURBS::NURBS_ERROR;
                return ENUM_NURBS::NURBS_SUCCESS;
            }
            ENUM_NURBS get_point(const Box<T, dim> &box, Eigen::Vector<T, dim> &point) const
            {
                for (int i = 0; i < dim; ++i)
                {
                    if (index[i] == 0)
                        point[i] = box.Min[i];
                    else
                        point[i] = box.Max[i];
                }             
                return ENUM_NURBS::NURBS_SUCCESS;
            }
        };
        Index idx;
        Eigen::Vector<T, dim> vec;
        idx.get_point(box, vec);
        vec -= point;
        vec.normalize();
        cone<T, dim> result(0.0, vec);
        while (idx.add_one())
        {
            idx.get_point(box, vec);
            vec -= point;
            result.merge_vector(vec);
        }
        return result;
    }

}

