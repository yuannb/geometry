#include "bezier.h"

//copy construct
bezier::bezier(const bezier &lhs)
{
    m_degree = lhs.m_degree;
    m_control_number = lhs.m_control_number;
    m_control_point = lhs.m_control_point;
}


//copy-assignment
bezier& bezier::operator=(const bezier &lhs)
{
    m_degree = lhs.m_degree;
    m_control_number = lhs.m_control_number;
    m_control_point = lhs.m_control_point;
    return *this;
}

//calculate B^i_n
double bezier::bernstein(unsigned i, double u) const   //NURBS book P12
    {
        assert (u + error >= 0 && u - error <= 1);

        std::vector<double> result(m_degree + 1, 0.0);
        result[m_degree - i] = 1.0;
        double u1 = 1 - u;
        for (unsigned k = 1; k <= m_degree; ++k)
        {
            for (unsigned j = m_degree; j >= k; --j)
            {
                result[j] = u1 * result[j] + u * result[j - 1];
            }
        }
        return result[m_degree];
    }

//evaluate all value of polynomial of bernstein of degree n
std::vector<double> bezier::all_bernstein(double u) const  //NURBS P13
{
    assert (u + error >= 0 && u - error <= 1);

    std::vector<double> result(0,  + 1);
    result[0] = 1.0;
    double u1 = 1 -u;

    for (unsigned j = 1; j <= m_degree; ++j)
    {
        double save = 0.0;
        for (int k = 0; k < j; ++k)
        {
            double temp = result[k];
            result[k] = u1 * temp + u * save;
            save = temp;
        }
        result[j] = u * save;
    }
    return result;
}

bool bezier::get_point(double v, point3d &point) const
    {
        std::vector<double> B = all_bernstein(v);
        point3d point(0.0, 0.0, 0.0);
        for (unsigned k = 0; k <= m_degree; ++k)
        {
            point = point + (B[k] * m_control_point[k]);
        }
        return true;
    }