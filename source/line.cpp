#include "line.h"
#include <assert.h>
double line::length() const
{
    range krange = get_range();
    assert(krange.isvalid());

    if (krange.get_range_type() == BOUNDED)
        return m_scale * (krange.range_length());
    else
        return infinite;
}

bool line::get_start_point(point3d &point) const
{
    range krange = get_range();
    assert(krange.isvalid());

    double u = krange.get_low();

    point = m_point + m_scale * u * m_vector; 
    //if krange is below unbounded or unbounded, return false
    if (krange.get_range_type() == UNBOUNDED || krange.get_range_type() == BELOW_UNBOUNDED)
        return false;
    
    return krange.isvalid();
}

bool line::get_end_point(point3d &point) const
{
    range krange = get_range();
    assert(krange.isvalid());

    double u = krange.get_height();

    point = m_point + m_scale * u * m_vector; 
    //if krange is below unbounded or unbounded, return false
    if (krange.get_range_type() == UNBOUNDED || krange.get_range_type() == ABOVE_UNBOUNDED)
        return false;
    
    return krange.isvalid();
}

bool line::eval_der(const double u, vector3d **der /*= nullptr*/, 
                            unsigned n_th /*= 0*/, limit_type side /*= unknown_side*/) const
{
    point3d point;
    bool flag = get_point(u, point);
    der[0][0] = point;

    //0-th represent position at param of u of curve
    if (n_th == 0)
        return flag;
    
    der[1][0] = m_scale * m_vector;
    return flag;

}

bool line::get_point(double u, point3d &point) const
{
    bool flag = get_range().isvalid();
    bool flag_1 = get_range().contain(u);
    flag = flag && flag_1;

    //param equation
    point = m_point + m_scale * u * m_vector;
    return flag;
}