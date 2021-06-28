#include "curve.h"

range::range(double low, double height, range_type type)
{
    m_low = low;
    m_height = height;
    m_type = type;
}

range::range(double low, double height)
{
    m_low = low;
    m_height = height;

    //if height < low, let m_type = UNKNOWE
    if ((height + error) < low)
        m_type = UNKNOWN_RANGE;
    else
        m_type = BOUNDED;
}

range::range(const range &krange)
{
    m_low = krange.m_low;
    m_height = krange.m_height;
    m_type = krange.m_type;
}

range& range::operator=(const range &krange)
{
    m_low = krange.m_low;
    m_height = krange.m_height;
    m_type = krange.m_type;

    return *this;
}

bool range::isvalid() const
{
    if (m_type == UNKNOWN_RANGE)
        return false;
    
    if (m_type != BOUNDED)
        return true;
    
    if (m_height + error < m_low)
        return false;
    
    return true;
}


bool range::set_height(double height)
{
    m_height = height;
    return isvalid();
}

bool range::set_low(double low)
{
    m_low = low;
    return isvalid();
}

bool range::contain(double u) const
{
    if (m_type == UNKNOWN_RANGE)
        return false;
    else if (m_type == BOUNDED)
        return (m_low - error < u) && (m_height + error > u);
    else if (m_type == ABOVE_UNBOUNDED)
        return (m_low - error < u);
    else
        return (m_height + error > u);
}