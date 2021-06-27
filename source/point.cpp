#include <iostream>
#include "point.h"

point3d::point3d(const vector3d &vector)
{
    m_x = vector.getX();
    m_y = vector.getY();
    m_z = vector.getZ();
}

point3d::point3d(const point3d &point)
{
    m_x = point.m_x;
    m_y = point.m_y;
    m_z = point.m_z;
}

point3d& point3d::operator=(const point3d &point)
{
    m_x = point.m_x;
    m_y = point.m_y;
    m_z = point.m_z;
    return *this;
}

point3d& point3d::operator=(const vector3d &vector)
{
    *this = point3d(vector);
    return *this;
}

point3d point3d::operator+(const point3d &point) const
{
    point3d result;

    result.setX(m_x + point.m_x);
    result.setY(m_y + point.m_y);
    result.setZ(m_z + point.m_z);

    return result;
}

point3d point3d::operator+(const vector3d &vector) const
{
    return(*this + point3d(vector));
}

point3d point3d::operator-(const point3d &point) const
{
    point3d result;

    result.setX(m_x - point.m_x);
    result.setY(m_y - point.m_y);
    result.setZ(m_z - point.m_z);

    return result;
}

point3d point3d::operator-(const vector3d &vector) const
{
    return (*this - point3d(vector));
}
bool point3d::set_point(const point3d &point)
{
    m_x = point.m_x;
    m_y = point.m_y;
    m_z = point.m_z;
    return true;
}

bool point3d::set_point(const vector3d &vector)
{
    this->set_point(point3d(vector));
    return true;
}


vector3d::vector3d(const vector3d &vector)
{
    m_x = vector.m_x;
    m_y = vector.m_y;
    m_z = vector.m_z;
}

vector3d::vector3d(const point3d &point)
{
    m_x = point.getX();
    m_y = point.getY();
    m_z = point.getZ();
}

vector3d& vector3d::operator=(const vector3d &vector)
{
    m_x = vector.m_x;
    m_y = vector.m_y;
    m_z = vector.m_z;
}

vector3d& vector3d::operator=(const point3d &point)
{
    *this = vector3d(point);
    return *this;
}

vector3d vector3d::operator+(const point3d &point) const
{
    return (*this + vector3d(point));
}

vector3d vector3d::operator+(const vector3d &vector) const
{
    vector3d result;
    result.m_x = m_x + vector.m_x;
    result.m_y = m_y + vector.m_y;
    result.m_z = m_z + vector.m_z;
    return result;
}

vector3d vector3d::operator-(const vector3d &vector) const
{
    vector3d result;
    result.m_x = m_x - vector.m_x;
    result.m_y = m_y - vector.m_y;
    result.m_z = m_z - vector.m_z;
    return result;
}

vector3d vector3d::operator-(const point3d &point) const
{
    return (*this - vector3d(point));
}

bool vector3d::set_vector(const point3d &point)
{
    *this = vector3d(point);
    return true;
}

bool vector3d::set_vector(const vector3d &vector)
{
    *this = vector;
    return true;
}

double vector3d::dot(const vector3d &vector) const
{
    double result;
    result += m_x * vector.m_x;
    result += m_y * vector.m_y;
    result += m_z * vector.m_z;
    return result;
}

vector3d vector3d::cross(const vector3d &vector) const
{
    vector3d result;
    result.m_x = m_y * vector.m_z - m_z * vector.m_z;
    result.m_y = m_z * vector.m_x - m_x * vector.m_z;
    result.m_z = m_x * vector.m_y - m_y * vector.m_x;
    return result;
}

