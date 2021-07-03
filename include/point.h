#if !defined(POINT_H)
#define POINT_H

#include <iostream>
#include "entity.h"

class vector3d;

//will become template
class point3d
{
friend point3d operator*(const double scale, const point3d &point);
private:
    double m_x;
    double m_y;
    double m_z;
public:
    //default construct
    point3d() : m_x(0.0), m_y(0.0), m_z(0.0) { };

    //construct
    point3d(const double x, const double y, const double z) :
        m_x(x), m_y(y), m_z(z) {}
    
    //construct : vector-->point
    point3d(const vector3d &vector);

    //destruct
    ~point3d();

    //copy construct
    point3d(const point3d &point);

    //copy-assignment
    point3d &operator=(const point3d &point);
    //copy-assignment : vector --> point
    point3d &operator=(const vector3d &vector);

    //overload
    point3d operator+(const point3d &point) const;
    point3d operator+(const vector3d &vector) const;

    point3d operator*(const double sacle) const;

    point3d operator-(const point3d &point) const;
    point3d operator-(const vector3d &vector) const;

    //set point
    bool set_point(const point3d &point);
    //set point : vector --> point
    bool set_point(const vector3d &vector);

    //set coordinate
    bool setX(const double x) { m_x = x; return true; }
    bool setY(const double y) { m_y = y; return true; }
    bool setZ(const double z) {m_z = z; return true; }

    //get coordinate
    double getX() const { return m_x; }
    double getY() const { return m_y; }
    double getZ() const { return m_z; }
};

class vector3d
{

friend vector3d operator*(const double u, const vector3d vec);

private:
    double m_x;
    double m_y;
    double m_z;
public:
    //default construct
    vector3d() : m_x(0.0), m_y(0.0), m_z(0.0) { };

    //construct
    vector3d(const double x, const double y, const double z) :
        m_x(x), m_y(y), m_z(z) { }

    //construct : point --> vector
    vector3d(const point3d &point);

    //copy construct
    vector3d(const vector3d &vector);
    
    //decontruct
    ~vector3d() { };

    //copy-assignment
    vector3d &operator=(const vector3d &vector);

    //copy-assignment
    vector3d &operator=(const point3d &point);

    //overload
    vector3d operator+(const point3d &point)const ;
    vector3d operator+(const vector3d &vector) const ;

    vector3d operator*(const double scale) const ;
    vector3d operator/(const double scale) const;

    vector3d operator-(const point3d &point) const;
    vector3d operator-(const vector3d &vector) const ;

    //set point
    bool set_vector(const point3d &point);
    //set point : vector --> point
    bool set_vector(const vector3d &vector);

    //set coordinate
    bool setX(const double x) { m_x = x; return true; }
    bool setY(const double y) { m_y = y; return true; }
    bool setZ(const double z) {m_z = z; return true; }

    //get coordinate
    double getX() const { return m_x; }
    double getY() const { return m_y; }
    double getZ() const { return m_z; }

    vector3d cross(const vector3d &vector) const;
    double dot(const vector3d &vector) const;
};

#endif // POINT