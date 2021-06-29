#if !defined(CURVE_H)
#define CURVE_H
#include <iostream>
#include "curvetype.h"
#include "point.h"


enum range_type
    {
        UNKNOWN_RANGE = 0,
        BOUNDED = 1,
        ABOVE_UNBOUNDED = 2,
        BELOW_UNBOUNDED = 3,
        UNBOUNDED = 4
    };

class range
{
private:
     double m_low;
     double m_height;
     range_type m_type;
public:
    //default construct
    range() : m_low(0.0), m_height(1.0), m_type(UNBOUNDED) { };

    //construct
    range(double low, double height, range_type type);

    //construct
    range(double low, double height);

    //copy construct
    range(const range &krange);
    //copy-assignment
    range& operator=(const range &krange);

    //deconstruct
    ~range() {};

    //valid
    bool isvalid() const;

    //get low of range
    double get_low() const { return m_low; };

    //if range is not valid, return false
    double get_height() const { return m_height; }

    //set low : if [low, m_hight] is not valid, return false
    bool set_low(double low);

    //set height : if [m_low, height] is not valid, return false
    bool set_height(double height);

    //set range type : if {[m_low, m_height], type} is invaid, return false
    bool set_type(range_type type);

    //get range type
    int get_range_type() { return m_type; }

    double range_length() { return m_height - m_low; }

    //if m_range contain u, return true; else return false;
    bool contain(double u) const;
};


enum limit_type
{
    both_side = 0,      //the limit on both side
    left_side = 1,      //the limit on left side
    right_side = 2,     //the limit on right side
    unknown_side = 3    //unknown side
};


class curve
{
private:
    unsigned m_type = 0;  //curve type
    range m_range;
public:
    //calculate length of curve
    virtual double length() const = 0;
    
    //get start point
    virtual bool get_start_point(point3d &point) const = 0;

    //get end point
    virtual bool get_end_point(point3d &point) const = 0;
    
    //get point at parameter v
    virtual bool get_point(double v, point3d &point) const = 0;

    //derivative
    virtual bool eval_der(const double u, vector3d **der = nullptr, 
                            unsigned n_th = 0, limit_type side = unknown_side) const = 0;

    //get range
    range get_range() const {return m_range; }

    //set range
    virtual bool set_range(const range &krange) = 0;

    //set_range_type
    virtual bool set_range_type(const range_type krange_type) = 0;

    //eval param u such that \alpha(u) = point
    virtual bool eval_param(const point3d &point, double &param) const = 0;
};






#endif // CURVE


