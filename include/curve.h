#if !defined(CURVE)
#define CURVE
#include <iostream>
#include "curvetype.h"
#include "point.h"


class range
{
public:
    enum range_type
    {
        UNKNOWN_RANGE = 0,
        BOUNDED = 1,
        ABOVE_UNBOUNDED = 2,
        BELOW_UNBOUNDED = 3,
        UNBOUNDED = 4
    };
private:
     double m_low;
     double m_height;
     range_type m_type;
public:
    //default construct
    range() : m_low(0.0), m_height(0.0), m_type(BOUNDED) { };

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

    //if range is not valid, return false
    bool get_low(double &low) const;

    //if range is not valid, return false
    bool get_height(double &height) const;

    //set low : if [low, m_hight] is not valid, return false
    bool set_low(double low);

    //set height : if [m_low, height] is not valid, return false
    bool set_height(double height);

    //set range type : if {[m_low, m_height], type} is invaid, return false
    bool set_type(range_type type);
};

class curve
{
private:
    unsigned m_type = 0;  //curve type
    range m_range;
public:
    //calculate length of curve
    virtual double length() = 0;
    
    //get start point
    virtual bool get_start_point(point3d &point) = 0;

    //get start point
    virtual bool get_end_point(point3d &point) = 0;
    
    //get point at parameter v
    virtual bool get_point(double v, point3d &point) = 0;

};






#endif // CURVE


