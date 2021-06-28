#if !defined(LINE)
#define LINE
#include "curve.h"

//  \alpha(t) = m_point + scale * param * m_vector

class line : public curve
{
private:

    //one position of line
    point3d m_point;

    //direction vector, unit vector
    vector3d m_vector;

    //scaling
    double m_scale;

public:
    //default construct
    line() : m_point(0, 0, 0), m_vector(1, 0, 0), m_scale(1.0) { };

    //construct
    line(const point3d &point, const vector3d &vector, double scale) { }

    //copy construct
    line(const line &kline);

    //copy-assignment
    line& operator=(const line &kline);

    //deconstruct
    ~line();

    //memeber function

    //lenght
    virtual double length() const override;

    //start point
    virtual bool get_start_point(point3d &point) const override;

    //get end point
    virtual bool get_end_point(point3d &point) const override;

    //get point at parameter v
    virtual bool get_point(double v, point3d &point) const override;

    //derivative
    virtual bool eval_der(const double u, vector3d **der = nullptr, 
                            unsigned n_th = 0, limit_type side = unknown_side) const override;

    //line type
    virtual bool curve_type() { return LINE; }

};



#endif // LINE
