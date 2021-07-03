#if !defined(BEZIER_H)
#define BEZIER_H

#include <vector>
#include "curve.h"
#include <assert.h>
#include "kMath.h"

class bezier : public curve
{
private:
    //degree of bezier curve
    unsigned m_degree;
    //the number of contorl point
    unsigned m_control_number;
    //control point
    std::vector<point3d> m_control_point;

public:
    //default construct
    bezier() = default;
    
    //construct : use degree and control point
    bezier(const unsigned degree, const std::vector<point3d> &control_point) :
            m_degree(degree), m_control_number(static_cast<unsigned> (control_point.size()))
            , m_control_point(control_point) { };

    //construct : use degree and contorl point , control point will default initialize
    bezier(const unsigned degree, const unsigned point_number) :
            m_degree(degree), m_control_number(point_number) {  m_control_point.resize(point_number); }

    //copy construct
    bezier(const bezier &lhs);

    //deconstruct
    ~bezier() { };

    //copy- 
    bezier &operator=(const bezier &lhs);

    //calculate B^i_n
    double bernstein(unsigned i, double u) const;   //NURBS book P12


    //evaluate all value of polynomial of bernstein of degree n
    std::vector<double> all_bernstein(double u) const;  //NURBS P13

    virtual bool get_point(double v, point3d &point) const override;

};



#endif // BEZIER_H
