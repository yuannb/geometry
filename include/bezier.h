#if !defined(BEZIER_H)
#define BEZIER_H

#include <vector>
#include "curve.h"
#include <assert.h>

class bezier : public curve
{
private:
    //degree of bezier curve
    unsigned degree;
    //the number of contorl point
    unsigned control_point_number;
    //control point
    std::vector<point3d> control_point;

public:
    //calculate B^i_n
    double betnstein(unsigned i, unsigned degree, double u)   //NURBS book P12
    {
        assert(degree >= i);
        assert (u + error >= 0 && u - error <= 1);

        std::vector<double> result(degree + 1, 0.0);
        result[degree - i] = 1.0;
        double u1 = 1 - u;
        for (unsigned k = 1; k <= degree; ++k)
        {
            for (unsigned j = degree; j >= k; --j)
            {
                result[j] = u1 * result[j] + u * result[j - 1];
            }
        }
        return result[degree];
    }
};



#endif // BEZIER_H
