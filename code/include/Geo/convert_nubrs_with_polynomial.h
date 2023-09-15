#pragma once

#include "nurbs_curve.h"
#include "polynomial_curve.h"

template<typename T, int dim, bool is_rational>
ENUM_NURBS convert_nubrs_to_polynomial(const nurbs_curve<T, dim, is_rational, -1, -1> &old_curve,
    polynomial_curve<T, dim, is_rational, -1, -1>)
{
    
}