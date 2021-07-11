#ifndef KMATH_H
#define KMATH_H

#include "point.h"
#include <vector>
#include <algorithm> 

//evaluate a_n * x^n + .... + a_1 * x + a_0;
point3d Hornerl(const std::vector<point3d> &point_set, const unsigned n, const double u)
{
    point3d result;
    for (int i = n - 1; i >= 0; --i)
    {
        result = result * u + point_set[i];
    }
    return result;
}

//evaluate point of ploynormal surface
point3d Horner2(const std::vector<std::vector<point3d>> &control_point, const unsigned n
                    , const unsigned m, const double u, const double v)
    {
        std::vector<point3d> b(control_point.size());
        for (unsigned i = 0; i <= n; ++i)
            b[i] = Hornerl(control_point[i], n, u);

        point3d result = Hornerl(b, m, v);
        return result;
    }


//evaluate point of bezier curve
point3d deCastejaul(const std::vector<point3d> &control_point, const unsigned degree, const double u)
    {
        std::vector<point3d> result = control_point;
        for (unsigned k = 0; k <= degree; ++k)
        {
            for (unsigned i = 0; i <= degree - k; ++i)
            {
                result[i] = (1 - u) * result[i] + u * result[i + 1];
            }
        }
        return result[0];
    }

//Nrubs book P26, evaluate point if bezier surface
point3d deCastejaul2(const std::vector<std::vector<point3d>> &control_point, const unsigned n
                    , const unsigned m, const double u, const double v)
    {
        if (n <= m)
        {
            std::vector<point3d> Q(m);
            for (unsigned j = 0; j <= m; ++j)
            {
                std::vector<point3d> P(n);
                for (unsigned i = 0; i <= n; ++i)
                {
                    P[i] = control_point[i][j];
                    Q[j] = deCastejaul(P, n, u);
                }
            }
        }
        else
        {
            std::vector<point3d> Q(n);
            for (unsigned i = 0; i <= n; ++i)
                Q[i] = deCastejaul(control_point[i], m, v);
        }
    }
//if knots is valid, return true, else return false
bool isValid(const std::vector<double> &knots, const unsigned p)
{
    //if knots is empty, return false
    if (knots.begin() == knots.end())
        return false;
    
    
    bool flag = false;
    unsigned count = 1;
    double temp = knots[0];
    for (auto it = knots.begin() + 1; it != knots.end(); ++it)
    {
        //konts should be monotone incresed
        if (*it < temp)
            return false;

        if (std::abs(*it - temp) > error)
        {
            if (flag == false)
            {
                if (count != p + 1)
                    return false;
            }
            else
            {
                if (count > p + 1)
                    return false;
            }
            count = 1;
            temp = *it;
        }
        else
            ++count;
    }
    return (count == p + 1);
}


/*
* @brief         evaluate index of range that contains param u  
* @param n       the upper boundary of search
* @param p       the degree of bspline
*/
bool find_span(const unsigned n, const unsigned p, const double u, std::vector<double> &knots, unsigned &index)
{
    //upper boundary should less konts's size
    assert(n <= knots.size());
    
    //konts should be valid
    if (isValid(knots, p) == false)
        return false;
    bool flag = u >= knots[p] && u <= knots[n];
    assert(flag);

    if (u == knots[n + 1])
    {
        index = n;
        return true;
    }
    //because the degree of repeat of first kont of spline is p + 1, so we search begin from pth knots
    unsigned lower = p;
    unsigned upper = n + 1;
    unsigned mid = (lower + upper) / 2;
    while (u < knots[mid] || u >= knots[mid + 1])
    {
        if (u < knots[mid])
            upper = mid;
        else
            lower = mid;
        
        mid = (upper + lower) / 2;
    }
    return mid;
    
}

//the document for this function in basisFuns.svg
bool basisFuns(const unsigned i, const double u, const unsigned p, const std::vector<double> &knots, std::vector<double> &N)
{
    //initial N;
    N.clear();
    N.resize(p + 1);
    N[0] = 1.0;

    //we do not use left[0], right[0]
    std::vector<double> left(p + 1);
    std::vector<double> right(p + 1);

    for (unsigned j = 1; j <= p; ++j)
    {
        left[j] = u - knots[i + 1 - j];
        right[j] = knots[i + j] - u;
        double saved = 0.0;
        for (unsigned r = 0; r <= j; ++r)
        {
            double temp = N[r] / (right[r + 1] + left[j - r]);
            N[r] = saved + right[r + 1] * N[r];
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
    return true;
}
/*
* @brief                                    evaluate nozero basis function and derivatives of bspline
* @param i                                  u contain in [ u_i, u_(i + 1) )
* @param u                                  param u
* @param p                                  degree
* @param unsigned n                         0_th - n_th derivative
* @const std::vector<double> &U             knots vector
* @std::vector<std::vector<double>> &der    0_th - n_th derivative of bspline
*/

//the document for this function in basisFuns.svg
bool DerBasisFuns(const unsigned i, const double u, const unsigned p, const unsigned n
                  , const std::vector<double> &U, std::vector<std::vector<double>> &der)
{
    //ndu[p + 1][p + 1] storage basis function and difference between the node
    std::vector<std::vector<double>> ndu(p + 1);
    for (auto it = ndu.begin(); it != ndu.end(); ++it)
    {
        it->resize(p + 1);
    }

    //a[2][p + 1] storage lasted a^k_j and a^(k - 1)_j
    std::vector<std::vector<double>> a(2);
    for (auto it = a.begin(); it != a.end(); ++it)
    {
        it->resize(p + 1);
    }

    //left[j] = u - U[i + 1 - j]
    //right[j] = u_(i + j) - u
    std::vector<double> left(p + 1);
    std::vector<double> right(p + 1);//will waste storage left[0] and right[0]



    for (unsigned j = 1; j <= p; ++j)
    {
        left[j] = u - U[i + 1 - j];
        right[j] = U[i + j] - u;
        double saved = 0.0;
        for (unsigned r = 0; r < j; ++r)
        {
            ndu[j][r] = right[r + 1] - left[j - r];
            double temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    //0_th der
    for (unsigned j = 0; j <= p; ++j)
    {
        der[0][j] = ndu[j][p];
    }

    //evaluate k_th der of N_p^(i - p + r)
    for (unsigned r = 0; r <= p; ++r)
    {
        unsigned s1 = 0, s2 = 1;
        a[0][0] = 1.0;

        //evaluate n_th der
        for (unsigned k = 1; k < n; ++k)
        {
            double d = 0.0;
            double rk = r - k;
            double pk = p - k;
            unsigned j1 = 0;
            unsigned j2 = 0;
            if (r >= k)
            {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            j1 = rk >= -1 ? 1 : -rk;
            j2 = r - 1 <= pk ? k - 1; p -r;
            for (unsigned j = j1; j <= j2; ++j)
            {
                a[s2][j] = (a[s1][j] - a[s1][j - 1] / ndu[pk + 1][rk + j]);
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk)
            {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            der[k][j] = d;
            std::swap(s1, s2);
        }
    }

    //coefficient
    unsigned r = p;
    for (unsigned k = 1; k <= n; ++k)
    {
        for (unsigned j = 0; j <= p; ++j)
            der[k][j] *= r;
        r *= (p - k);
    }
    return true;
}

#endif  //KMATH_H

