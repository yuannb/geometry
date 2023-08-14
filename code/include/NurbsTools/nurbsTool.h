#pragma once
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <array>
#include <vector>

#include <iostream>

//曲线的DeCasteljaul算法
template<int dim, int degree>
bool DeCasteljaul(double u, const Eigen::Vector2d &m_interval, Eigen::Matrix<double, dim, degree + 1> m_control_points, Eigen::Vector<double, dim> &point)
{
    Eigen::Matrix<double, dim, degree + 1> control_points(m_control_points);
    for (int cureentDegree = 1; cureentDegree < degree + 1; ++cureentDegree)
    {
        for (int col = 0; col <= degree - cureentDegree; ++col)
            control_points.col(col) = (m_interval[1] - u) * control_points.col(col) + u * control_points.col(col + 1);
    }
    point = control_points.col(0);
    return true;
}

template<int degree>
bool AllBernstein(double u, const Eigen::Vector2d &m_interval, Eigen::Vector<double, degree + 1> &bs)
{
    bs[0] = 1.0;
    double u1 = m_interval[1] - u;
    for (int j = 1; j <= degree; ++j)
    {
        double saved = 0.0;
        for (int k = 0; k < j; ++k)
        {
            double temp = bs[k];
            bs[k] = saved + u1 * temp;
            saved = u * temp;
        }
        bs[j] = saved;
    }
    return true;
}


//曲面的DeCasteljaul算法
template<int dim, int du, int dv>
bool DeCasteljaul2(double u, double v, const Eigen::Vector2d &u_interval, const Eigen::Vector2d &v_interval, 
                    std::array<Eigen::Matrix<double, dim, du + 1>, dv + 1> control_points, Eigen::Vector<double, dim> &point)
{
    Eigen::Matrix<double, dim, dv + 1> points;
    for (int j = 0; j <= dv; ++j)
    {
        // // Eigen::Matrix<double, dim, du + 1> mat = control_points[j].transpose();
        Eigen::Vector<double, dim> tpoint;
        DeCasteljaul<dim, du>(u, u_interval, control_points[j], tpoint);
        points.col(j) = tpoint;
    }
    DeCasteljaul<dim, dv>(v, v_interval, points, point);
    return true;
}

//返回区间左端点的序号(用户保证u必须在knots之中)
bool FindSpan(int degree, double u, const Eigen::VectorX<double> &knots, int &index)
{
    //degree == 0时，不处理这种情况
    if (degree == 0) return false;
    int knotsSize = knots.size();
    if (knotsSize < 2)  return false;
    if (knotsSize < 2 * (degree + 1)) return false;
    if (u > knots[knotsSize - 1] || u < knots[0]) return false;
    if (u == knots[knotsSize - 1])
    {
        index = knotsSize - 2;
        return true;
    }
    int low = degree;
    int high = knotsSize - degree;
    int mid = (low + high) / 2;

    while (u < knots[mid] || u >= knots[mid + 1])
    {
        if (u < knots[mid]) high = mid;
        else                low = mid;
        mid = (low + high) / 2;
    }
    index = mid;
    return true;
}

bool BasisFunc(int i, double u, int p, const Eigen::VectorX<double> &knots, Eigen::VectorX<double> &basisValues)
{
    basisValues[0] = 1.0;
    std::vector<double> left(p + 1, 0);
    std::vector<double> right(p + 1, 0);
    for (int j = 1; j <= p; ++j)
    {
        left[j] = u - knots[i + 1 - j];
        right[j] = knots[i + j] - u;
        double saved = 0.0;
        
        //可以使用稀疏矩阵优化
        for (int r = 0; r < j; ++r)
        {
            double temp = basisValues[r] / (right[r + 1] + left[j - r]);
            basisValues[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        basisValues[j] = saved;
    }
    return true;
}

//用户保证dersMat的空间
bool DersBasisFuns(int i, double u, int n, int p, Eigen::VectorX<double> &knots, Eigen::MatrixX<double> &dersMat)
{
    //将函数值和节点值存到数组中
    Eigen::MatrixX<double> ndu(p + 1, p + 1);
    ndu(0, 0) = 1.0;
    Eigen::VectorX<double> left(p + 1);
    Eigen::VectorX<double> right(p + 1);
    Eigen::Vector<double, 1> uv(u);
    left.block(1, 0, p, 1) = (knots(Eigen::seq(i, i + 1 - p, Eigen::fix<-1>)) * -1).rowwise() + uv;
    right.block(1, 0, p, 1) = knots(Eigen::seq(i + 1, i + p)).rowwise() - uv;
    for (int r = 0; r < p; ++r)
    {
        ndu.block(r + 1, r, p - r, 1) =  left(Eigen::seq(1, p - r)).rowwise() + Eigen::Vector<double, 1>(right[r + 1]);
    }
    for (int j = 1; j <= p; ++j)
    {
        double saved = 0.0;
        for (int r = 0; r < j; ++r)
        {
            double temp = ndu(r, j - 1) / ndu(j, r);
            ndu(r, j) = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu(j, j) = saved;
    }

    //计算导数
    dersMat.block(0, 0, p + 1, 1) = ndu.block(0, p, p + 1, 1);

    for (int r = 0; r <= p; ++r) //函数下标
    {
        Eigen::MatrixX<double> Akj(p + 1, 2);
        Akj.setConstant(0.0);
        Akj(0, 0) = 1.0;
        int s1 = 0, s2 = 1;

        //循环计算k阶导数
        for (int k = 1; k <= n; ++k)
        {
            double d = 0.0;
            int rk = r - k;
            int pk = p - k;
            if (r >= k)
            {
                Akj(0, s2) = Akj(0, s1) / ndu(pk + 1, rk);
                d = Akj(0, s2) * ndu(rk, pk);
            }
            int j1 = rk >= -1 ? 1 : -rk;
            int j2 = r - 1 <= pk ? k - 1 : p - r; 
            int jj = j2 - j1 + 1;
            jj = std::max(0, jj);
            Akj.block(j1, s2, jj, 1) = (Akj.block(j1, s1, jj, 1) - Akj.block(j1 - 1, s1, jj, 1)).array() / ndu.block(pk + 1, rk + j1, jj, 1).array();
            d += (ndu.block(rk + j1, pk, 1, jj) * Akj.block(j1, s2, jj, 1))(0, 0);
            if (r <= pk)
            {
                Akj(k, s2) = -Akj(k - 1, s1) / ndu(pk + 1, r);
                d += Akj(k, s2) * ndu(r, pk);
            }
            dersMat(r, k) = d;
            std::swap(s1, s2);
        }

    }
    int r = p;
    for (int k = 1; k <= n; ++k)
    {
        for (int j = 0; j <= p; ++j)
            dersMat(j, k) *= r;
        r *= (p - k);
    }
    return true;
}

//用户保证ndu的空间
bool AllBasisFuns(int i, double u, int p, Eigen::VectorX<double> &knots, Eigen::MatrixX<double> &ndu)
{
    //将函数值和节点值存到数组中
    ndu(0, 0) = 1.0;
    Eigen::VectorX<double> left(p + 1);
    Eigen::VectorX<double> right(p + 1);
    Eigen::Vector<double, 1> uv(u);
    left.block(1, 0, p, 1) = (knots(Eigen::seq(i, i + 1 - p, Eigen::fix<-1>)) * -1).rowwise() + uv;
    right.block(1, 0, p, 1) = knots(Eigen::seq(i + 1, i + p)).rowwise() - uv;
    for (int r = 0; r < p; ++r)
    {
        ndu.block(r + 1, r, p - r, 1) =  left(Eigen::seq(1, p - r)).rowwise() + Eigen::Vector<double, 1>(right[r + 1]);
    }
    for (int j = 1; j <= p; ++j)
    {
        double saved = 0.0;
        for (int r = 0; r < j; ++r)
        {
            double temp = ndu(r, j - 1) / ndu(j, r);
            ndu(r, j) = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu(j, j) = saved;
    }
    return true;
}


//用户确保p和i给的是合理的
bool OneBasisFun(int p, const Eigen::VectorX<double> &knots, int i, double u, double &Nip)
{
    int knotsSize = knots.size();

    //保证是nurbs, 前后节点重复p + 1次
    if (knotsSize < 2 * p + 2)
        return false;

    if ( (i == 0 && u == knots[0]) ||
        (i == knotsSize - p - 2 && u == knots[knotsSize - 1]))
    {
        Nip = 1.0;
        return true;
    }

    if (u < knots[i] || u >= knots[i + p + 1])
    {
        Nip = 0.0;
        return true;
    }

    Eigen::VectorX<double> N(p + 1);
    for (int j = 0; j <= p; ++j)    //初始化0次基函数
    {
        if (u >= knots[i + j] && u < knots[i + j + 1])
            N[j] = 1.0;
        else
            N[j] = 0.0;
    }

    Eigen::Vector<double, 1> uv(u);
    Eigen::VectorX<double> left = (- 1 * knots.block(i, 0, p + 1, 1)).rowwise() + uv;
    Eigen::VectorX<double> right = knots.block(i + 1, 0, p + 1, 1).rowwise() - uv;
    Eigen::VectorX<double> ndu(p + 1);
    for (int k = 1; k <= p; ++k)
    {  
        int nduCount = p + 2 - k;
        for (int index = 0; index < nduCount; ++index)
        {
            ndu[index] = right[index + k - 1] + left[index];
            if (ndu[index] == 0.0)
                ndu[index] = 10;
        }
        N.block(0, 0, p + 1 - k, 1) = left.block(0, 0, p + 1 - k, 1).array() * N.block(0, 0, p + 1 - k, 1).array() / ndu.block(0, 0, p + 1 - k, 1).array() + 
                                      right.block(k, 0, p + 1 - k, 1).array() * N.block(1, 0, p + 1 - k, 1).array() / ndu.block(1, 0, p + 1 - k, 1).array();
    }
    Nip = N[0];
    return true;
}

bool DersOneBasisFun(int p, const Eigen::VectorX<double> &knots, int i, double u, int n, Eigen::VectorX<double> &ders)
{
    int knotsSize = knots.size();
    //保证是nurbs, 前后节点重复p + 1次
    if (knotsSize < 2 * p + 2)
        return false;

    if (u < knots[i] || u >= knots[i + p + 1])
    {
        ders.setConstant(0.0);
        return true;
    }

    Eigen::MatrixX<double> N(p + 1, p + 1);
    //初始化0次基函数
    for (int j = 0; j <= p; ++j)
    {
        if (u >= knots[i + j] && u < knots[i + j + 1])
            N(j, 0) = 1.0;
        else
            N(j, 0) = 0.0;
    }
    
    //迭代基函数的第一列到第p列
    Eigen::Vector<double, 1> uv(u);
    Eigen::VectorX<double> left = (- 1 * knots.block(i, 0, p + 1, 1)).rowwise() + uv;
    Eigen::VectorX<double> right = knots.block(i + 1, 0, p + 1, 1).rowwise() - uv;
    Eigen::VectorX<double> ndu(p + 1);
    for(int k = 1; k <= p; ++k)
    {
        int nduCount = p + 2 - k;
        for (int index = 0; index < nduCount; ++index)
        {
            ndu[index] = right[index + k - 1] + left[index];
            if (ndu[index] == 0.0)
                ndu[index] = 10;
        }
        N.block(0, k, p + 1 - k, 1) = left.block(0, 0, p + 1 - k, 1).array() * N.block(0, k - 1, p + 1 - k, 1).array() / ndu.block(0, 0, p + 1 - k, 1).array() + 
                                      right.block(k, 0, p + 1 - k, 1).array() * N.block(1, k - 1, p + 1 - k, 1).array() / ndu.block(1, 0, p + 1 - k, 1).array();
    }
    ders[0] = N(0, p);

    std::cout << N << std::endl;
    //计算各阶导数
    for (int k = 1; k <= n; ++k)
    {
        Eigen::VectorX<double> ND = N.block(0, p - k, k + 1, 1);
        //计算宽度为k的表, 循环列
        for (int jj = 1; jj <= k; ++jj)
        {
            Eigen::VectorX<double> denominator(k + 2 - jj);
            int count = k + 2 - jj;
            for (int index = 0; index < count; ++index)
            {
                denominator[index] = left[index] + right(p - k + jj + index - 1);
                if (denominator[index] == 0.0)
                    denominator[index] = 1.0;
            }
            ND.block(0, 0, k + 2 - jj, 1) = (ND.block(0, 0, k + 2 - jj, 1).array() / denominator.array()).rowwise() * Eigen::Array<double, 1, 1>(p - k + jj) ;
            ND.block(0, 0, k + 1 - jj, 1) = ND.block(0, 0, k + 1 - jj, 1) - ND.block(1, 0, k + 1 - jj, 1);
        }
        ders[k] = ND[0];
    }
    return true;
}

