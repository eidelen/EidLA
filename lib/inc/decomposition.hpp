/****************************************************************************
** Copyright (c) 2017 Adrian Schneider
**
** Permission is hereby granted, free of charge, to any person obtaining a
** copy of this software and associated documentation files (the "Software"),
** to deal in the Software without restriction, including without limitation
** the rights to use, copy, modify, merge, publish, distribute, sublicense,
** and/or sell copies of the Software, and to permit persons to whom the
** Software is furnished to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in
** all copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
** FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
** DEALINGS IN THE SOFTWARE.
**
*****************************************************************************/

#ifndef MY_DECOMPOSITION_H
#define MY_DECOMPOSITION_H

#include "matrix.hpp"

class Decomposition
{
public:
    struct LUResult
    {
        LUResult(Matrix<double> l, Matrix<double> u)
        : L(l), U(u)
        {
        }
        const Matrix<double> L; // Lower triangle matrix
        const Matrix<double> U; // Upper triangle matrix
    };

    struct EigenPair
    {
        EigenPair(Matrix<double> v, double l)
        : V(v), L(l)
        {
        }
        const Matrix<double> V; // Eigen vector
        const double         L; // Eigen value
    };

public:
    /**
     * LU decomposition of the matrix mat.
     * @param mat
     * @return LUResult
     */
    template <class T>
    static LUResult luDecomposition(const Matrix<T>& mat);

    /**
     * Eigen decomposition of the matrix mat.
     * @param mat
     * @return Eigen pairs in desending eigenvalue order.
     */
    template <class T>
    static std::vector<EigenPair> eigen(const Matrix<T>& mat);

    /**
     * Compute Rayleigh quotient of a matrix and a vector. This can be
     * used to find the Eigenvalue to a corresponding Eigenvector and
     * matrix.
     * @param mat
     * @param vec
     * @return Rayleigh quotient
     */
    template <class T>
    static double rayleighQuotient(const Matrix<T>& mat, const Matrix<T> vec);

private:
    template <class T>
    static LUResult doolittle(const Matrix<T>& a);
};

// Infos from:
//  http://mathonline.wikidot.com/the-algorithm-for-doolittle-s-method-for-lu-decompositions

template <class T>
Decomposition::LUResult Decomposition::luDecomposition(const Matrix<T>& mat)
{
    if (mat.rows() != mat.cols())
    {
        std::cout << "Square matrix required";
        std::exit(-1);
    }

    Matrix<double> a(mat);
    return doolittle(a);
}

template <class T>
Decomposition::LUResult Decomposition::doolittle(const Matrix<T>& a)
{
    size_t         n = a.rows();
    Matrix<double> l = Matrix<double>(n, n);
    Matrix<double> u = Matrix<double>(n, n);

    l.setToIdentity();
    u.fill(0.0);

    for (size_t k = 0; k < n; k++)
    {
        for (size_t m = k; m < n; m++)
        {
            // compute u(k,m)
            double pSum = 0;

            if (k > 0)
            {
                for (size_t j = 0; j < k; j++)
                {
                    double lkj = l(k, j);
                    double ujm = u(j, m);
                    pSum += l(k, j) * u(j, m);
                }
            }

            // debug
            double akm = a(k, m);
            double ukm = a(k, m) - pSum;

            u(k, m) = a(k, m) - pSum;
        }

        for (size_t i = k + 1; i < n; i++)
        {
            // compute l(i,k)
            if (k == i)
            {
                l(i, k) = 1.0;
            }
            else
            {
                double pSum = 0;

                if (k > 0)
                {
                    for (size_t j = 0; j < k; j++)
                    {
                        pSum += l(i, j) * u(j, k);
                    }
                }

                double ukk = u(k, k);
                if (std::abs(ukk) > std::numeric_limits<double>::min()) // check if divider > 0
                {
                    l(i, k) = (a(i, k) - pSum) / ukk;
                }
                else
                {
                    // this is a singular matrix, therefore l(i,k) can be freely chosen -> lu decomposition is not unique
                    // info: https://math.stackexchange.com/questions/2010470/doolittle-transformation-is-non-unique-for-singular-matrices
                    l(i, k) = 0;
                }
            }
        }
    }

    LUResult ret(l, u);
    return ret;
}

template <class T>
std::vector<Decomposition::EigenPair> Decomposition::eigen(const Matrix<T>& mat)
{
    std::vector<EigenPair> ret;

    // Power iteration : https://en.wikipedia.org/wiki/Power_iteration
    // Find biggest eigenvalue
    Matrix<double> matD      = mat;
    double         initRange = 1.0 / std::sqrt(mat.cols());
    Matrix<double> ev        = Matrix<double>::random(matD.cols(), 1, -initRange, initRange);
    Matrix<double> ev_before = ev;

    bool go = true;
    while (go)
    {
        Matrix<double> prod = matD * ev;
        ev                  = prod.normalizeColumns();

        // check stopping criteria -> stop if all entries almost equal
        go        = !(ev.compare(ev_before));
        ev_before = ev;
    }

    double eigenVal = rayleighQuotient(matD, ev);
    ret.push_back(EigenPair(ev, eigenVal));

    return ret;
}

// info: https://www.mathematik.uni-wuerzburg.de/~borzi/RQGradient_Chapter_10.pdf
template <class T>
double Decomposition::rayleighQuotient(const Matrix<T>& m, const Matrix<T> v)
{
    Matrix<T> vT = v.transpose();
    return static_cast<double>((vT * m * v)(0, 0)) / static_cast<double>((vT * v)(0, 0));
}

#endif //MY_DECOMPOSITION_H
