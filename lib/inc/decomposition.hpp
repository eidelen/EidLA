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
        LUResult(Matrix<double> l, Matrix<double> u) :
                L(l), U(u) {}
        const Matrix<double> L;
        const Matrix<double> U;
    };

public:

    template <class T>
    static LUResult luDecomposition(const Matrix<T>& mat);

private:

    static LUResult doolittle(const Matrix<double>& a);

};


// Infos from:
//  http://mathonline.wikidot.com/the-algorithm-for-doolittle-s-method-for-lu-decompositions

template <class T>
Decomposition::LUResult Decomposition::luDecomposition(const Matrix<T>& mat)
{
    if( mat.rows() != mat.cols() )
    {
        std::cout << "Square matrix required";
        std::exit(-1);
    }

    Matrix<double> a(mat);
    return doolittle(a);
}

Decomposition::LUResult Decomposition::doolittle(const Matrix<double>& a)
{
    size_t n = a.rows();
    Matrix<double> l = Matrix<double>(n,n);
    Matrix<double> u = Matrix<double>(n,n);

    l.setToIdentity(); u.fill(0.0);

    for( size_t k = 0; k < n; k++ )
    {
        for( size_t m = k; m < n; m++ )
        {
            // compute u(k,m)
            double pSum = 0;

            if( k > 0 )
            {
                for (size_t j = 0; j < k; j++ )
                {
                    pSum += l(k,j)*u(j,m);
                }
            }
            u(k,m) = a(k,m) - pSum;
        }

        for( size_t i = k+1; i < n; i++ )
        {
            // compute l(i,k)
            if( k == i )
            {
                l(i,k) = 1.0;
            }
            else
            {
                double pSum = 0;

                if( k > 0 )
                {
                    for (size_t j = 0; j < k; j++ )
                    {
                        pSum += l(i,j)*u(j,k);
                    }
                }
                l(i,k) = (a(i,k) - pSum) / u(k,k);
            }
        }
    }

    LUResult ret(l,u);
    return ret;
}

#endif //MY_DECOMPOSITION_H

