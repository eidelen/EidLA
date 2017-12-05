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

    size_t n = mat.rows();

    Matrix<double> a(mat);
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

