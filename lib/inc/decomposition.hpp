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


template <class T>
Decomposition::LUResult Decomposition::luDecomposition(const Matrix<T>& mat)
{
    if( mat.rows() != mat.cols() )
    {
        std::cout << "Square matrix required";
        std::exit(-1);
    }

    size_t s = mat.rows();

    Matrix<double> l = Matrix<double>(s,s);
    Matrix<double> u = Matrix<double>(s,s);

    l.fill(0.0); u.fill(0.0);

    u(0,0) = mat(0,0);
    l.setToIdentity();

    for( size_t m = 0; m < s; m++ )
    {
        for(size_t n = 0; n < s; m++ )
        {
            // compute l_mn
            if( m > n )
            {
                double accum = 0;
                for (size_t p = 0; p < n - 1; p++) {
                    accum += l(m, p) * u(p, n);
                }


                l(m, n) = (mat(m, n) - accum) / u(n, n);
            }



            double accum = 0;
            for( size_t p = 0; p < m-1; p++ )
            {
                accum += l(m,p)*u(p,n);
            }

            u(m,n) = mat(m,n) - accum;

        }
    }



    LUResult ret(l,u);

    return ret;
}

#endif //MY_DECOMPOSITION_H

