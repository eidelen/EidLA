#ifndef MY_TRANSFORMATION_H
#define MY_TRANSFORMATION_H

#include "matrix.hpp"

class Transformation
{
public:

    template <class T>
    static Matrix<double> echelon(const Matrix<T>& mat);

    template <class T>
    static Matrix<double> echelonReduced(const Matrix<T>& mat);
};



// Algorithm described in
// http://stattrek.com/matrix-algebra/echelon-transform.aspx#MatrixA

template <class T>
Matrix<double> Transformation::echelon(const Matrix<T>& mat)
{
    Matrix<double> ret(mat);

    for( size_t n = 0; n < ret.cols(); n++ )
    {
        // find first non zero entry in col n
        size_t pivotRow = 0;
        bool foundNonZeroPivot = false;
        for( size_t m = n; m < ret.rows(); m++ )
        {
            if( std::fabs( ret(m,n) ) > 1e-5 )
            {
                pivotRow = m;
                foundNonZeroPivot = true;
                break;
            }
        }

        if( foundNonZeroPivot )
        {
            if( pivotRow > n )
                ret.swapRows(n,pivotRow); // move pivot up

            // adapt pivot line
            double pivotElement = ret(n,n);
            auto pivotRow = ret.row(n);
            auto scaledPivotRow = pivotRow*(1/pivotElement);
            ret.setRow(n, scaledPivotRow);

            // add scaled pivot line to below rows so that elements in column n become zero
            for( size_t q = n+1; q < ret.rows(); q++ )
            {
                double localPivotElement = ret(q,n);
                double localPivotFactor = localPivotElement / pivotElement * (-1);
                auto newLocalRow = ( pivotRow*localPivotFactor) + ret.row(q);
                ret.setRow(q, newLocalRow);
            }
        }
    }

    return ret;
}

template <class T>
Matrix<double> Transformation::echelonReduced(const Matrix<T>& mat)
{

}



#endif //MY_TRANSFORMATION_H

