#ifndef MY_TRANSFORMATION_H
#define MY_TRANSFORMATION_H

#include "matrix.hpp"

class Transformation
{
public:

    template <class T>
    static Matrix<double> echelon(const Matrix<T>& mat);

};


// Algorithm described in
// http://stattrek.com/matrix-algebra/echelon-transform.aspx#MatrixA

template <class T>
Matrix<double> Transformation::echelon(const Matrix<T>& mat)
{
    Matrix<double> ret(mat);

    size_t processingRow = 0;
    for( size_t n = 0; n < ret.cols(); n++ )
    {
        // find first non zero entry in col n from row processingRow on
        size_t pivotRow = 0;
        bool foundNonZeroPivot = false;
        for( size_t m = processingRow; m < ret.rows(); m++ )
        {
            if( std::fabs( ret(m,n) ) > std::numeric_limits<double>::min() )
            {
                pivotRow = m;
                foundNonZeroPivot = true;
                break;
            }
        }

        if( foundNonZeroPivot )
        {
            // if the found pivot row is not equal the current processing row,
            // swap these two row.
            if( pivotRow > processingRow )
                ret.swapRows(processingRow,pivotRow); // move pivot up

            // adapt pivot line
            double pivotElement = ret(processingRow,n);
            auto pivotRow = ret.row(processingRow);
            auto scaledPivotRow = pivotRow*(1/pivotElement);
            ret.setRow(processingRow, scaledPivotRow);

            // add scaled pivot line to below rows so that elements in column n become zero
            for( size_t q = processingRow+1; q < ret.rows(); q++ )
            {
                double localPivotElement = ret(q,n);
                double localPivotFactor = localPivotElement / pivotElement * (-1);
                auto newLocalRow = ( pivotRow*localPivotFactor) + ret.row(q);
                ret.setRow(q, newLocalRow);
            }

            processingRow++;
        }
    }

    return ret;
}


#endif //MY_TRANSFORMATION_H

