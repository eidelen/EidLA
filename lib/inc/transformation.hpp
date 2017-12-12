#ifndef MY_TRANSFORMATION_H
#define MY_TRANSFORMATION_H

#include "matrix.hpp"

class Transformation
{
public:

/**
 * Computes the Echelon form of matrix mat.
 * @param mat
 * @return Echelon form.
 */
    template <class T>
    static Matrix<double> echelon(const Matrix<T>& mat);

/**
 * Computes the Echelon form of matrix mat.
 * @param mat
 * @param rowOperations List of required row operations at return.
 * @return Echelon form.
 */
    template <class T>
    static Matrix<double> echelon(const Matrix<T>& mat, std::vector<Matrix<double>>& rowOperations);

/**
 * Computes the reduced Echelon form of matrix mat.
 * @param mat
 * @return Echelon form.
 */
    template <class T>
    static Matrix<double> reduced_echelon(const Matrix<T>& mat);

/**
 * Computes the reduced Echelon form of matrix mat.
 * @param mat
 * @param rowOperations List of required row operations at return.
 * @return Echelon form.
 */
template <class T>
static Matrix<double> reduced_echelon(const Matrix<T>& mat, std::vector<Matrix<double>>& rowOperations);

};



// Algorithm described in
// http://stattrek.com/matrix-algebra/echelon-transform.aspx#MatrixA

template <class T>
Matrix<double> Transformation::echelon(const Matrix<T>& mat)
{
    std::vector<Matrix<double>> ops; // not used
    return Transformation::echelon(mat, ops);
}

template <class T>
Matrix<double> Transformation::echelon(const Matrix<T>& mat, std::vector<Matrix<double>>& rowOperations)
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
            {
                ret.swapRows(processingRow, pivotRow); // move pivot up
            }

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


template <class T>
Matrix<double> Transformation::reduced_echelon(const Matrix<T>& mat)
{
    std::vector<Matrix<double>> rowOps; // not used
    return reduced_echelon(mat, rowOps);
}


template <class T>
Matrix<double> Transformation::reduced_echelon(const Matrix<T>& mat, std::vector<Matrix<double>>& rowOperations)
{
    Matrix<double> echMat = Transformation::echelon(mat, rowOperations);

    if( echMat.rows() < 2 )
    {
        return echMat;
    }

    // going from down up
    size_t processingRow = echMat.rows()-1; // rows is min = 2!
    while( processingRow > 0 ) // when processing row is 1, the last row edited is 0 (first row)
    {
        // find first non-zero entry searching from left to right in row processingRow
        // nonZeroCol will be the kinda pivot element
        size_t nonZeroCol = 0; bool nonZeroEntryFound = false;
        for( size_t i = 0; i < echMat.cols(); i++ )
        {
            if( std::abs(echMat(processingRow,i)) > 0.0 )
            {
                nonZeroCol = i;
                nonZeroEntryFound = true;
                break;
            }
        }

        if( nonZeroEntryFound )
        {
            // modify rows above
            double pivotElement = echMat(processingRow, nonZeroCol);

            for( size_t m = 0; m < processingRow; m++ )
            {
                double rowFactor = - echMat(m, nonZeroCol) / pivotElement;
                echMat.setRow(m,  echMat.row(m) + (echMat.row(processingRow) * rowFactor)  );
            }
        }

        processingRow--;
    }

    return echMat;
}



#endif //MY_TRANSFORMATION_H

