#include <gtest/gtest.h>
#include "matrix.hpp"

TEST(Matrix, InitializationCheckSizes)
{
    size_t nRows = 5;
    size_t nCols = 10;
    size_t nElements = nRows * nCols;

    MatrixISP m( new Matrix<int>(nRows,nCols) );

    ASSERT_EQ(m->rows(), nRows);
    ASSERT_EQ(m->cols(), nCols);
    ASSERT_EQ(m->getNbrOfElements(), nElements);
}

TEST(Matrix, SetGetSingleValue)
{
    size_t nRows = 5; size_t nCols = 10; int currentVal = 2;
    Matrix<int> mat(nRows, nCols);

    for( size_t m = 0; m < nRows; m++ )
    {
        for (size_t n = 0; n < nCols; n++)
        {
            currentVal++;
            mat.setValue(m,n,currentVal);
            ASSERT_EQ(mat.getValue(m, n), currentVal);

            currentVal++;  //operators
            mat(m,n) = currentVal;
            ASSERT_EQ(mat(m, n), currentVal);
        }
    }
}


TEST(Matrix, SetElementsToSameValue)
{
    size_t nRows = 5; size_t nCols = 10; int val = 99;

    MatrixISP mat( new Matrix<int>(nRows,nCols) );
    mat->setValue(val);

    for( size_t m = 0; m < nRows; m++ )
        for( size_t n = 0; n < nCols; n++ )
            ASSERT_EQ(mat->getValue(m, n), val);
}


TEST(Matrix, Identity)
{
    size_t s = 5;
    Matrix<int> eyeM = Matrix<int>::eye(s);

    for(size_t m = 0; m < s; m++)
        for(size_t n = 0; n < s; n++)
            if( n == m )
                ASSERT_EQ(eyeM(m,n), 1);
            else
                ASSERT_EQ(eyeM(m,n), 0);
}

TEST(Matrix, Row)
{
    auto mat = Matrix<int>::eye(3);

    auto r3 = mat.row(2); // r3 -> 0,0,1

    auto s3 = Matrix<int>(1,3);
    s3(0,0) = 0; s3(0,1) = 0; s3(0,2) = 1;

    ASSERT_TRUE(s3.compare(r3));
}

TEST(Matrix, Column)
{
    auto mat = Matrix<int>::eye(3);

    auto r3 = mat.column(2); // r3 -> 0;0;1

    auto s3 = Matrix<int>(3,1);
    s3(0,0) = 0; s3(1,0) = 0; s3(2,0) = 1;

    ASSERT_TRUE(s3.compare(r3));
}
