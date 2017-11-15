#include <gtest/gtest.h>
#include "matrix.hpp"

TEST(Matrix, Addition)
{
    size_t nRows = 5;
    size_t nCols = 10;

    auto m1 = Matrix<int>(nRows,nCols);
    auto m2 = Matrix<int>(nRows,nCols);

    m1.fill(2);
    m2.fill(4);


    auto res = m1 + m2;
    for(size_t m = 0; m < nRows; m++)
        for(size_t n = 0; n < nCols; n++)
            ASSERT_EQ(res(m,n),6);

    auto res2 = m2 + m1;
    for(size_t m = 0; m < nRows; m++)
        for(size_t n = 0; n < nCols; n++)
            ASSERT_EQ(res2(m,n),6);
}

TEST(Matrix, Subtraction)
{
    size_t nRows = 5;
    size_t nCols = 10;

    auto m1 = Matrix<int>(nRows,nCols);
    auto m2 = Matrix<int>(nRows,nCols);

    m1.fill(2);
    m2.fill(4);

    auto res = m1 - m2;
    for(size_t m = 0; m < nRows; m++)
        for(size_t n = 0; n < nCols; n++)
            ASSERT_EQ(res(m,n),-2);

    auto res2 = m2 - m1;
    for(size_t m = 0; m < nRows; m++)
        for(size_t n = 0; n < nCols; n++)
            ASSERT_EQ(res2(m,n),2);
}

TEST(Matrix, Comparisson)
{
    auto m0 = Matrix<int>::identity(4);
    auto m1 = Matrix<int>::identity(4);

    ASSERT_TRUE(m0.compare(m1));
    ASSERT_TRUE(m1.compare(m0));

    m1(0,0) = 2;

    ASSERT_FALSE(m0.compare(m1));
    ASSERT_FALSE(m1.compare(m0));


    auto f0 = Matrix<double>::identity(4);
    auto f1 = Matrix<double>::identity(4);

    ASSERT_TRUE(f0.compare(f1,0.05));
    ASSERT_TRUE(f1.compare(f0,0.05));

    f0(0,0) = 1.01;


    ASSERT_FALSE(f0.compare(f1));
    ASSERT_FALSE(f1.compare(f0));

}