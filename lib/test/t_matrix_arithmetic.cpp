#include <gtest/gtest.h>
#include "matrix.hpp"

TEST(Matrix, Addition)
{
    size_t nRows = 5;
    size_t nCols = 10;

    auto m1 = Matrix<int>(nRows, nCols);
    auto m2 = Matrix<int>(nRows, nCols);

    m1.fill(2);
    m2.fill(4);

    auto res = m1 + m2;
    for (size_t m = 0; m < nRows; m++)
        for (size_t n = 0; n < nCols; n++)
            ASSERT_EQ(res(m, n), 6);

    auto res2 = m2 + m1;
    for (size_t m = 0; m < nRows; m++)
        for (size_t n = 0; n < nCols; n++)
            ASSERT_EQ(res2(m, n), 6);
}

TEST(Matrix, Subtraction)
{
    size_t nRows = 5;
    size_t nCols = 10;

    auto m1 = Matrix<int>(nRows, nCols);
    auto m2 = Matrix<int>(nRows, nCols);

    m1.fill(2);
    m2.fill(4);

    auto res = m1 - m2;
    for (size_t m = 0; m < nRows; m++)
        for (size_t n = 0; n < nCols; n++)
            ASSERT_EQ(res(m, n), -2);

    auto res2 = m2 - m1;
    for (size_t m = 0; m < nRows; m++)
        for (size_t n = 0; n < nCols; n++)
            ASSERT_EQ(res2(m, n), 2);
}

TEST(Matrix, Comparisson)
{
    auto m0 = Matrix<int>::identity(4);
    auto m1 = Matrix<int>::identity(4);

    ASSERT_TRUE(m0.compare(m1));
    ASSERT_TRUE(m1.compare(m0));

    m1(0, 0) = 2;

    ASSERT_FALSE(m0.compare(m1));
    ASSERT_FALSE(m1.compare(m0));

    auto f0 = Matrix<double>::identity(4);
    auto f1 = Matrix<double>::identity(4);

    ASSERT_TRUE(f0.compare(f1));
    ASSERT_TRUE(f1.compare(f0));

    f0(0, 0) = 1.1;

    ASSERT_FALSE(f0.compare(f1));
    ASSERT_FALSE(f1.compare(f0));

    ASSERT_TRUE(f0.compare(f1,true,0.2));
    ASSERT_TRUE(f1.compare(f0,true,0.2));
}

TEST(Matrix, MaxElement)
{
    int         inData[6] = {1, 3, 4, 3, 2, 2}; // max is 4 at pos 0,2
    Matrix<int> in(2, 3, inData);

    auto   max = in.max();
    size_t m   = std::get<0>(max);
    size_t n   = std::get<1>(max);
    int    val = std::get<2>(max);

    ASSERT_EQ(m, 0);
    ASSERT_EQ(n, 2);
    ASSERT_EQ(val, 4);
}

TEST(Matrix, MinElement)
{
    int         inData[6] = {1, 3, 4, -5, 2, 2}; // max is -5 at pos 1,0
    Matrix<int> in(2, 3, inData);

    auto   min = in.min();
    size_t m   = std::get<0>(min);
    size_t n   = std::get<1>(min);
    int    val = std::get<2>(min);

    ASSERT_EQ(m, 1);
    ASSERT_EQ(n, 0);
    ASSERT_EQ(val, -5);
}

TEST(Matrix, RandomInt)
{
    auto rand = Matrix<int>::random(4, 6, 0, 10);

    ASSERT_EQ(4, rand.rows());
    ASSERT_EQ(6, rand.cols());

    for (size_t s = 0; s < rand.getNbrOfElements(); s++)
    {
        int val = rand.data()[s];
        ASSERT_GE(val, 0);
        ASSERT_LE(val, 10);
    }
}

TEST(Matrix, RandomDouble)
{
    auto rand = Matrix<double>::random(4, 6, 0.0, 1.0);

    ASSERT_EQ(4, rand.rows());
    ASSERT_EQ(6, rand.cols());

    for (size_t s = 0; s < rand.getNbrOfElements(); s++)
    {
        double val = rand.data()[s];
        ASSERT_GE(val, 0.0);
        ASSERT_LE(val, 1.0);
    }
}

TEST(Matrix, Norm)
{
    // column vector
    int data[] = {1,2,3};
    auto vec = Matrix<int>(3,1,data);

    ASSERT_NEAR(vec.norm(), 3.7417, 0.001);

    // row vector
    ASSERT_NEAR(vec.transpose().norm(), 3.7417, 0.001);

    auto mat = Matrix<double>::random(2,2,1.0,2.0);
    ASSERT_FLOAT_EQ(mat.norm(), 0.0);
}

TEST(Matrix, NormSquare)
{
    // column vector
    int data[] = {1,2,3};
    auto vec = Matrix<int>(3,1,data);

    ASSERT_EQ(vec.normSquare(), 1*1 + 2*2 + 3*3 );

    // row vector
    ASSERT_EQ(vec.transpose().normSquare(), 1*1 + 2*2 + 3*3);

    // 0 for matrix
    auto mat = Matrix<double>::random(2,2,1.0,2.0);
    ASSERT_FLOAT_EQ(mat.normSquare(), 0.0);
}

TEST(Matrix, SetToZero)
{
    int data[] = {1,2, 3,-1};
    auto mat = Matrix<int>(2,2,data);

    int sooldata[] = {0,2, 3,0};
    auto soll = Matrix<int>(2,2,sooldata);

    mat.setToZero(2);
    ASSERT_TRUE( mat.compare(soll));
}

TEST(Matrix, SetToDouble)
{
    double data[] = {0.01, 0.2,   -0.5, -0.04};
    auto mat = Matrix<double>(2,2,data);

    double sooldata[] = {0.0, 0.2,   -0.5, 0.0};
    auto soll = Matrix<double>(2,2,sooldata);

    mat.setToZero(0.06);
    ASSERT_TRUE( mat.compare(soll));
}

TEST(Matrix, Absolute)
{
    double data[] = {-1.0, 2.0,   5.0, -3.0};
    auto mat = Matrix<double>(2,2,data);

    double sollData[] = {1.0, 2.0,   5.0, 3.0};
    auto soll = Matrix<double>(2,2,sollData);

    soll.compare( mat.absolute() );
}

TEST(Matrix, Mean)
{
    double data[] = {-1, -3, 2, 4};
    auto mat = Matrix<double>(2,2,data);

    ASSERT_FLOAT_EQ(0.5, mat.mean());
}

TEST(Matrix, Normalize)
{
    double data[] = {0,1,2,3};
    auto mat = Matrix<double>(2,2,data);

    double sollData[] = {-1, -1.0/3.0, 1.0/3.0, 1};
    auto soll = Matrix<double>(2,2,sollData);

    double scale, mean;

    auto normMat = mat.normalize(mean,scale);

    ASSERT_TRUE( soll.compare( normMat));
    ASSERT_FLOAT_EQ( 1.5, mean );
    ASSERT_FLOAT_EQ( 1.0/3.0 * 2.0 , scale);

    // recover original matrix
    auto denormMat = normMat.denormalize(mean,scale);
    ASSERT_TRUE( mat.compare(denormMat));
}

TEST(Matrix, NormalizeSpecial)
{
    // same values
    double data[] = {5, 5, 5, 5};
    auto mat = Matrix<double>(2,2,data);

    double sollData[] = {0, 0, 0, 0};
    auto soll = Matrix<double>(2,2,sollData);

    double scale, mean;

    auto normMat = mat.normalize(mean,scale);

    ASSERT_TRUE( soll.compare( normMat));
    ASSERT_FLOAT_EQ( 5, mean );
    ASSERT_FLOAT_EQ( 1.0 , scale);

    // recover original matrix
    auto denormMat = normMat.denormalize(mean,scale);
    ASSERT_TRUE( mat.compare(denormMat));
}

TEST(Matrix, NormalizeBatch)
{
    for( int i = 0; i < 100; i++ )
    {
        // same values
        auto mat = Matrix<double>::random(5,5,-100,100);

        double scale, mean;
        auto normMat = mat.normalize(mean, scale);

        // recover original matrix
        auto denormMat = normMat.denormalize(mean, scale);
        ASSERT_TRUE(mat.compare(denormMat, true, 0.0001));
    }
}
