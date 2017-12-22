#include <gtest/gtest.h>
#include "matrix.hpp"

TEST(Multiplier, Swap)
{
    auto mat = Matrix<double>::identity(3);

    auto soll = mat;
    soll.swapRows(1, 2);

    auto swapOpMat = Multiplier::swapRow(mat, 1, 2);
    auto swapRes   = swapOpMat * mat;

    ASSERT_TRUE(soll.compare(swapRes));
}

TEST(Multiplier, SwapNonSymetric)
{
    int  matData[] = {1, 2, 3, 4, 5, 6};
    auto mat       = Matrix<int>(3, 2, matData);

    auto soll = mat;
    soll.swapRows(0, 2);

    auto swapOpMat = Multiplier::swapRow(mat, 0, 2);
    auto swapRes   = swapOpMat * mat;

    ASSERT_TRUE(soll.compare(swapRes));
}

TEST(Multiplier, MulAdd)
{
    auto mat = Matrix<int>(4, 4);
    mat.setToIdentity();

    // add 9*row 2 to row 3.
    auto soll = mat;
    soll.setRow(3, soll.row(3) + (mat.row(2) * 9));

    auto res = Multiplier::addProductOfRow(mat, 9, 2, 3) * mat;

    ASSERT_TRUE(soll.compare(res));
}

TEST(Multiplier, MulAddNonSymetric)
{
    int  matData[] = {1, 2, 3, 4, 5, 6};
    auto mat       = Matrix<int>(3, 2, matData);

    // add -8*row 0 to row 2.
    auto soll = mat;
    soll.setRow(2, soll.row(2) + (mat.row(0) * (-8)));

    auto res = Multiplier::addProductOfRow(mat, -8, 0, 2) * mat;

    ASSERT_TRUE(soll.compare(res));
}

TEST(Multiplier, MultiplyRow)
{
    int  matData[] = {1, 2, 3, 4, 5, 6};
    auto mat       = Matrix<int>(3, 2, matData);

    // row 1 x 6
    auto soll = mat;
    soll.setRow(1, soll.row(1) * 6);

    auto res = Multiplier::multiplyRow(mat, 6, 1) * mat;

    ASSERT_TRUE(soll.compare(res));
}
