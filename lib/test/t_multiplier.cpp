#include <gtest/gtest.h>
#include "matrix.hpp"
#include "multiplier.hpp"

TEST(Multiplier, Swap)
{
    auto mat = Matrix<double>::identity(3);

    auto soll = mat;
    soll.swapRows(1, 2);

    auto swapOpMat = Multiplier::swapRow(mat,1,2);
    auto swapRes = swapOpMat*mat;

    ASSERT_TRUE(soll.compare(swapRes));
}

TEST(Multiplier, SwapNonSymetric)
{
    int matData[] = {1,2,  3,4,  5,6};
    auto mat = Matrix<int>(3,2,matData);

    auto soll = mat;
    soll.swapRows(0, 2);

    auto swapOpMat = Multiplier::swapRow(mat,0,2);
    auto swapRes = swapOpMat*mat;

    ASSERT_TRUE(soll.compare(swapRes));
}