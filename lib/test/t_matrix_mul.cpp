#include <gtest/gtest.h>
#include "matrix.hpp"

TEST(MatMul, Scalar)
{
    size_t s = 5;
    auto mat = Matrix<int>(s,s);
    mat.setValue(3);

    auto res1 = mat * 2; // res should be 2*3 = 6

    auto soll1 = Matrix<int>(s,s);
    soll1.setValue(6);

    ASSERT_TRUE(res1.compare(soll1));

    auto res2 = 2 * mat; // res should be 2*3 = 6
    ASSERT_TRUE(res2.compare(soll1));
}

/*
TEST(MatMul, MatrixMatrixSquare)
{
    size_t s = 4;

    auto m1 = Matrix<int>(s,s); m1.setValue(2);
    auto m2 = Matrix<int>(s,s); m2.setValue(3);

    auto res1 = m1 * m2; // 2*3*s = 24

    auto soll1 = Matrix<int>(s,s); soll1.setValue(24);
    ASSERT_TRUE(res1.compare(soll1));


    m2(3,3) = 1;
    auto res2 = m1 * m2; // 2*3*s = 24 and res2(3,3) = 2*3*3 + 1*2 = 20;
    std::cout << m1;
    std::cout << m2;
    std::cout << res2;

    auto soll2 = Matrix<int>(s,s); soll2.setValue(24); soll2(0,3) = 20; soll2(1,3) = 20; soll2(2,3) = 20; soll2(3,3) = 20;
    ASSERT_TRUE(res2.compare(soll2));

    std::cout << res2;
}*/




