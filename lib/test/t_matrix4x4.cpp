#include <gtest/gtest.h>
#include "matrix4x4.hpp"

TEST(Matrix4x4, Constructor)
{
    Matrix4x4 mat = Matrix4x4();

    ASSERT_EQ(mat.cols(), 4);
    ASSERT_EQ(mat.rows(), 4);
    ASSERT_TRUE(Matrix<double>::identity(4).compare(mat));
}

TEST(Matrix4x4, RotZ)
{
    Matrix4x4 mat = Matrix4x4();
    mat.rotZ(M_PI/2.0); // rotate 90

    double v1Data[] = {2,1,7};
}