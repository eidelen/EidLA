#include <gtest/gtest.h>
#include "matrix4x4.hpp"

TEST(Matrix4x4, Constructor)
{
    Matrix4x4 mat = Matrix4x4();

    ASSERT_EQ(mat.cols(), 4);
    ASSERT_EQ(mat.rows(), 4);
}