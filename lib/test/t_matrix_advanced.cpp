#include <gtest/gtest.h>
#include "matrix.hpp"

// http://stattrek.com/matrix-algebra/how-to-find-inverse.aspx
TEST(MatAdvanced, Inverse)
{
    double inData[] = {1, 2, 2, 2, 2, 2, 2, 2, 1};
    auto   mat      = Matrix<double>(3, 3, inData);

    double sollData[] = {-1, 1, 0, 1, -1.5, 1, 0, 1, -1};
    auto   soll       = Matrix<double>(3, 3, sollData);

    bool successfull     = false;
    auto computedInverse = mat.inverted(&successfull);

    ASSERT_TRUE(soll.compare(computedInverse));
    ASSERT_TRUE(successfull);

    // inv(mat) * mat = identity
    auto identRes  = computedInverse * mat;
    auto sollIdent = Matrix<double>::identity(3);
    ASSERT_TRUE(sollIdent.compare(identRes));
}

TEST(MatAdvanced, NonInvertable)
{
    double inData[] = {1, 2, 2, 2, 2, 2};
    auto   mat      = Matrix<double>(3, 2, inData);

    bool successfull     = true;
    auto computedInverse = mat.inverted(&successfull);

    ASSERT_FALSE(successfull);
}

TEST(MatAdvanced, InvertSingularMatrix)
{
    auto mat = Matrix<int>(3, 3);
    mat.fill(2);

    bool successfull     = true;
    auto computedInverse = mat.inverted(&successfull);

    ASSERT_FALSE(successfull);
}

// checked with octave
TEST(MatAdvanced, Determinant1)
{
    double inData[] = {1, 2, 2, 2, 0, -1, -2, 1, 3};
    auto   mat      = Matrix<double>(3, 3, inData);

    bool   ok;
    double det = mat.determinant(&ok);

    ASSERT_NEAR(det, -3.0, 0.00001);
    ASSERT_TRUE(ok);
}

// example from https://www.mathsisfun.com/algebra/matrix-determinant.html
TEST(MatAdvanced, Determinant2)
{
    double inData[] = {6, 1, 1, 4, -2, 5, 2, 8, 7};
    auto   mat      = Matrix<double>(3, 3, inData);

    bool   ok;
    double det = mat.determinant(&ok);

    ASSERT_NEAR(det, -306.0, 0.00001);
    ASSERT_TRUE(ok);
}

TEST(MatAdvanced, DeterminantIdentity)
{
    auto mat = Matrix<double>::identity(10);
    bool ok;
    ASSERT_NEAR(mat.determinant(&ok), 1.0, 0.00001);
}

TEST(MatAdvanced, DeterminantSingularMatrix)
{
    auto mat = Matrix<int>(3, 3);
    mat.fill(2);
    bool ok;
    ASSERT_NEAR(mat.determinant(&ok), 0.0, 0.00001);
}

TEST(MatAdvanced, ColumnNormalization)
{
    auto mat = Matrix<int>(1, 1);
    mat(0, 0) = 4;

    ASSERT_DOUBLE_EQ(1.0, mat.normalizeColumns()(0, 0));
}

TEST(MatAdvanced, ColumnNormalization2)
{
    auto mat = Matrix<int>(3, 2);
    mat.fill(2);

    auto soll = Matrix<double>(3, 2);
    soll.fill(2.0 / std::sqrt(12.0));

    ASSERT_TRUE(soll.compare(mat.normalizeColumns()));
}