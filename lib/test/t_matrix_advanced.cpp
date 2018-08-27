#include <gtest/gtest.h>
#include "matrix.hpp"

// http://stattrek.com/matrix-algebra/how-to-find-inverse.aspx
TEST(MatAdvanced, Inverse)
{
    double inData[] = {1, 2, 2, 2, 2, 2, 2, 2, 1};
    auto   mat      = Matrix<double>(3, 3, inData);

    double sollData[] = {-1, 1, 0, 1, -1.5, 1, 0, 1, -1};
    auto   soll       = Matrix<double>(3, 3, sollData);

    auto computedInverse = mat.inverted();

    ASSERT_TRUE(soll.compare(computedInverse));

    // inv(mat) * mat = identity
    auto identRes  = computedInverse * mat;
    auto sollIdent = Matrix<double>::identity(3);
    ASSERT_TRUE(sollIdent.compare(identRes));
}

TEST(MatAdvanced, NonInvertable)
{
    double inData[] = {1, 2, 2, 2, 2, 2};
    auto   mat      = Matrix<double>(3, 2, inData);

    EXPECT_ANY_THROW(auto computedInverse = mat.inverted());
}

TEST(MatAdvanced, InvertSingularMatrix)
{
    auto mat = Matrix<int>(3, 3);
    mat.fill(2);

    EXPECT_ANY_THROW(auto computedInverse = mat.inverted());
}

// checked with octave
TEST(MatAdvanced, Determinant1)
{
    double inData[] = {1, 2, 2, 2, 0, -1, -2, 1, 3};
    auto   mat      = Matrix<double>(3, 3, inData);

    double det = mat.determinant();

    ASSERT_NEAR(det, -3.0, 0.00001);
}

// example from https://www.mathsisfun.com/algebra/matrix-determinant.html
TEST(MatAdvanced, Determinant2)
{
    double inData[] = {6, 1, 1, 4, -2, 5, 2, 8, 7};
    auto   mat      = Matrix<double>(3, 3, inData);

    double det = mat.determinant();
    ASSERT_NEAR(det, -306.0, 0.00001);
}

TEST(MatAdvanced, DeterminantIdentity)
{
    auto mat = Matrix<double>::identity(10);
    ASSERT_NEAR(mat.determinant(), 1.0, 0.00001);
}

TEST(MatAdvanced, DeterminantSingularMatrix)
{
    auto mat = Matrix<int>(3, 3);
    mat.fill(2);
    ASSERT_NEAR(mat.determinant(), 0.0, 0.00001);
}

TEST(MatAdvanced, Determinant2x2)
{
    double matData[] = {0, -2,   1, 1};
    auto mat = Matrix<double>(2,2, matData);

    ASSERT_NEAR(mat.determinant(), 2.0, 0.00001);
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

// example from https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html
TEST(MatAdvanced, FirstMinors)
{
    double matData[] = {3,0,2,  2,0,-2,  0,1,1};
    auto mat = Matrix<double>(3,3,matData);

    double sollData[] = {2,2,2,  -2,3,3,  0,-10,0};
    auto soll = Matrix<double>(3,3,sollData);

    auto firstMinors = mat.firstMinors();

    ASSERT_TRUE( soll.compare(firstMinors) );
}

// example from https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html
TEST(MatAdvanced, Cofactors)
{
    double matData[] = {3,0,2,  2,0,-2,  0,1,1};
    auto mat = Matrix<double>(3,3,matData);

    double sollData[] = {2,-2,2,  2,3,-3,  0,10,0};
    auto soll = Matrix<double>(3,3,sollData);

    auto cf = mat.cofactors();

    ASSERT_TRUE( soll.compare(cf) );
}

// example from https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html
TEST(MatAdvanced, Adjugate)
{
    double matData[] = {3,0,2,  2,0,-2,  0,1,1};
    auto mat = Matrix<double>(3,3,matData);

    // adjugate is the transpose of cofactors
    double cfData[] = {2,-2,2,  2,3,-3,  0,10,0};
    auto cf = Matrix<double>(3,3,cfData);
    auto soll = cf.transpose();

    auto aj = mat.adjugate();

    ASSERT_TRUE( soll.compare(aj) );
}

TEST(MatAdvanced, IsOrthogonal)
{
    Matrix<double> mat = Matrix<double>::identity(4);
    ASSERT_TRUE(mat.isOrthogonal());

    mat(0,1) = 1.0;
    ASSERT_FALSE(mat.isOrthogonal());
}

TEST(MatAdvanced, L1NormMatrix)
{
    double matData[] = {1,-2,-3,  -4,5,6,  7,-8,9};
    auto mat = Matrix<double>(3,3,matData);

    double l1NormShould = 3 + 6 + 9; // third column abs sum

    ASSERT_EQ( l1NormShould, mat.normL1() );
}

TEST(MatAdvanced, L1NormVector)
{
    double matData[] = {1,-2,-3};
    auto mat = Matrix<double>(3,1,matData);

    double l1NormShould = 1 + 2 + 3;

    ASSERT_EQ( l1NormShould, mat.normL1() );
}

TEST(MatAdvanced, InfNormMatrix)
{
    double matData[] = {1,-2,-3,  -4,5,6,  7,-8,9};
    auto mat = Matrix<double>(3,3,matData);

    double infNormShould = 7 + 8 + 9;

    ASSERT_EQ( infNormShould, mat.normInf() );
}

TEST(MatAdvanced, InfNormVector)
{
    double matData[] = {1,-2,-3};
    auto mat = Matrix<double>(3,1,matData);

    double infNormShould = 3;

    ASSERT_EQ( infNormShould, mat.normInf() );
}

// source https://ch.mathworks.com/help/matlab/ref/norm.html
TEST(MatAdvanced, L2NormMatrix)
{
    double matData[] = {2, 0, 1,    -1, 1, 0,   -3, 3, 0};
    auto mat = Matrix<double>(3,3,matData);

    double l2NormShould = 4.7234;

    ASSERT_NEAR( l2NormShould, mat.normL2(), 0.001 );
}

// from documents/norms.pdf
TEST(MatAdvanced, L2NormVector)
{
    double matData[] = {3,-4};
    auto mat = Matrix<double>(2,1,matData);

    double l2Norm = 5;

    ASSERT_EQ( l2Norm, mat.normL2() );
}

// from documents/norms.pdf
TEST(MatAdvanced, ConditionNumbers)
{
    double matData[] = {2, -1, 1,   1, 0, 1,   3, -1, 4};
    auto mat = Matrix<double>(3,3, matData);

    ASSERT_NEAR( 6 * 4.5, mat.conditionNumberL1(), 0.001);
    ASSERT_NEAR( 8 * 3.5, mat.conditionNumberInf(), 0.001);
}
