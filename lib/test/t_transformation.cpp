#include <gtest/gtest.h>
#include "matrix.hpp"
#include "transformation.hpp"

TEST(Transformation, EchelonForm)
{
    // Example matrix taken from
    // http://stattrek.com/matrix-algebra/echelon-transform.aspx#MatrixA

    double         inData[9] = {0, 1, 2, 1, 2, 1, 2, 7, 8};
    Matrix<double> in(3, 3, inData);

    double         sollData[9] = {1, 2, 1, 0, 1, 2, 0, 0, 0};
    Matrix<double> soll(3, 3, sollData);

    Matrix<double> res = Transformation::echelon(in);

    ASSERT_TRUE(soll.compare(res, 0.01));
}

TEST(Transformation, EchelonFormIdent)
{
    auto in  = Matrix<double>::identity(4);
    auto ech = Transformation::echelon(in);
    ASSERT_TRUE(in.compare(ech, std::numeric_limits<double>::min()));
}

TEST(Transformation, Echelon0)
{
    auto in  = Matrix<double>(5,5);
    in.fill(0.0);
    auto ech = Transformation::echelon(in);
    ASSERT_TRUE(ech.compare(in, std::numeric_limits<double>::min()));
}

TEST(Transformation, ReducedEchelon)
{
    // Example matrix taken from
    // http://stattrek.com/matrix-algebra/echelon-transform.aspx#MatrixA

    double         inData[9] = {0, 1, 2, 1, 2, 1, 2, 7, 8};
    Matrix<double> in(3, 3, inData);

    double         sollData[9] = {1,0,-3,  0,1,2, 0,0,0};
    Matrix<double> soll(3, 3, sollData);

    Matrix<double> res = Transformation::reduced_echelon(in);

    ASSERT_TRUE(soll.compare(res, 0.01));
}


TEST(Transformation, ReducedEchelon2)
{
    //https://en.wikipedia.org/wiki/Gaussian_elimination
    double inData[] = {2,1,-1,8,  -3,-1,2,-11,  -2,1,2,-3};
    auto mat = Matrix<double>(3,4,inData);

    double sollReducedEchelon[] = {1,0,0,2,   0,1,0,3,  0,0,1,-1};
    auto soll = Matrix<double>(3,4,sollReducedEchelon);

    auto redEch = Transformation::reduced_echelon(mat);

    ASSERT_TRUE(soll.compare(redEch));
}

TEST(Transformation, MatrixRank)
{
    double         inData[9] = {0, 1, 2, 1, 2, 1, 2, 7, 8};
    Matrix<double> in(3, 3, inData);

    ASSERT_EQ(in.getRank(), 2);
}

TEST(Transformation, MatrixRankIdent) {
    auto in = Matrix<double>::identity(3);
    ASSERT_EQ(in.getRank(), 3);

    in.setRow(1, in.row(2)); // set row 8 identical with row 9 -> linear dependant

    ASSERT_EQ(in.getRank(), 2);
}

TEST(Transformation, MatrixRank0) {
    auto in = Matrix<double>(5,5);
    in.fill(0.0);
    ASSERT_EQ(in.getRank(), 0);
}

/*

//https://math.dartmouth.edu/archive/m23s06/public_html/handouts/row_reduction_examples.pdf
TEST(Transformation, GaussElimination)
{
    double inData[] = {0, 2, 1, -8,   1,-2,-3,0,    -1, 1, 2, 3};
    auto in = Matrix<double>(3,4,inData);

    std::cout << "Gauss" << Transformation::echelon(in);



}*/