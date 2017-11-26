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
    ASSERT_TRUE(in.compare(ech, 0.0000001));
}

TEST(Transformation, MatrixRank)
{
    double         inData[9] = {0, 1, 2, 1, 2, 1, 2, 7, 8};
    Matrix<double> in(3, 3, inData);

    ASSERT_EQ(in.getRank(), 2);
}

TEST(Transformation, MatrixRankIdent)
{
    auto in = Matrix<double>::identity(3);
    ASSERT_EQ(in.getRank(), 3);

    in.setRow(1, in.row(2)); // set row 8 identical with row 9 -> linear dependant

    ASSERT_EQ(in.getRank(), 2);
}
