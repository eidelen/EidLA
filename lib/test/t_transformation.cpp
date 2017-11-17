#include <gtest/gtest.h>
#include "matrix.hpp"
#include "transformation.hpp"

TEST(Transformation, EchelonForm)
{
    // Example matrix taken from
    // http://stattrek.com/matrix-algebra/echelon-transform.aspx#MatrixA

    double inData[9] = {0,1,2,  1,2,1,  2,7,8};
    Matrix<double> in(3,3,inData);

    double sollData[9] = {1,2,1,  0,1,2,  0,0,0};
    Matrix<double> soll(3,3, sollData);

    Matrix<double> res = Transformation::echelon(in);

    ASSERT_TRUE(soll.compare(res,0.01));
}

