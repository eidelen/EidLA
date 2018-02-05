#include <gtest/gtest.h>
#include "matrix.hpp"
#include "solve.hpp"


// https://www.mathsisfun.com/algebra/systems-linear-equations-matrices.html
TEST(Solve, SimpleLinearSystem)
{
    double c_data[] = {1, 1, 1,
                       0, 2, 5,
                       2, 5, -1};
    auto c = Matrix<double>(3, 3, c_data);

    double b_data[] = {6, -4, 27};
    auto b = Matrix<double>(3, 1, b_data);

    double x_data[] = {5, 3, -2};
    auto x_soll = Matrix<double>(3, 1, x_data);

    // check that example is correct
    ASSERT_TRUE( b.compare( c * x_soll ));

    auto x_comp = Solve::solve_lseq(c,b);

    ASSERT_TRUE( x_soll.compare( x_comp ));
}

TEST(Solve, BatchLinearSystem)
{
    for( int i = 0; i < 1000; i++ )
    {
        auto c = Matrix<double>::random(10,10, -100, +100);
        auto b = Matrix<double>::random(10,1, -100, +100);

        auto x = Solve::solve_lseq(c,b);

        ASSERT_TRUE( b.compare(c*x, true, 0.000001) );
    }
}


