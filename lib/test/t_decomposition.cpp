#include <gtest/gtest.h>
#include "matrix.hpp"
#include "decomposition.hpp"

TEST(Decomposition, LU)
{
    // https://en.wikipedia.org/wiki/LU_decomposition

    double m_data[4] = {4.0, 3.0, 6.0, 3.0};
    auto m = Matrix<double>(2, 2, m_data);

    double soll_L_data[4] = {1.0, 0.0, 1.5, 1.0};
    auto soll_L = Matrix<double>(2, 2, soll_L_data);

    double soll_U_data[4] = {4.0, 3.0, 0.0, -1.5};
    auto soll_U = Matrix<double>(2, 2, soll_U_data);

    Decomposition::LUResult res = Decomposition::luDecomposition(m);

    ASSERT_TRUE(soll_L.compare(res.L, 0.000001));
    ASSERT_TRUE(soll_U.compare(res.U, 0.000001));
}