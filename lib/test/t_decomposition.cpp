#include <gtest/gtest.h>
#include "matrix.hpp"

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
    Matrix<double> lowerTriangle = res.L;
    Matrix<double> upperTriangle = res.U;

    ASSERT_TRUE(soll_L.compare(lowerTriangle));
    ASSERT_TRUE(soll_U.compare(upperTriangle));
}

TEST(Decomposition, LU3x3)
{
    // https://en.wikipedia.org/wiki/LU_decomposition

    double m_data[9] = {4.0, 1.0, 7.0,     3.0, 5.0, 10.0,     1.0, -4.0, 2.0};
    auto m = Matrix<double>(3, 3, m_data);

    Decomposition::LUResult res = Decomposition::luDecomposition(m);

    ASSERT_TRUE(m.compare(res.L*res.U));
}

TEST(Decomposition, LUIdent)
{
    // https://en.wikipedia.org/wiki/LU_decomposition
    auto in = Matrix<double>(6,6);
    in.setToIdentity();

    Decomposition::LUResult res = Decomposition::luDecomposition(in);

    // Both should be identity
    ASSERT_TRUE(in.compare(res.L));
    ASSERT_TRUE(in.compare(res.U));
}

TEST(Decomposition, Eigenvalue)
{
    double in_data[] = {1,2,  2,4};
    auto in = Matrix<double>(2,2,in_data);

    std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(in);

    ASSERT_GE(eig.size(), 1);
    Decomposition::EigenPair fEig = eig.at(0);

    ASSERT_DOUBLE_EQ(5.0, fEig.L);

    double soll_data[] = {1,2};
    auto soll = Matrix<double>(2,1,soll_data);
    auto sollN = soll.normalizeColumns();

    std::cout << sollN << fEig.V ;

    // Depending on initialization, Eigenvector can point in opposite direction too.
    ASSERT_TRUE( sollN.compare(fEig.V)  ||  sollN.compare( (fEig.V)*(-1)) );

}