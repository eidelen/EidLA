#include <gtest/gtest.h>
#include "matrix.hpp"

TEST(Decomposition, LUNoPivoting)
{
    // https://en.wikipedia.org/wiki/LU_decomposition

    double m_data[4] = {4.0, 3.0, 6.0, 3.0};
    auto   m         = Matrix<double>(2, 2, m_data);

    double soll_L_data[4] = {1.0, 0.0, 1.5, 1.0};
    auto   soll_L         = Matrix<double>(2, 2, soll_L_data);

    double soll_U_data[4] = {4.0, 3.0, 0.0, -1.5};
    auto   soll_U         = Matrix<double>(2, 2, soll_U_data);

    Decomposition::LUResult res           = Decomposition::luDecomposition(m,false);
    Matrix<double>          lowerTriangle = res.L;
    Matrix<double>          upperTriangle = res.U;

    ASSERT_TRUE(soll_L.compare(lowerTriangle));
    ASSERT_TRUE(soll_U.compare(upperTriangle));
}

TEST(Decomposition, LU3x3)
{
    // https://en.wikipedia.org/wiki/LU_decomposition

    double m_data[9] = {4.0, 1.0, 7.0, 3.0, 5.0, 10.0, 1.0, -4.0, 2.0};
    auto   m         = Matrix<double>(3, 3, m_data);

    Decomposition::LUResult res = Decomposition::luDecomposition(m,false);

    ASSERT_TRUE(m.compare(res.L * res.U));
}

TEST(Decomposition, LUIdent)
{
    // https://en.wikipedia.org/wiki/LU_decomposition
    auto in = Matrix<double>(6, 6);
    in.setToIdentity();

    Decomposition::LUResult res = Decomposition::luDecomposition(in,true);

    // Both should be identity
    ASSERT_TRUE(in.compare(res.L));
    ASSERT_TRUE(in.compare(res.U));
}


// example from https://math.stackexchange.com/questions/485513/what-are-pivot-numbers-in-lu-decomposition-please-explain-me-in-an-example
TEST(Decomposition, LUPivoting)
{
    double matData[] = {2,1,1,0,  4,3,3,1,  8,7,9,5,  6,7,9,8};
    auto mat = Matrix<double>(4,4,matData);

    double l_data[] = {1,0,0,0,  0.75,1,0,0,  0.5,-2.0/7.0,1,0,  0.25,-3.0/7.0,1.0/3.0,1 };
    auto l = Matrix<double>(4,4, l_data);

    double u_data[] = {8,7,9,5,  0,1.75,2.25,4.25,  0,0,-6.0/7.0,-2.0/7.0,  0,0,0,2.0/3.0};
    auto u = Matrix<double>(4,4,u_data);

    double p_data[] = {0,0,1,0,  0,0,0,1,  0,1,0,0,  1,0,0,0};
    auto p = Matrix<double>(4,4,p_data);

    Decomposition::LUResult res = Decomposition::luDecomposition(mat,true);

    ASSERT_TRUE( (p*mat).compare(l*u) );
    ASSERT_TRUE( (res.P*mat).compare(res.L*res.U) );
    ASSERT_TRUE( l.compare(res.L) );
    ASSERT_TRUE( u.compare(res.U) );
    ASSERT_TRUE( p.compare(res.P) );
    ASSERT_EQ(res.NbrRowSwaps, 3 );
}


TEST(Decomposition, LU2x2Pivoting)
{
    double matData[] = {0, -2,   1, 1};
    auto mat = Matrix<double>(2,2, matData);

    double l_data[] = {1,0, 0,1};
    auto l = Matrix<double>(2,2, l_data);
    double u_data[] = {1,1, 0,-2};
    auto u = Matrix<double>(2,2,u_data);
    double p_data[] = {0,1, 1,0};
    auto p = Matrix<double>(2,2,p_data);

    Decomposition::LUResult res = Decomposition::luDecomposition(mat,true);

    ASSERT_TRUE( (p*mat).compare(l*u) );
    ASSERT_TRUE( l.compare(res.L) );
    ASSERT_TRUE( u.compare(res.U) );
    ASSERT_TRUE( p.compare(res.P) );
}

TEST(Decomposition, EigenvalueSymmetricDiagonal)
{
    auto m = Matrix<double>::identity(3);
    for (size_t s = 0; s < 3; s++)
        m(s, s) = static_cast<double>(s + 1);

    // soll values
    auto v0 = Matrix<double>::identity(3).column(0);
    auto v1 = Matrix<double>::identity(3).column(1);
    auto v2 = Matrix<double>::identity(3).column(2);
    double l0 = 1;
    double l1 = 2;
    double l2 = 3;

    std::vector<Decomposition::EigenPair> soll;
    soll.push_back( Decomposition::EigenPair(v2,l2,true));
    soll.push_back( Decomposition::EigenPair(v1,l1,true));
    soll.push_back( Decomposition::EigenPair(v0,l0,true));

    std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m);

    ASSERT_EQ(soll.size(), eig.size());

    for( size_t i = 0; i < soll.size(); i++ )
    {
        auto cEP = eig.at(i);
        auto sEP = soll.at(i);

        ASSERT_TRUE(cEP.Valid);
        ASSERT_TRUE(cEP.V.compare(sEP.V,true,0.0001));
        ASSERT_NEAR(cEP.L, sEP.L, 0.001);
    }
}

// Example from https://en.wikipedia.org/wiki/Rayleigh_quotient_iteration
TEST(Decomposition, EigenvalueNonSymmetricLargest)
{
    double mat_data[] = {1,2,3,  1,2,1,  3,2,1};
    auto m = Matrix<double>(3,3,mat_data);

    double sollEigenVal = 3 + std::sqrt(5.0);
    double sollEigenVecData[] = {1, (1+std::sqrt(5.0))/2.0 - 1  , 1};
    auto sollEigenVec = Matrix<double>(3,1,sollEigenVecData).normalizeColumns();

    Decomposition::EigenPair sigEigenPair = Decomposition::powerIteration(m,20,std::numeric_limits<double>::epsilon());

    ASSERT_TRUE( sigEigenPair.Valid );
    ASSERT_FLOAT_EQ(sigEigenPair.L,sollEigenVal);
    ASSERT_TRUE(sollEigenVec.compare(sigEigenPair.V));
}

/*
// Example from https://en.wikipedia.org/wiki/Rayleigh_quotient_iteration
TEST(Decomposition, EigenvalueAllNonSymmetric)
{
    double mat_data[] = {1,2,3,  1,2,1,  3,2,1};
    auto m = Matrix<double>(3,3,mat_data);

    std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m);

    for (const Decomposition::EigenPair ep : eig)
    {
        if(ep.Valid)
            std::cout << "Valid" << std::endl;
        else
            std::cout << "Invalid" << std::endl;

        std::cout << "Eigenvalue: " << ep.L << std::endl;
        std::cout << "Eigenvector: " << ep.V.transpose() ;

        if( (m*ep.V).compare(ep.L*ep.V) )
        {
            std::cout << "True" << std::endl;
        }
        else
        {
            auto eucDist = ((m*ep.V) - (ep.L*ep.V));
            std::cout << "False: err = " << eucDist.transpose() << " length: " <<std::sqrt( (eucDist.transpose() * eucDist)(0,0)) << std::endl;
        }

        std::cout << "------" << std::endl << std::endl;
    }
}
*/

TEST(Decomposition, QRDecomposition)
{
    double matData[] = {0.8147, 0.0975, 0.1576,
                        0.9058, 0.2785, 0.9706,
                        0.1270, 0.5469, 0.9572,
                        0.9134, 0.9575, 0.4854,
                        0.6324, 0.9649, 0.8003 };

    auto mat = Matrix<double>(5, 3, matData);


}