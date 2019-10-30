#include <gtest/gtest.h>
#include <inc/matrix4x4.hpp>
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

TEST(Decomposition, EigenvalueSymmetricDiagonalQR)
{
    auto m = Matrix<double>::identity(3);
    for (size_t s = 0; s < 3; s++)
        m(s, s) = static_cast<double>(s + 1);

    std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m, Decomposition::QRAlgorithm);

    ASSERT_EQ(3, eig.size());

    for( size_t i = 0; i < eig.size(); i++ )
    {
        auto cEP = eig.at(i);

        ASSERT_TRUE( cEP.Valid );
        ASSERT_TRUE( (m * cEP.V).compare( cEP.L * cEP.V, true, 0.00001 ) );
    }
}

TEST(Decomposition, EigenvalueSymmetricDiagonalPower)
{
    auto m = Matrix<double>::identity(3);
    for (size_t s = 0; s < 3; s++)
        m(s, s) = static_cast<double>(s + 1);

    std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m, Decomposition::PowerIterationAndHotellingsDeflation);

    ASSERT_EQ(3, eig.size());

    for( size_t i = 0; i < eig.size(); i++ )
    {
        auto cEP = eig.at(i);

        ASSERT_TRUE( cEP.Valid );
        ASSERT_TRUE( (m * cEP.V).compare( cEP.L * cEP.V, true, 0.001 ) );
    }
}


TEST(Decomposition, EigenvalueSymmetric)
{
    double data[] = { 6,10,11,  10,17,21,  11,21,42};
    auto m = Matrix<double>(3,3,data);

    // qr is called
    std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m);

    ASSERT_EQ(3, eig.size());

    // just check if properties are fullfilled.
    for( size_t i = 0; i < eig.size(); i++ )
    {
        auto cEP = eig.at(i);

        ASSERT_TRUE( cEP.Valid );
        ASSERT_TRUE( (m * cEP.V).compare( cEP.L * cEP.V, true, 0.00001 ) );
    }
}

TEST(Decomposition, EigenvalueSymmetricBatchTesting)
{
    for( int k = 0; k < 100; k++ )
    {
        auto m = Matrix<double>::random(5, 5, -10.0, 10.0);
        m = m * m.transpose();

        std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m);
        for( size_t i = 0; i < eig.size(); i++ )
        {
            auto cEP = eig.at(i);

            if( cEP.Valid )
            {
                ASSERT_TRUE((m * cEP.V).compare(cEP.L * cEP.V, true, 0.00001));
            }
        }
    }
}

TEST(Decomposition, Eigenvalue2x2)
{
    // test different cases

    // m01 = 0;
    for( int k = 0; k < 1000; k++ )
    {
        auto m = Matrix<double>::random(2, 2, -10.0, 10.0);
        m(0,1) = 0.0;

        std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m);
        for (size_t i = 0; i < eig.size(); i++)
        {
            auto cEP = eig.at(i);
            ASSERT_TRUE((m * cEP.V).compare(cEP.L * cEP.V, true, 0.00001));
        }
    }

    // m10 = 0;
    for( int k = 0; k < 1000; k++ )
    {
        auto m = Matrix<double>::random(2, 2, -10.0, 10.0);
        m(1,0) = 0.0;

        std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m);
        for (size_t i = 0; i < eig.size(); i++)
        {
            auto cEP = eig.at(i);
            ASSERT_TRUE((m * cEP.V).compare(cEP.L * cEP.V, true, 0.00001));
        }
    }

    // m10 = 0 & m01 = 0;
    for( int k = 0; k < 1000; k++ )
    {
        auto m = Matrix<double>::random(2, 2, -10.0, 10.0);
        m(1,0) = 0.0;
        m(0,1) = 0.0;

        std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m);
        for (size_t i = 0; i < eig.size(); i++)
        {
            auto cEP = eig.at(i);
            ASSERT_TRUE((m * cEP.V).compare(cEP.L * cEP.V, true, 0.00001));
        }
    }
}

TEST(Decomposition, Eigenvalue2x2BatchTesting)
{
    for( int k = 0; k < 500; k++ )
    {
        auto m = Matrix<double>::random(2, 2, -10.0, 10.0);
        m = m * m.transpose();

        std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(m);
        for( size_t i = 0; i < eig.size(); i++ )
        {
            auto cEP = eig.at(i);

            if( cEP.Valid )
            {
                ASSERT_TRUE((m * cEP.V).compare(cEP.L * cEP.V, true, 0.00001));
            }
        }
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

    Decomposition::EigenPair sigEigenPair = Decomposition::powerIteration(m,40,std::numeric_limits<double>::epsilon());

    ASSERT_TRUE( sigEigenPair.Valid );
    ASSERT_FLOAT_EQ(sigEigenPair.L,sollEigenVal);
    ASSERT_TRUE(sollEigenVec.compare(sigEigenPair.V));


    // This calls also the power iteration
    std::vector<Decomposition::EigenPair> eigenPair = Decomposition::eigen(m);

    ASSERT_TRUE( eigenPair.at(0).Valid );
    ASSERT_FLOAT_EQ(eigenPair.at(0).L,sollEigenVal);
    ASSERT_TRUE(sollEigenVec.compare(eigenPair.at(0).V));
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

// example from document /documents/qr_decomposition.pdf
TEST(Decomposition, QRDecompositionHouseholder)
{
    double matData[] = {0.8147, 0.0975, 0.1576,
                        0.9058, 0.2785, 0.9706,
                        0.1270, 0.5469, 0.9572,
                        0.9134, 0.9575, 0.4854,
                        0.6324, 0.9649, 0.8003 };

    auto mat = Matrix<double>(5, 3, matData);

    double rData[] = {1.6536,1.1405,1.2569,
                      0, 0.9661, 0.6341,
                      0, 0, 0.8816,
                      0, 0, 0,
                      0, 0, 0};
    auto rSoll = Matrix<double>(5,3,rData);

    Decomposition::QRResult res = Decomposition::qr_householder(mat);

    ASSERT_TRUE( res.Q.isOrthogonal(0.01)); // directly comparing Q does not make sense, since signs can change
    ASSERT_TRUE( rSoll.compare(res.R, true, 0.01));
    ASSERT_TRUE( mat.compare( res.Q * res.R, true, 0.01));
}

// example from document /documents/qr_decomposition.pdf
TEST(Decomposition, QRDecompositionGivens)
{
    double matData[] = {0.8147, 0.0975, 0.1576,
                        0.9058, 0.2785, 0.9706,
                        0.1270, 0.5469, 0.9572,
                        0.9134, 0.9575, 0.4854,
                        0.6324, 0.9649, 0.8003 };

    auto mat = Matrix<double>(5, 3, matData);

    Decomposition::QRResult res = Decomposition::qr_givens(mat);

    ASSERT_TRUE( res.Q.isOrthogonal(0.001)); // directly comparing Q does not make sense, since signs can change
    ASSERT_TRUE( mat.compare( res.Q * res.R, true, 0.001));
}

TEST(Decomposition, QRDecompositionSpecialMatrixWhichFailed)
{
    double matData[] = {2.0/3.0, -1.0/3.0, -1.0/3.0,
                        -1.0/3.0, 2.0/3.0, -1.0/3.0,
                        -1.0/3.0, -1.0/3.0, 2.0/3.0};
    Matrix<double> mat(3,3, matData);
    Decomposition::QRResult res = Decomposition::qr(mat);
}

TEST(Decomposition, QRDecompositionSpecialMatrix)
{
    // Rotation Matrix
    Matrix4x4 rotMat;
    rotMat.rotZ(0.2); rotMat.rotX(0.4);

    Decomposition::QRResult res = Decomposition::qr(rotMat);

    ASSERT_TRUE( rotMat.compare( res.Q * res.R, true, 0.0001 ));

    // check if Q orthogonal
    for( size_t i = 0; i < 4; i++ )
    {
        ASSERT_FLOAT_EQ( res.Q.column(i).norm(), 1.0 );
        ASSERT_FLOAT_EQ( res.Q.row(i).norm(), 1.0 );
    }
}

TEST(Decomposition, QRBatchTesting)
{
    for( int k = 0; k < 5000; k++ )
    {
        auto m = Matrix<double>::random(5, 5, -100.0, 100.0);

        // check householder qr algorithm
        Decomposition::QRResult res_h = Decomposition::qr(m,true,Decomposition::Householder);
        ASSERT_TRUE( res_h.Q.isOrthogonal(0.00001) );
        ASSERT_TRUE( m.compare(res_h.Q*res_h.R, true, 0.00001) );

        // check givens qr algorithm
        Decomposition::QRResult res_g = Decomposition::qr(m,true,Decomposition::Givens);
        ASSERT_TRUE( res_g.Q.isOrthogonal(0.0001) );
        ASSERT_TRUE( m.compare(res_g.Q*res_g.R, true, 0.00001) );
    }
}

TEST(Decomposition, QRSignChanger)
{
    double matData[] = {0.8147, 0.0975, 0.1576,
                        0.9058, 0.2785, 0.9706,
                        0.1270, 0.5469, 0.9572,
                        0.9134, 0.9575, 0.4854,
                        0.6324, 0.9649, 0.8003 };

    auto mat = Matrix<double>(5, 3, matData);

    Decomposition::QRResult res = Decomposition::qr(mat);

    for( size_t k = 0; k < res.R.rows(); k++ )
    {
        Decomposition::QRResult nextRes = Decomposition::qrSignModifier(res.Q, res.R, k);
        ASSERT_TRUE(mat.compare(nextRes.Q * nextRes.R, true, 0.001));
    }
}

TEST(Decomposition, RQBatchTesting)
{
    for( int k = 0; k < 500; k++ )
    {
        auto m = Matrix<double>::random(6, 6, -10.0, 10.0);

        Decomposition::QRResult res = Decomposition::rq(m);

        auto R = res.R;
        auto Q = res.Q;

        ASSERT_TRUE( Q.isOrthogonal(0.001) );
        ASSERT_TRUE( m.compare(R*Q, true, 0.001) );
    }
}


TEST(Decomposition, HouseHolderBatchTesting)
{
    size_t rows = 5;

    for( int k = 0; k < 100; k++ )
    {
        auto x = Matrix<double>::random(rows, 1, -100.0, 100.0);

        Decomposition::HouseholderResult res = Decomposition::householder(x);
        auto p = Decomposition::householderMatrix(res.V, res.B);
        ASSERT_TRUE( p.isOrthogonal(0.0001) );
        ASSERT_TRUE( p.isSymmetric() );

        // rx should be a multiple of e1 with length x -> rx = [ +/- norm(x), 0, 0, ... 0 ]
        auto rx = p * x;

        auto rxSoll = Matrix<double>(rows, 1);
        rxSoll.fill(0.0);
        rxSoll(0,0) = x.norm();

        ASSERT_TRUE(rxSoll.compare(rx, true, 0.00001) || rxSoll.compare(rx*(-1.0),true, 0.00001) );
    }
}

bool testDiagonalizationResult( const Matrix<double>& a, Decomposition::DiagonalizationResult res)
{
    bool ok = true;

    ok = ok && res.V.isOrthogonal(0.0001) ;
    ok = ok && res.U.isOrthogonal(0.0001) ;
    ok = ok && res.D.compare(  res.U.transpose() * a * res.V, true, 0.0001 );
    ok = ok && a.compare(  res.U * res.D * res.V.transpose(), true, 0.0001 );

    return ok;
}

TEST(Decomposition, BidiagonalizationSymmetric)
{
/* Octave example
A =
   6   8   8   7
   7   8   6   5
   4   2   3   7
   9   1   4   8

U =
   0.4447496  -0.6879265  -0.0015212   0.5735439
   0.5188745  -0.1063169   0.6637464  -0.5281161
   0.2964997  -0.2746100  -0.7223026  -0.5612097
   0.6671244   0.6633575  -0.1942097   0.2778209

B =
   13.49074  -18.90026    0.00000    0.00000
    0.00000    6.65700   -4.13058    0.00000
    0.00000    0.00000    5.03030    0.19847
    0.00000    0.00000    0.00000    1.02931

V =

   1.00000   0.00000   0.00000   0.00000
   0.00000  -0.47455   0.74583  -0.46748
   0.00000  -0.54122   0.17160   0.82318
   0.00000  -0.69418  -0.64365  -0.32223
*/

    // Compared to the example, signs might change
    double aData[] = {6,8,8,7,  7,8,6,5,   4,2,3,7,   9,1,4,8};
    auto a = Matrix<double>(4,4,aData);
    Decomposition::DiagonalizationResult res = Decomposition::bidiagonalization(a);

    ASSERT_TRUE(testDiagonalizationResult(a, res));
}

TEST(Decomposition, BidiagonalizationNonSymmetric)
{
/* Octave example
A =
           1           4           1           6
          10           0           6           4
           7          10           2           7
           8           4           5           3
           7          10           0           1


U =
   -0.061663      0.6189     0.44282     0.16385     0.62468
    -0.61663    -0.50111     0.26481     0.53573     0.10737
    -0.43164     0.58508   -0.019502     0.19214    -0.65884
     -0.4933   -0.078335     0.32425     -0.8033   0.0097606
    -0.43164      0.1319    -0.79263    -0.06282     0.40507

B =
      16.217     -15.114           0           0
           0      8.8768     -6.5868           0
           0           0      4.7876      5.4597
           0           0           0      2.3784
           0           0           0           0

V =
          -1          -0          -0          -0
          -0     0.71806    -0.44909    -0.53171
          -0     0.46918     0.87662    -0.10678
          -0     0.51406    -0.17279     0.84017
*/

    // Compared to the example, signs might change
    double aData[] = {1,4,1,6,    10,0,6,4,    7,10,2,7,    8,4,5,3,     7,10,0,1 };
    auto a = Matrix<double>(4,4,aData);
    Decomposition::DiagonalizationResult res = Decomposition::bidiagonalization(a);

    ASSERT_TRUE(testDiagonalizationResult(a, res));
}


TEST(Decomposition, BidiagonalizationBatchTest)
{
    size_t rows = 6;

    for( int k = 0; k < 200; k++ )
    {
        auto a = Matrix<double>::random(rows, rows, -100.0, 100.0);

        Decomposition::DiagonalizationResult res = Decomposition::bidiagonalization(a);
        ASSERT_TRUE(testDiagonalizationResult(a, res));
    }
}

// http://mysite.science.uottawa.ca/phofstra/MAT2342/SVDproblems.pdf
TEST(Decomposition, SVD)
{
    double matData[] = {0.0, 1.0, 1.0,
                        std::sqrt(2.0), 2.0, 0,
                        0.0, 1.0, 1.0};

    auto mat = Matrix<double>(3, 3, matData);

    Matrix<double> s_soll = Matrix<double>(3,3);
    s_soll.fill(0.0);
    s_soll(0,0) = 2.0 * std::sqrt(2.0);
    s_soll(1,1) = std::sqrt(2.0);

    /*
    double u_soll_data[] = {1.0/std::sqrt(6.0),   -1.0/std::sqrt(3.0),   1/std::sqrt(2.0),
                       2.0/std::sqrt(6.0),   1.0/std::sqrt(3.0),       0.0,
                       1.0/std::sqrt(6.0),   -1.0/std::sqrt(3.0),   -1.0/std::sqrt(2.0)};
    auto u_soll = Matrix<double>(3,3,u_soll_data);

    double v_soll_data[] = {1.0/std::sqrt(6.0),   1.0/std::sqrt(3.0),    1/std::sqrt(2.0),
                            3.0/std::sqrt(12.0),   0.0,           -0.5,
                            1.0/std::sqrt(12.0),   -2.0/std::sqrt(6.0),   0.5};
    auto v_soll = Matrix<double>(3,3,v_soll_data); */


    Decomposition::SVDResult res = Decomposition::svd(mat);

    // U and V cannot be compared directly, because it can change in signs.
    ASSERT_TRUE( res.U.isOrthogonal(0.001) );
    ASSERT_TRUE( res.V.isOrthogonal(0.001) );

    std::cout << res.S << std::endl;
    std::cout << s_soll << std::endl;

    ASSERT_TRUE( s_soll.compare(res.S, true, 0.001) );

    ASSERT_TRUE( mat.compare(res.U * res.S * res.V.transpose(), true, 0.001 ) );
}

TEST(Decomposition, SVDIdentity)
{
    Matrix<double> ident = Matrix<double>::identity(3);
    Decomposition::SVDResult res = Decomposition::svd(ident);
    ASSERT_TRUE( res.S.compare(ident, true, 0.001) );
}

TEST(Decomposition, SVDProblemMatrix)
{
    // this is from 3d3d registration. svd does not find the right rotation matrix.
    double h_data[] = {2, 1, -1, 1, 2, 0, -1, 0, 2.0/3.0};
    Matrix<double> h = Matrix<double>(3,3,h_data);

    Decomposition::SVDResult res = Decomposition::svd(h);

    std::cout << "U ortho = " << res.U.isOrthogonal(0.01) << std::endl;
    std::cout << "V ortho = " << res.V.isOrthogonal(0.01) << std::endl;
    std::cout << "SVD correct = " << h.compare(res.U * res.S * res.V.transpose(), true, 0.0001 ) << std::endl;

    std::cout << "dec.V = " << std::endl << res.V<< std::endl;
    std::cout << "dec.S = " << std::endl << res.S<< std::endl;
    std::cout << "dec.U = " << std::endl << res.U<< std::endl;

    // rotation should be identity
    Matrix<double> r1 = res.V * res.U.transpose();

    double d1 = r1.determinant();
    std::cout << "d1 = " << d1 << std::endl;
    std::cout << "r1 = " << std::endl << r1 << std::endl;

    Matrix<double> vModified = res.V;
    vModified.setColumn(2, vModified.column(2) * (-1.0) );
    Matrix<double> r2 = vModified * res.U.transpose();
    double d2 = r2.determinant();
    std::cout << "d2 = " << d2 << std::endl;
    std::cout << "r2 = " << std::endl << r2 << std::endl;
}

TEST(Decomposition, SVDBatch)
{
    size_t rows = 5;

    for( int k = 0; k < 100; k++ )
    {
        auto a = Matrix<double>::random(rows, rows, -1.0, 1.0);

        Decomposition::SVDResult res = Decomposition::svdEigen(a);

        ASSERT_TRUE( res.U.isOrthogonal(0.05) );
        ASSERT_TRUE( res.V.isOrthogonal(0.05) );
        ASSERT_TRUE( a.compare(res.U * res.S * res.V.transpose(), true, 0.1 ) );
    }
}

TEST(Decomposition, SVDCompression)
{
    // third line is almost linear combination of both upper
    // therefore, using first two singular values should lead to a accurate representation of the matrix.
    double data[] = {1, 0, 0,   0, 1, 0,   1, 1, 0.01}; // third line is almost linear combination of both upper
    Matrix<double> mat(3,3, data);

    Decomposition::SVDResult res = Decomposition::svd(mat);

    ASSERT_TRUE( res.U.isOrthogonal(0.05) );
    ASSERT_TRUE( res.V.isOrthogonal(0.05) );
    ASSERT_TRUE( mat.compare(res.U * res.S * res.V.transpose(), true, 0.05 ) );

    auto compressed = res.U.subMatrix(0,0,3,2) * res.S.subMatrix(0,0,2,2) * res.V.transpose().subMatrix(0,0,2,3);

    ASSERT_TRUE( mat.compare(compressed, true, 0.01));
}

TEST(Decomposition, SortEigenPair)
{
    Matrix<double> ev1 = Matrix<double>::random(3,1,0.0,10.0);
    Decomposition::EigenPair ep1(ev1, 1.0, true);

    Matrix<double> ev2 = Matrix<double>::random(3,1,0.0,10.0);
    Decomposition::EigenPair ep2(ev2, 2.0, true);

    Matrix<double> ev3 = Matrix<double>::random(3,1,0.0,10.0);
    Decomposition::EigenPair ep3(ev3, 3.0, true);

    std::vector<Decomposition::EigenPair> vec;
    vec.push_back(ep3);
    vec.push_back(ep1);
    vec.push_back(ep2);

    Decomposition::sortDescending(vec);

    ASSERT_FLOAT_EQ( ep3.L, vec.at(0).L );
    ASSERT_FLOAT_EQ( ep2.L, vec.at(1).L );
    ASSERT_FLOAT_EQ( ep1.L, vec.at(2).L );

    ASSERT_TRUE( ep1.V.compare(vec.at(2).V)) ;
    ASSERT_TRUE( ep2.V.compare(vec.at(1).V)) ;
    ASSERT_TRUE( ep3.V.compare(vec.at(0).V)) ;

    Decomposition::sortAscending(vec);

    ASSERT_FLOAT_EQ( ep3.L, vec.at(2).L );
    ASSERT_FLOAT_EQ( ep2.L, vec.at(1).L );
    ASSERT_FLOAT_EQ( ep1.L, vec.at(0).L );

    ASSERT_TRUE( ep1.V.compare(vec.at(0).V)) ;
    ASSERT_TRUE( ep2.V.compare(vec.at(1).V)) ;
    ASSERT_TRUE( ep3.V.compare(vec.at(2).V)) ;
}

TEST(Decomposition, FilterEigenPair)
{
    Matrix<double> ev1 = Matrix<double>::random(3,1,0.0,10.0);
    Decomposition::EigenPair ep1(ev1, 1.0, true);

    Matrix<double> ev2 = Matrix<double>::random(3,1,0.0,10.0);
    Decomposition::EigenPair ep2(ev2, 2.0, true);

    Matrix<double> ev3 = Matrix<double>::random(3,1,0.0,10.0);
    Decomposition::EigenPair ep3(ev3, 3.0, true);

    std::vector<Decomposition::EigenPair> vec;
    vec.push_back(ep3);
    vec.push_back(ep1);
    vec.push_back(ep2);

    Decomposition::filter(vec, 2.1);

    ASSERT_FLOAT_EQ( ep3.L, vec.at(0).L );
    ASSERT_TRUE( ep3.V.compare(vec.at(0).V)) ;

    ASSERT_EQ(1, vec.size());
}

TEST(Decomposition, GivensRotation)
{
    {
        Matrix<double> d1(2, 1);
        d1(0, 0) = 1;
        d1(1, 0) = 1;

        Decomposition::GivensRotation r1 = Decomposition::givensRotation(d1(0, 0), d1(1, 0));

        Matrix<double> g1(2, 2);
        g1(0, 0) = r1.C;
        g1(1, 1) = r1.C;
        g1(0, 1) = -r1.S;
        g1(1, 0) = r1.S;

        Matrix<double> n1 = g1 * d1;

        ASSERT_NEAR(0.0, n1(1, 0), 0.00000001);
        ASSERT_NEAR(std::sqrt(2), n1(0, 0), 0.000001);
        ASSERT_NEAR(std::sqrt(2), r1.R, 0.000001);
    }

    {
        Matrix<double> d1(2, 1);
        d1(0, 0) = -1;
        d1(1, 0) = -1;

        Decomposition::GivensRotation r1 = Decomposition::givensRotation(d1(0, 0), d1(1, 0));

        Matrix<double> g1(2, 2);
        g1(0, 0) = r1.C;
        g1(1, 1) = r1.C;
        g1(0, 1) = -r1.S;
        g1(1, 0) = r1.S;

        Matrix<double> n1 = g1 * d1;

        ASSERT_NEAR(0.0, n1(1, 0), 0.00000001);
        ASSERT_NEAR(std::sqrt(2), std::abs(n1(0, 0)), 0.000001);
        ASSERT_NEAR(n1(0,0), r1.R, 0.000001);
    }
}

TEST(Decomposition, GivensRotationBatch)
{
    int n = 10000;

    for( int k = 0; k < n; k++ )
    {
        Matrix<double> input = Matrix<double>::random(2,1,-100.0, 100.0);
        Decomposition::GivensRotation res = Decomposition::givensRotation(input(0, 0), input(1, 0));

        Matrix<double> g(2, 2);
        g(0, 0) = res.C;
        g(1, 1) = res.C;
        g(0, 1) = -res.S;
        g(1, 0) = res.S;

        Matrix<double> n = g * input;

        ASSERT_NEAR(0.0, n(1, 0), 0.00000001);
        ASSERT_NEAR(n(0,0), res.R, 0.000001);
    }
}

// https://en.wikipedia.org/wiki/Givens_rotation#Triangularization
TEST(Decomposition, GivensRotationMatrix)
{
    double data[] = {6,5,0 ,5,1,4, 0,4,3};
    Matrix<double> a = Matrix<double>(3,3,data);

    Matrix<double> g0 = Decomposition::givensRotatioColumnDirection(a, 0, 0, 1);
    Matrix<double> a1 = g0 * a;

    ASSERT_NEAR( 0.0, a1(1,0), 0.0001 );

    Matrix<double> g1 = Decomposition::givensRotatioColumnDirection(a1, 1, 1, 2);
    Matrix<double> a2 = g1 * a1;

    ASSERT_NEAR( 0.0, a2(2,1), 0.0001 );

    Matrix<double> g = g1 * g0;
    Matrix<double> r = g * a;

    ASSERT_TRUE( r.compare(a2,true,0.001) );
}

TEST(Decomposition, GivensRotationMatrixBatchCol)
{
    for(int k = 0; k < 1000; k++ )
    {
        // Make uper triangle matrix

        Matrix<double> a = Matrix<double>::random(4, 4, -2.0, 2.0);

        a = Decomposition::givensRotatioColumnDirection(a, 0, 2, 3) * a;
        a = Decomposition::givensRotatioColumnDirection(a, 0, 1, 2) * a;
        a = Decomposition::givensRotatioColumnDirection(a, 0, 0, 1) * a;

        a = Decomposition::givensRotatioColumnDirection(a, 1, 1, 2) * a;
        a = Decomposition::givensRotatioColumnDirection(a, 1, 1, 3) * a;

        a = Decomposition::givensRotatioColumnDirection(a, 2, 2, 3) * a;

        ASSERT_NEAR( a(1,0), 0.0, 0.000001 );
        ASSERT_NEAR( a(2,0), 0.0, 0.000001 );
        ASSERT_NEAR( a(3,0), 0.0, 0.000001 );

        ASSERT_NEAR( a(2,1), 0.0, 0.000001 );
        ASSERT_NEAR( a(3,1), 0.0, 0.000001 );

        ASSERT_NEAR( a(3,2), 0.0, 0.000001 );
    }
}

TEST(Decomposition, GivensRotationMatrixZeroUpperRowElement)
{
    double data[] = {1,2,3,4,
                     5,6,7,8,
                     9,10,11,12,
                     13,14,15,16};
    Matrix<double> a = Matrix<double>(4,4,data);
    Matrix<double> orig = Matrix<double>(4,4,data);

    a = Decomposition::givensRotatioColumnDirection(a, 1, 1, 0) * a;

    ASSERT_NEAR( a(0,1), 0.0, 0.000001 );
    a = Decomposition::givensRotatioColumnDirection(a, 2, 1, 0) * a;
    ASSERT_NEAR( a(0,2), 0.0, 0.000001 );
    a = Decomposition::givensRotatioColumnDirection(a, 3, 1, 0) * a;
    ASSERT_NEAR( a(0,3), 0.0, 0.000001 );
}

TEST(Decomposition, GivensRotationMatrixCol)
{
    double data[] = {1,2,3, 4,5,6, 7,8,3};
    Matrix<double> mat(3,3,data);

    Matrix<double> g1 = Decomposition::givensRotationRowDirection(mat, 0, 1, 2);
    Matrix<double> a = mat * g1;
    Matrix<double> g2 = Decomposition::givensRotationRowDirection(a, 0, 0, 1);
    a = a * g2;
    Matrix<double> g3 = Decomposition::givensRotationRowDirection(a, 1, 1, 2);
    a = a * g3;

    ASSERT_NEAR( a(0,1), 0.0, 0.000001 );
    ASSERT_NEAR( a(0,2), 0.0, 0.000001 );
    ASSERT_NEAR( a(1,2), 0.0, 0.000001 );

    ASSERT_TRUE( mat.compare( a * g3.inverted() * g2.inverted() * g1.inverted(), true, 0.00001) );
}

TEST(Decomposition, GivensRotationMatrixRowBatch)
{
    for(int k = 0; k < 1000; k++ )
    {
        Matrix<double> a = Matrix<double>::random(3, 3, -2.0, 2.0);

        // Make lower triangle matrix
        a = a * Decomposition::givensRotationRowDirection(a, 0, 1, 2);
        a = a * Decomposition::givensRotationRowDirection(a, 0, 0, 1);
        a = a * Decomposition::givensRotationRowDirection(a, 1, 1, 2);

        ASSERT_NEAR(a(0, 1), 0.0, 0.000001);
        ASSERT_NEAR(a(0, 2), 0.0, 0.000001);
        ASSERT_NEAR(a(1, 2), 0.0, 0.000001);
    }
}

TEST(Decomposition, SVDGolubKahan)
{ /*
    >> a = [1,2,1; -1,1,2; 4, 2, 2]
    >> [u,s,v] = svd(a)
    u =

        -0.41178  -0.28554  -0.86539
        -0.13290  -0.92067   0.36702
        -0.90154   0.26614   0.34116

    s =

        Diagonal Matrix

    5.37236         0         0
    0   2.52038         0
    0         0   0.88624

    v =

        -0.72315   0.67438   0.14920
        -0.51365  -0.38069  -0.76892
        -0.46174  -0.63268   0.62169
    */

    double matData[] = {1,2,1, -1,1,2, 4, 2, 2};
    auto mat = Matrix<double>(3, 3, matData);

    Matrix<double> s_soll = Matrix<double>(3,3);
    s_soll.fill(0.0);
    s_soll(0,0) = 5.37236;
    s_soll(1,1) = 2.52038;
    s_soll(2,2) = 0.88624;

    Decomposition::SVDResult res = Decomposition::svdGolubKahan(mat);

    // U and V cannot be compared directly, because it can change in signs.
    ASSERT_TRUE( res.U.isOrthogonal(0.00001) );
    ASSERT_TRUE( res.V.isOrthogonal(0.00001) );
    ASSERT_TRUE( s_soll.compare(res.S, true, 0.001) );
    ASSERT_TRUE( mat.compare(res.U * res.S * res.V.transpose(), true, 0.00001 ) );
}

TEST(Decomposition, SVDGolubKahanWithNullSingularValue)
{
    double matData[] = {0.0, 1.0, 1.0,
                        std::sqrt(2.0), 2.0, 0,
                        0.0, 1.0, 1.0};

    auto a = Matrix<double>(3, 3, matData);

    //double matdata[] = {1,2,1,3,4,  4,3,7,1,2, 8,9,1,3,3, 6,6,6,3,2, 4,4,6,3,1};
    //size_t rows = 5;
    //Matrix<double> a(rows,rows, matdata);

    Decomposition::SVDResult res = Decomposition::svdGolubKahan(a);

    ASSERT_TRUE( res.U.isOrthogonal(0.05) );
    ASSERT_TRUE( res.V.isOrthogonal(0.05) );
    ASSERT_TRUE( a.compare(res.U * res.S * res.V.transpose(), true, 0.1 ) );
}

TEST(Decomposition, SVDGolubKahanMultipleNullSingularValues)
{
    Matrix<double> m = Matrix<double>(3, 3);
    m.fill(1);

    Decomposition::SVDResult res = Decomposition::svdGolubKahan(m);

    std::cout << "Computed Results" << std::endl << std::endl;


    std::cout << "Singualar values " << std::endl;
    std::cout << res.S << std::endl << std::endl;

    std::cout << "U (left) " << std::endl;
    std::cout << res.U << std::endl << std::endl;

    std::cout << "V (right) " << std::endl;
    std::cout << res.V << std::endl << std::endl;

    ASSERT_TRUE( res.U.isOrthogonal(0.00001) );
    ASSERT_TRUE( res.V.isOrthogonal(0.00001) );
    ASSERT_TRUE( m.compare(res.U * res.S * res.V.transpose(), true, 0.00001 ) );


    /* From octave
U =
    -0.57735      0.8165  -8.7561e-17
    -0.57735    -0.40825    -0.70711
    -0.57735    -0.40825     0.70711

S =

           3           0           0
           0   1.338e-16           0
           0           0  1.3194e-49

V =
    -0.57735      0.8165          -0
    -0.57735    -0.40825     0.70711
    -0.57735    -0.40825    -0.70711
     */
}

TEST(Decomposition, SVDGolubKahanBatch)
{
    size_t rows = 5;
    size_t n_test = 10000;

    for( size_t k = 0; k < n_test; k++ )
    {
        auto a = Matrix<double>::random(rows, rows, -100.0, 100.0);

        Decomposition::SVDResult res = Decomposition::svdGolubKahan(a);


        ASSERT_TRUE( res.U.isOrthogonal(0.0000001) );
        ASSERT_TRUE( res.V.isOrthogonal(0.0000001) );

        bool isAccurate = a.compare(res.U * res.S * res.V.transpose(), true, 0.00001 );

        if(!isAccurate)
        {
            a.toMatlab();
            std::cout << std::endl;
        }


        ASSERT_TRUE( isAccurate );

        if(k % 1000 == 0 )
            std::cout << "SVD Batch testing: " << (double)(k) / n_test * 100.0 << " %" << std::endl;
    }
}

TEST(Decomposition, SVDGolubKahanMatrixCheckIdentity)
{
    size_t cols = 5;
    Matrix<double> m = Matrix<double>::identity(cols);

    size_t p, q;
    Decomposition::svdCheckMatrixGolubKahan(m,p,q);

    ASSERT_EQ(q,cols);
    ASSERT_EQ(p,0);


    // one entry
    cols = 1;
    m = Matrix<double>::identity(cols);

    Decomposition::svdCheckMatrixGolubKahan(m,p,q);

    ASSERT_EQ(q,cols);
    ASSERT_EQ(p,0);
}

TEST(Decomposition, SVDGolubKahanMatrixCheckLowerIdentUpperBiag)
{
    size_t cols = 6;
    Matrix<double> m = Matrix<double>::identity(cols);
    m(0,1) = 2;
    m(1,2) = 3;

    size_t p, q;
    Decomposition::svdCheckMatrixGolubKahan(m,p,q);

    ASSERT_EQ(q,cols-3);
    ASSERT_EQ(p,0);
}

TEST(Decomposition, SVDGolubKahanMatrixCheckLowerIdentMiddleBiagUppestDiag)
{
    size_t cols = 6;
    Matrix<double> m = Matrix<double>::identity(cols);
    m(1,2) = 3;

    size_t p, q;
    Decomposition::svdCheckMatrixGolubKahan(m,p,q);

    ASSERT_EQ(q,cols-3);
    ASSERT_EQ(p,1);
}

TEST(Decomposition, SVDGolubKahanMatrixCheckAllBiag)
{
    size_t cols = 4;
    Matrix<double> m = Matrix<double>::identity(cols);
    m(0,1) = 3;
    m(1,2) = 3;
    m(2,3) = 3;

    size_t p, q;
    Decomposition::svdCheckMatrixGolubKahan(m,p,q);

    ASSERT_EQ(q,0);
    ASSERT_EQ(p,0);
}

TEST(Decomposition, SVDGolubKahanMatrixCheckUpDiag)
{
    size_t cols = 4;
    Matrix<double> m = Matrix<double>::identity(cols);
    m(1,2) = 3;
    m(2,3) = 3;

    size_t p, q;
    Decomposition::svdCheckMatrixGolubKahan(m,p,q);

    ASSERT_EQ(q,0);
    ASSERT_EQ(p,1);
}

TEST(Decomposition, SVDGolubKahanPaddingGeneral)
{
    for( size_t p = 0; p < 5; p++ )
    {
        for( size_t q = 0; q < 5; q++ )
        {
            for( size_t n = 1; n < 5; n++ )
            {
                Matrix<double> mat = Matrix<double>::identity(n);
                Matrix<double> pad = Decomposition::svdPaddingRotation(mat,p,q);
                Matrix<double> should = Matrix<double>::identity(n+q+p);

                ASSERT_TRUE( pad.compare(should) );
            }
        }
    }
}

TEST(Decomposition, SVDGolubKahanPaddingMat)
{
    double shouldData[] = {1.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 1.5, 2.0, 0.0, 0.0,
                           0.0, 3.0, 4.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 1.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 1.0};

    double matData[] = { 1.5, 2.0,
                         3.0, 4.0};

    Matrix<double> should = Matrix<double>(5,5,shouldData);
    Matrix<double> mat = Matrix<double>(2,2,matData);

    Matrix<double> padded = Decomposition::svdPaddingRotation(mat, 1 , 2);

    ASSERT_TRUE( padded.compare(should) );
}

TEST(Decomposition, SVDGolubKahanZeroRow)
{
    double matData[] =    {1.0, 3.2, 0.0, 0.0, 0.0,
                           0.0, 2.0, 1.8, 0.0, 0.0,
                           0.0, 0.0, 0.0, 1.2, 0.0,
                           0.0, 0.0, 0.0, 1.0, 1.4,
                           0.0, 0.0, 0.0, 0.0, 1.0};

    Matrix<double> mat = Matrix<double>(5,5,matData);
    Matrix<double> orig = mat;
    Matrix<double> givens = Decomposition::svdZeroRow(mat,2);

    // check that the 3rd row is zero
    for( size_t i = 0; i < 5; i++ )
        ASSERT_FLOAT_EQ( mat(2,i), 0.0 );

    ASSERT_TRUE( mat.compare(givens * orig, true, 0.000001) );

    std::cout << "SVDGolubKahanZeroRow" << std::endl << mat << std::endl;
}

TEST(Decomposition, SVDGolubKahanZeroColumnStepByStep)
{
    Matrix<double> mat = Matrix<double>(4,4, {1.0, 1.0, 0.0, 0.0,
                                              0.0, 1.0, 1.0, 0.0,
                                              0.0, 0.0, 1.0, 1.0,
                                              0.0, 0.0, 0.0, 0.0}); //<<<<<
    Matrix<double> original = mat;

    auto g0 = Decomposition::givensRotationRowDirection(mat, 2,  2, 3); // sets element 2,3 to zero
    mat = mat * g0;

    auto g1 = Decomposition::givensRotationRowDirection(mat, 1,  1, 3); // sets element 1,3 to zero
    mat = mat * g1;

    auto g2 = Decomposition::givensRotationRowDirection(mat, 0,  0, 3); // sets element 0,3 to zero
    mat = mat * g2;

    // assure last row and last column are zero
    for(size_t k = 0; k < 4; k++)
    {
        ASSERT_NEAR(mat(3,k), 0.0, 0.00001);
        ASSERT_NEAR(mat(k,3), 0.0, 0.00001);
    }
}

TEST(Decomposition, SVDGolubKahanZeroColumnLast)
{
    Matrix<double> mat = Matrix<double>(4,4, {1.0, 1.0, 0.0, 0.0,
                                              0.0, 1.0, 1.0, 0.0,
                                              0.0, 0.0, 1.0, 1.0,
                                              0.0, 0.0, 0.0, 0.0}); //<<<<<

    Matrix<double> orig = mat;
    Matrix<double> givens = Decomposition::svdZeroColumn(mat,3);

    // assure last row and last column are zero
    for(size_t k = 0; k < 4; k++)
    {
        ASSERT_NEAR(mat(3,k), 0.0, 0.00001);
        ASSERT_NEAR(mat(k,3), 0.0, 0.00001);
    }

    ASSERT_TRUE( mat.compare( orig * givens, true, 0.000001) );

    std::cout << "SVDGolubKahanZeroColumnLast" << std::endl << mat << std::endl;
}

TEST(Decomposition, SVDSort)
{
    Matrix<double> s_unsorted =  Matrix<double>(4,4, {1.0, 0.0, 0.0, 0.0,
                                                      0.0, 2.0, 0.0, 0.0,
                                                      0.0, 0.0, 3.0, 0.0,
                                                      0.0, 0.0, 0.0, 0.5});

    Matrix<double> u_unsorted = Matrix<double>::identity(4);
    Matrix<double> v_unsorted = Matrix<double>::identity(4);

    Decomposition::SVDResult unsorted(u_unsorted, s_unsorted, v_unsorted);

    Decomposition::SVDResult sorted = Decomposition::sortSingularValues(unsorted);


    std::cout << "Singualar values " << std::endl;
    std::cout << sorted.S << std::endl << std::endl;

    std::cout << "U (left) " << std::endl;
    std::cout << sorted.U << std::endl << std::endl;

    std::cout << "V (right) " << std::endl;
    std::cout << sorted.V << std::endl << std::endl;

    Matrix<double> a_unsorted = unsorted.U * unsorted.S * unsorted.V.transpose();
    Matrix<double> a_sorted = sorted.U * sorted.S * sorted.V.transpose();
    ASSERT_TRUE(a_unsorted.compare(a_sorted, true, 0.00001));

    ASSERT_TRUE( sorted.U.isOrthogonal(0.00001) );
    ASSERT_TRUE( sorted.V.isOrthogonal(0.00001) );

    // check ascending s
    double s = 4.0; // biggest value is 3
    for(size_t k = 0; k < 4; k++)
    {
        ASSERT_GE(s, sorted.S(k,k));
        s = sorted.S(k,k);
    }
}

TEST(Decomposition, SVDGolubKahanInaccurateMat)
{
    // Golub Kahan method shows accuracy problems with this very matrix

    /* Octave

 a =

  -0.102025  -0.723523   0.891811   0.481759  -0.633212
  -0.348421   0.711962   0.531478   0.935420  -0.097372
  -0.018528  -0.464990   0.812909  -0.601147  -0.149676
  -0.195534   0.625486   0.594395  -0.914214   0.648131
  -0.042660   0.641828  -0.261626  -0.363415   0.391454

>>
>> [U, S, V] = svd (a)
U =

  -0.683046   0.311696   0.171736   0.049374   0.635895
  -0.195915  -0.334253   0.870418  -0.137152  -0.271025
  -0.138227   0.706147   0.015195  -0.521594  -0.458212
   0.525910   0.534706   0.443223   0.468542   0.146728
   0.446517  -0.081030   0.127288  -0.697968   0.539161

S =

Diagonal Matrix

   1.895558          0          0          0          0
          0   1.544089          0          0          0
          0          0   1.354855          0          0
          0          0          0   0.107890          0
          0          0          0          0   0.014395

V =

   0.0098270  -0.0191181  -0.3049551  -0.0873793  -0.9481061
   0.5457635  -0.3299046   0.6253882  -0.4239857  -0.1497695
  -0.3322823   0.6563000   0.6334728   0.0763477  -0.2274687
  -0.5656895  -0.6776752   0.3220629   0.3183853  -0.1251318
   0.5211809   0.0287058   0.1043058   0.8398827  -0.1061316

     */

    Matrix<double> m = Matrix<double>(5, 5, {-0.10202521252727214662,-0.72352319888331184661,0.89181119050052326536,0.48175927978607036017,-0.63321242927972276604,-0.34842097643748359825,0.71196246658029571641,0.53147834081151068553,0.93542036356107161055,-0.097372463543626719407,-0.018527628599318224367,-0.46498953560260436468,0.81290852214630815453,-0.60114650914488110267,-0.14967587343376720366,-0.19553402337527159283,0.62548581871254160802,0.59439451823066846714,-0.91421383569353364962,0.64813108152028542364,-0.042660400203598625168,0.64182817888695598008,-0.26162605672833538772,-0.36341526642480959097,0.39145355245384694243});


    auto eigen = Decomposition::svdEigen(m);
    auto golub = Decomposition::svdGolubKahan(m);

    std::cout << "Eigen SVD: " << std::endl << eigen.U << std::endl << eigen.S << std::endl<< eigen.V << std::endl;
    std::cout << "Golub SVD: " << std::endl << golub.U << std::endl << golub.S << std::endl << golub.V << std::endl;
}
