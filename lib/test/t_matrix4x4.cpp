#include <gtest/gtest.h>
#include "matrix4x4.hpp"

TEST(Matrix4x4, Constructor)
{
    Matrix4x4 mat = Matrix4x4();

    ASSERT_EQ(mat.cols(), 4);
    ASSERT_EQ(mat.rows(), 4);
    ASSERT_TRUE(Matrix<double>::identity(4).compare(mat));
}

TEST(Matrix4x4, ConstructorData)
{
    Matrix4x4 mat = Matrix4x4(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    for( size_t m = 0; m < 4; m++ )
    {
        for( size_t n = 0; n < 4; n++ )
        {
            double soll = static_cast<double>( m*4 + n );
            ASSERT_FLOAT_EQ(soll, mat(m,n) );
        }
    }
}

TEST(Matrix4x4, CopyConstructorMatrix)
{
    Matrix<double> matD = Matrix<double>(4,4);
    matD.fill(5.0);

    Matrix4x4 mat = matD;

    ASSERT_TRUE(mat.compare(matD));
}

TEST(Matrix4x4, CopyConstructor)
{
    Matrix4x4 a = Matrix4x4();
    a.fill(6.0);

    Matrix4x4 b = a;

    ASSERT_TRUE(b.compare(a));
}


TEST(Matrix4x4, AssignMatrixDouble)
{
    auto mat = Matrix4x4();

    auto matD = Matrix<double>(4,4);
    matD.fill(3.0);

    mat = matD;

    ASSERT_TRUE(mat.compare(matD));
}


TEST(Matrix4x4, Mulitplication4x4)
{
    Matrix4x4 a = Matrix4x4();
    Matrix4x4 b = Matrix4x4();
    for (size_t i = 0; i < 16; i++)
    {
        a.data()[i] = 1.0 * i;
        b.data()[i] = 2.0 * i;
    }

    Matrix4x4 soll = a.matMulR(b); // generic
    Matrix4x4 res = a * b;

    ASSERT_TRUE(soll.compare(res));
}

TEST(Matrix4x4, Mulitplication4x4Performance)
{
    auto r = Matrix4x4();
    r.rotZ(M_PI);
    auto soll = r;

    auto d = Matrix4x4();
    d.setToIdentity();

    // test Matrix4x4 multiplication
    for( int i = 0; i < 100000; i++ )
    {
        r = d*r;
    }
}

TEST(Matrix4x4, MulitplicationGenericPerformance)
{
    auto r = Matrix4x4();
    r.rotZ(M_PI);
    auto soll = r;

    auto d = Matrix4x4();
    d.setToIdentity();

    // test general multiplication
    for( int i = 0; i < 100000; i++ )
    {
        r = d.matMulR(r);
    }
}

TEST(Matrix4x4, RotZ)
{
    Matrix4x4 mat = Matrix4x4();
    mat.rotZ(M_PI/2.0); // rotate 90

    double v1Data[] = {2,1,7,1};
    auto v1 = Matrix<double>(4,1,v1Data);

    auto r1 = mat * v1;

    double r1SollData[] = {-1,2,7,1};
    auto r1Soll = Matrix<double>(4,1,r1SollData);

    ASSERT_TRUE(r1Soll.compare(r1));

    mat.rotZ(3*M_PI/2.0); // With above rotation, this makes a whole 360 degree
    r1 = mat * v1;
    ASSERT_TRUE( v1.compare(r1) );
}

