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

TEST(Matrix4x4, RotY)
{
    Matrix4x4 mat = Matrix4x4();
    mat.rotY(M_PI/2.0); // rotate 90

    double v1Data[] = {2,1,7,1};
    auto v1 = Matrix<double>(4,1,v1Data);

    auto r1 = mat * v1;

    double r1SollData[] = {7,1,-2,1};
    auto r1Soll = Matrix<double>(4,1,r1SollData);

    ASSERT_TRUE(r1Soll.compare(r1));

    mat.rotY(3*M_PI/2.0); // With above rotation, this makes a whole 360 degree
    r1 = mat * v1;
    ASSERT_TRUE( v1.compare(r1, true, 0.00001) );
}

TEST(Matrix4x4, RotX)
{
    Matrix4x4 mat = Matrix4x4();
    mat.rotX(M_PI/2.0); // rotate 90

    double v1Data[] = {2,1,7,1};
    auto v1 = Matrix<double>(4,1,v1Data);

    auto r1 = mat * v1;

    double r1SollData[] = {2,-7,1,1};
    auto r1Soll = Matrix<double>(4,1,r1SollData);

    ASSERT_TRUE(r1Soll.compare(r1));

    mat.rotX(3*M_PI/2.0); // With above rotation, this makes a whole 360 degree
    r1 = mat * v1;

    ASSERT_TRUE( v1.compare(r1, true, 0.00001) );
}

TEST(Matrix4x4, RigidInverse)
{
    Matrix4x4 t;
    t.rotZ(M_PI/2.0); // set z-rot
    t.rotY(M_PI/3.0); // set y-rot
    t(0,3) = 2.0; // set x-translation

    bool invertable;
    Matrix4x4 t_slowComputedInv = t.inverted(&invertable); // true
    ASSERT_TRUE(invertable);

    Matrix4x4 t_fastInv = t.inverted_rg();

    ASSERT_TRUE(t_slowComputedInv.compare(t_fastInv,true,0.0000001));

    Matrix4x4 ident;
    ASSERT_TRUE(ident.compare(t_fastInv * t, true, 0.0000001));
}


TEST(Matrix4x4, SetAndGetRotation)
{
    double rotData[] = {1,2,3,  4,5,6,  7,8,9};
    auto rot = Matrix<double>(3,3, rotData);

    Matrix4x4 t = Matrix4x4();
    t.setRotation(rot);

    auto rotGet = t.getRotation();

    ASSERT_TRUE(rot.compare(rotGet));
}

TEST(Matrix4x4, SetAndGetTranslation)
{
    double transData[] = {1,2,3};
    auto trans = Matrix<double>(3,1, transData);

    Matrix4x4 t = Matrix4x4();
    t.setTranslation(trans);
    auto transGet = t.getTranslation();

    ASSERT_TRUE(trans.compare(transGet));

    Matrix4x4 t2 = Matrix4x4();
    t2.setTranslation(1,2,3);
    auto getTrans2 = t2.getTranslation();

    ASSERT_TRUE(trans.compare(getTrans2));
}
