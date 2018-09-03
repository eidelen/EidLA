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

    Matrix4x4 t_slowComputedInv = t.inverted();

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

    double transData4[] = {1,2,3,1.0};
    auto trans4 = Matrix<double>(4,1, transData4);

    Matrix4x4 t = Matrix4x4();
    t.setTranslation(trans);
    auto transGet = t.getTranslation();

    ASSERT_TRUE(trans.compare(transGet));

    t.setTranslation(trans4);
    transGet = t.getTranslation();

    ASSERT_TRUE(trans.compare(transGet));

    Matrix4x4 t2 = Matrix4x4();
    t2.setTranslation(1,2,3);
    auto getTrans2 = t2.getTranslation();

    ASSERT_TRUE(trans.compare(getTrans2));

    // exceptions
    auto throwM1 = Matrix<double>(2,1);
    auto throwM2 = Matrix<double>(5,1);
    auto throwM3 = Matrix<double>(3,2);
    EXPECT_ANY_THROW( t2.setTranslation(throwM1));
    EXPECT_ANY_THROW( t2.setTranslation(throwM2));
    EXPECT_ANY_THROW( t2.setTranslation(throwM3));
}

TEST(Matrix4x4, FindAffineTransformation)
{
    // Generate test point cloud
    size_t nInput = 100;
    Matrix<double> p_A = Matrix<double>::random(4,nInput,-1,1);
    for(size_t k = 0; k < nInput; k++)
        p_A(3,k) = 1.0;

    // Make example transformation
    Matrix4x4 t = Matrix4x4();
    t.setTranslation( 1, 2, 3 );
    t.rotX(0.2);
    t.rotY(0.5);
    t.rotZ(-0.1);

    // Transform p_A -> corresponding point cloud
    Matrix<double> p_B = t * p_A;

    // this function needs to recover transformation t
    double err = 0.0;
    auto t_res = Matrix4x4::findRigidTransformation(p_A.subMatrix(0,0,3,nInput), p_B.subMatrix(0,0,3,nInput), err);

    ASSERT_TRUE(t_res.compare(t, true, 0.01));
    ASSERT_LT(err, 0.01);
}

TEST(Matrix4x4, FindAffineTransformation_3Points)
{
    // Generate test point cloud
    double p_A_data[] = {1 , 0, 0,
                         0 , 1, 0,
                         0,  0, 1,
                         1,  1, 1};
    Matrix<double> p_A = Matrix<double>(4,3, p_A_data);

    // Make example transformation
    Matrix4x4 t = Matrix4x4();
    t.setToIdentity();
    t.setTranslation(0.0, 0.3, 0.5);
    t.rotY(0.6);

    // Transform p_A -> corresponding point cloud
    Matrix<double> p_B = t * p_A;

    // this function needs to recover transformation t
    double err = 0.0;
    auto t_res = Matrix4x4::findRigidTransformation(p_A.subMatrix(0,0,3,3), p_B.subMatrix(0,0,3,3), err);
    ASSERT_TRUE(t_res.compare(t, true, 0.01));
    ASSERT_LT(err, 0.001);
}


TEST(Matrix4x4, FindAffineTransformation_Error)
{
    // Generate test point cloud
    double p_A_data[] = {1 , 2, 3, 2,
                         1 , 0, 1, 1,
                         0,  1, 0, 3,
                         1,  1, 1, 1};
    Matrix<double> p_A = Matrix<double>(4,4, p_A_data);

    // Make example transformation
    Matrix4x4 t = Matrix4x4();
    t.setTranslation( 1, 0, 0 );
    t.rotY(0.5);

    // Transform p_A -> corresponding point cloud
    Matrix<double> p_B = t * p_A;
    p_B(0,0) = p_B(0,0) + 0.1; // modify a bit

    // this function needs to recover transformation t
    double err = 0.0;
    auto t_res = Matrix4x4::findRigidTransformation(p_A.subMatrix(0,0,3,4), p_B.subMatrix(0,0,3,4), err);
    ASSERT_GT(err, 0.01);
}

TEST(Matrix4x4, FindAffineTransformation_Batch)
{
    size_t runs = 1000;

    Matrix<size_t> nbrInputPoints = Matrix<size_t>::random(runs,1,3,20);

    for(size_t run = 0; run < runs; run++)
    {
        size_t n = nbrInputPoints(run,0);

        Matrix<double> inputA = Matrix<double>::random( 4, n, -10, +10 );
        Matrix<double> h = Matrix<double>(1,n); h.fill(1);
        inputA.setRow(3, h );

        Matrix4x4 t = Matrix4x4();
        t.setToIdentity();
        Matrix<double> t_rand = Matrix<double>::random( 6, 1, -3, 3 );
        t.rotX(t_rand(0,0)); t.rotY(t_rand(1,0)); t.rotZ(t_rand(2,0));
        t.setTranslation(t_rand(3,0),t_rand(4,0),t_rand(5,0));

        Matrix<double> inputB = t * inputA;

        double error;
        auto t_res = Matrix4x4::findRigidTransformation(inputA.subMatrix(0,0,3,n), inputB.subMatrix(0,0,3,n), error);

        if( error > 0.05 )
            std::cout << "Error exceeded 0.05: " << error << std::endl;

        if( ! t.compare(t_res, true, 0.01) )
            std::cout << "T differs by more than 0.01 " << std::endl;


        // todo: This values need to be lower! -> Decrease after implementing a new SVD
        ASSERT_TRUE( t.compare(t_res, true, 0.02) );
        ASSERT_LT( error, 0.05 );

        if( run % 100 == 0 )
            std::cout << "Run: " << run << " of " << runs << std::endl;
    }
}

TEST(Matrix4x4, ComputeTransError)
{
    double p_A_data[] = {1 , 2,
                         1 , 2,
                         0,  0};
    Matrix<double> p_A = Matrix<double>(3,2, p_A_data);

    double p_B_data[] = {1 , 2,
                         1 , 2,
                         0,  0};
    Matrix<double> p_B = Matrix<double>(3,2, p_B_data);

    Matrix4x4 t;
    t.setToIdentity();

    ASSERT_LT(Matrix4x4::computeTransformationError(p_A, p_B, t), 0.00001);


    p_B(0,0) = 1.5; // error = sqrt(0.5^2)/2 = 0.25
    ASSERT_NEAR( Matrix4x4::computeTransformationError(p_A, p_B, t), 0.25, 0.001);
}

