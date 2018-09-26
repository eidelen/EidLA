#include <gtest/gtest.h>
#include "interpolation.hpp"

TEST(Interpolation, LinearInterpolation)
{
    double data[] = {0,0,
                     1,1,
                     2,3,
                     3,0};

    Matrix<double> mat(4,2,data);

    LinearInterpolation* li = new LinearInterpolation(mat);

    ASSERT_FLOAT_EQ( li->interpolate(-100), -100 );
    ASSERT_FLOAT_EQ( li->interpolate(0.5), 0.5 );
    ASSERT_FLOAT_EQ( li->interpolate(-0.5), -0.5 );
    ASSERT_FLOAT_EQ( li->interpolate(1.0), 1.0 );
    ASSERT_FLOAT_EQ( li->interpolate(1.5), 2.0 );
    ASSERT_FLOAT_EQ( li->interpolate(2.5), 1.5 );
    ASSERT_FLOAT_EQ( li->interpolate(3.0), 0.0 );
    ASSERT_FLOAT_EQ( li->interpolate(3.5), -1.5 );

    delete li;
}

TEST(Interpolation, LinearInterpolationUnorderedInput)
{
    double data[] = {1,1,
                     0,0,
                     3,0,
                     2,3};

    Matrix<double> mat(4,2,data);

    LinearInterpolation* li = new LinearInterpolation(mat);

    ASSERT_FLOAT_EQ( li->interpolate(-100), -100 );
    ASSERT_FLOAT_EQ( li->interpolate(0.5), 0.5 );
    ASSERT_FLOAT_EQ( li->interpolate(-0.5), -0.5 );
    ASSERT_FLOAT_EQ( li->interpolate(1.0), 1.0 );
    ASSERT_FLOAT_EQ( li->interpolate(1.5), 2.0 );
    ASSERT_FLOAT_EQ( li->interpolate(2.5), 1.5 );
    ASSERT_FLOAT_EQ( li->interpolate(3.0), 0.0 );
    ASSERT_FLOAT_EQ( li->interpolate(3.5), -1.5 );

    delete li;
}

TEST(Interpolation, CubicInterpolation)
{
    double data[] = {0,0,
                     1,1,
                     2,0,
                     3,-1,
                     4,0};

    Matrix<double> mat(5,2,data);

    CubicSplineInterpolation* ci = new CubicSplineInterpolation(mat);

    // check that interpolation works for stützpunkte
    for( size_t m = 0; m < mat.rows(); m++ )
    {
        double x = mat(m,0);
        double y = mat(m,1);
        ASSERT_FLOAT_EQ(y, ci->interpolate(x));
    }

    // check that function is smooth - changes in derivation

    size_t n = 1000000;
    Matrix<double> xVals(n,1);
    double stepSize = 4.0 / n;
    for( size_t k = 0; k < n; k++ )
        xVals(k,0) = 0.0 + stepSize * k;

    Matrix<double> functionValues = ci->interpolate(xVals);
    Matrix<double> derivation(n-1,1);
    for(size_t k = 0; k < n - 1; k++ )
        derivation(k,0) = ( functionValues(k+1,1) - functionValues(k,1) ) / ( functionValues(k+1,0) - functionValues(k,0) );

    // no jump in derivation
    for(size_t k = 0; k < n - 2; k++ )
        ASSERT_LT( std::abs(derivation(k,0) - derivation(k+1,0)), 0.0001 );

    delete ci;
}

TEST(Interpolation, IsInRange)
{
    class A: public Interpolation1D
    {
    public:
        A(){};
        void testRange()
        {
            range1D r0; r0.startR = -1.0; r0.endR = 0.0;

            ASSERT_FALSE(isInRange(r0,-1.1));
            ASSERT_TRUE(isInRange(r0,-1.0));
            ASSERT_TRUE(isInRange(r0,-0.5));
            ASSERT_TRUE(isInRange(r0,0.0));
            ASSERT_FALSE(isInRange(r0,0.1));
        }
        double interpolate(double x) const override {return 0.0;}
    };

    A* a = new A();
    a->testRange();
    delete a;
}
