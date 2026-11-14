#include <gtest/gtest.h>
#include "matrix.hpp"

#include <cmath>
#include <limits>


TEST(Accuracy, GivensRotation)
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double tolerance = 16.0 * eps;

    // Scaling a vector must not change whether its second component is eliminated.
    // The tiny case has b = eps / 2: the absolute epsilon cutoff incorrectly
    // returns the identity rotation, leaving a relative residual of 4 / 5.
    for (const double scale : {1.0, eps / 8.0})
    {
        SCOPED_TRACE(::testing::Message() << "scale = " << scale);
        const double a = 3.0 * scale;
        const double b = 4.0 * scale;
        const double magnitude = std::hypot(a, b);
        const auto rotation = Decomposition::givensRotation(a, b);

        // Apply the rotation using the library's sign convention.
        const double rotatedA = rotation.C * a - rotation.S * b;
        const double rotatedB = rotation.S * a + rotation.C * b;

        EXPECT_NEAR(rotatedB / magnitude, 0.0, tolerance);
        EXPECT_NEAR(rotatedA / magnitude, 1.0, tolerance);
        EXPECT_NEAR(rotation.R / magnitude, 1.0, tolerance);
        EXPECT_NEAR(rotation.C * rotation.C + rotation.S * rotation.S,
                    1.0, tolerance);
    }
}

TEST(Accuracy, GivensRotationExtremeScales)
{
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon();
    for (const double scale : {std::numeric_limits<double>::denorm_min(),
                               std::numeric_limits<double>::min(),
                               std::numeric_limits<double>::max() / 8.0})
    {
        SCOPED_TRACE(::testing::Message() << "scale = " << scale);
        const auto rotation = Decomposition::givensRotation(3.0 * scale, 4.0 * scale);

        // Check coefficients without subnormal rounding in the rotated vector.
        EXPECT_NEAR(rotation.C, 0.6, tolerance);
        EXPECT_NEAR(rotation.S, -0.8, tolerance);
        EXPECT_NEAR(rotation.R / scale, 5.0, tolerance);
    }
}

TEST(Accuracy, GivensRotationZerosAndSigns)
{
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon();
    for (const double a : {-3.0, -0.0, 0.0, 3.0})
    {
        for (const double b : {-4.0, 0.0, 4.0})
        {
            SCOPED_TRACE(::testing::Message() << "a = " << a << ", b = " << b);
            const auto rotation = Decomposition::givensRotation(a, b);
            EXPECT_NEAR(rotation.S * a + rotation.C * b, 0.0, tolerance);
            EXPECT_NEAR(rotation.C * a - rotation.S * b, rotation.R, tolerance);
            EXPECT_NEAR(rotation.R, std::copysign(std::hypot(a, b), a), tolerance);
            EXPECT_NEAR(rotation.C * rotation.C + rotation.S * rotation.S,
                        1.0, tolerance);
            if (b == 0.0)
            {
                EXPECT_EQ(rotation.C, 1.0);
                EXPECT_EQ(rotation.S, 0.0);
                EXPECT_EQ(rotation.R, a);
            }
        }
    }
}
