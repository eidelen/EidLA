#include <gtest/gtest.h>
#include "matrix.hpp"
#include "sample_generation.hpp"

TEST(SampleGen, Basic)
{
    Matrix<double> mat = Matrix<double>::random(10, 1, -100, +100);

    Sample s_labeled = Sample(mat, 33);
    ASSERT_EQ(s_labeled.m_label, 33);
    ASSERT_EQ(s_labeled.m_type, Sample::Labeled);
    ASSERT_TRUE(s_labeled.m_data.compare(mat, true, 0.0000001));

    Sample s_unlabeled = Sample(mat);
    ASSERT_EQ(s_unlabeled.m_type, Sample::Unlabeled);
    ASSERT_TRUE(s_unlabeled.m_data.compare(mat, true, 0.0000001));
}
