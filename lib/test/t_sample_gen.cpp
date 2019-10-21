#include <gtest/gtest.h>
#include "matrix.hpp"
#include "sample_generation.hpp"

TEST(SampleGen, SampleStorage)
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

TEST(SampleGen, UniformCtor)
{
    Matrix<double> mean_invalid(2,2);
    Matrix<double> mean(2,1);

    std::vector<std::pair<double,Matrix<double>>> dist_invalid;
    dist_invalid.emplace_back(0, Matrix<double>(4,1));

    std::vector<std::pair<double,Matrix<double>>> dist;
    dist.emplace_back(0, Matrix<double>(2,1));

    {
        ASSERT_ANY_THROW( SampleGenerator *uniform = new UniformSampleGenerator(mean_invalid, dist, 0 ) );
    }

    {
        ASSERT_ANY_THROW( SampleGenerator *uniform = new UniformSampleGenerator(mean, dist_invalid, 0 ) );
    }

    {
        ASSERT_NO_THROW( SampleGenerator *uniform = new UniformSampleGenerator(mean, dist, 0 ) );
    }
}

TEST(SampleGen, UniformNoVariation)
{
    Matrix<double> mean(2, 1, {0.3, 0.7});
    std::vector<std::pair<double, Matrix<double>>> noVar;

    SampleGenerator* noVarDist = new UniformSampleGenerator(mean, noVar, 999);

    // always mean returned, as there is not variation
    for(int k = 0; k < 100; k++)
    {
        SamplePtr sample = noVarDist->get();

        ASSERT_TRUE(mean.compare(sample->m_data, 0.0));
        ASSERT_EQ(sample->m_label, 999);
        ASSERT_EQ(sample->m_type, Sample::Labeled);
    }

    delete noVarDist;
}

TEST(SampleGen, UniformOneVariation)
{
    Matrix<double> mean(2, 1, {0.0, 0.0});
    std::vector<std::pair<double, Matrix<double>>> var;

    double direction = 1.0 / sqrt(2.0);
    Matrix<double> xy(2, 1, {direction, direction});
    var.push_back({1.0, xy});

    SampleGenerator* varDist = new UniformSampleGenerator(mean, var, 111);

    double max = -1;
    double min = +1;
    double meanSamples = 0.0;

    int n = 1000;
    for(int k = 0; k < n; k++)
    {
        SamplePtr sample = varDist->get();

        double x = sample->m_data(0,0);
        double y = sample->m_data(1,0);

        meanSamples += x;
        if( x < min )
            min = x;
        if( x > max )
            max = x;

        //xy are same, as variation goes into direction 1,1
        ASSERT_NEAR(x, y, 0.00001);
        ASSERT_EQ(sample->m_label, 111);
        ASSERT_EQ(sample->m_type, Sample::Labeled);
    }

    ASSERT_NEAR(min, -direction, 0.05);
    ASSERT_NEAR(max, direction, 0.05);

    // mean should be close to 0.0
    double compMean = meanSamples / n;
    ASSERT_NEAR(compMean, 0.0, 0.05);

    delete varDist;
}

TEST(SampleGen, UniformMultipleVariation)
{
    Matrix<double> mean(2, 1, {0.0, 0.0});
    std::vector<std::pair<double, Matrix<double>>> var;

    double direction = 1.0 / sqrt(2.0);
    Matrix<double> xy(2, 1, {direction, direction});
    Matrix<double> nxy(2, 1, {direction, -direction});
    var.push_back({1.0, xy});
    var.push_back({1.0, nxy});

    SampleGenerator* varDist = new UniformSampleGenerator(mean, var, 111);

    double max = -1;
    double min = +1;
    double meanSamples = 0.0;

    int n = 1000;
    for(int k = 0; k < n; k++)
    {
        SamplePtr sample = varDist->get();

        double x = sample->m_data(0,0);
        double y = sample->m_data(1,0);

        ASSERT_GE(x, -2*direction);
        ASSERT_GE(y, -2*direction);
        ASSERT_LE(x, 2*direction);
        ASSERT_LE(y, 2*direction);

        ASSERT_EQ(sample->m_label, 111);
        ASSERT_EQ(sample->m_type, Sample::Labeled);
    }

    // mean should be close to 0.0
    double compMean = meanSamples / n;
    ASSERT_NEAR(compMean, 0.0, 0.05);

    delete varDist;
}

TEST(SampleGen, NormalOneVariation)
{
    Matrix<double> mean(2, 1, {0.0, 0.0});
    std::vector<std::pair<double, Matrix<double>>> var;

    double direction = 1.0 / sqrt(2.0);
    Matrix<double> xy(2, 1, {direction, direction});
    var.push_back({1.0, xy});

    SampleGenerator* varDist = new NormalSampleGenerator(mean, var, 111);

    double meanSamples = 0.0;
    std::vector<int> distEval(10);

    int n = 1000;
    for(int k = 0; k < n; k++)
    {
        SamplePtr sample = varDist->get();

        double x = sample->m_data(0,0);
        double y = sample->m_data(1,0);

        meanSamples += x;

        //check segment from -2 to +2
        for(size_t i = 0; i < 10; i++)
        {
            double lowerSeg = -2.0 + i*0.4;
            double upperSeg = lowerSeg + 0.4;

            if( x > lowerSeg && x < upperSeg )
            {
                distEval.at(i) = distEval.at(i) + 1;
                break;
            }
        }

        //xy are same, as variation goes into direction 1,1
        ASSERT_NEAR(x, y, 0.00001);
        ASSERT_EQ(sample->m_label, 111);
        ASSERT_EQ(sample->m_type, Sample::Labeled);
    }

    // mean should be close to 0.0
    double compMean = meanSamples / n;
    ASSERT_NEAR(compMean, 0.0, 0.05);

    delete varDist;

    // test rising
    for(size_t i = 0; i < 4; i++)
    {
        ASSERT_GT(distEval.at(i+1), distEval.at(i));
    }

    // test falling
    for(size_t i = 5; i < 9; i++)
    {
        ASSERT_LT(distEval.at(i+1), distEval.at(i));
    }

    for( int d : distEval )
    {
        std::cout << std::string(d/5, '*') << std::endl;
    }
}
