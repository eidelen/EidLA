#include <gtest/gtest.h>
#include "matrix.hpp"
#include "sample.hpp"
#include "sample_generation.hpp"
#include "sample_statistic.hpp"

TEST(SampleStatistics, Mean)
{
    // no variation samples
    Matrix<double> mean(2, 1, {1.0, 3.0});
    std::vector<std::pair<double, Matrix<double>>> noVar;

    SampleGenerator* noVarDist = new NormalSampleGenerator(mean, noVar, 0);
    std::vector<SamplePtr> noVarSamples = noVarDist->getN(100);
    delete noVarDist;

    Matrix<double> compMean = SampleStatistic::mean(noVarSamples);
    ASSERT_TRUE(compMean.compare(mean, true, 0.00001));


    // mix with another sample vector
    Matrix<double> mean2(2, 1, {3.0, 1.0});

    SampleGenerator* noVarDist2 = new NormalSampleGenerator(mean2, noVar, 0);
    std::vector<SamplePtr> noVarSamples2 = noVarDist2->getN(100);
    delete noVarDist2;

    std::vector<SamplePtr> combinedSamples;
    std::move(noVarSamples.begin(), noVarSamples.end(), std::back_inserter(combinedSamples));
    std::move(noVarSamples2.begin(), noVarSamples2.end(), std::back_inserter(combinedSamples));

    ASSERT_EQ(combinedSamples.size(), 200);
    auto compMean2 = SampleStatistic::mean(combinedSamples);
    ASSERT_FLOAT_EQ(compMean2(0,0), 2.0);
    ASSERT_FLOAT_EQ(compMean2(1,0), 2.0);
}
