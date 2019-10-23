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

}
