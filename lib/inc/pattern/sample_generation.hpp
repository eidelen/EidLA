/****************************************************************************
** Copyright (c) 2017 Adrian Schneider
**
** Permission is hereby granted, free of charge, to any person obtaining a
** copy of this software and associated documentation files (the "Software"),
** to deal in the Software without restriction, including without limitation
** the rights to use, copy, modify, merge, publish, distribute, sublicense,
** and/or sell copies of the Software, and to permit persons to whom the
** Software is furnished to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in
** all copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
** FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
** DEALINGS IN THE SOFTWARE.
**
*****************************************************************************/

#ifndef MY_SAMPLEGENERATION_H
#define MY_SAMPLEGENERATION_H

#include "matrix.hpp"
#include "exceptions.hpp"

#include <vector>
#include <memory>
#include <random>


/**
 * This class holds one sample.
 */
class Sample
{
public:
    enum Type
    {
        Labeled,    ///< Labeled sample
        Unlabeled   ///< Unlabeled sample
    };

public:

    /**
     * Construct a sample with given data and label.
     * @param data Sample data
     * @param label Sample label
     */
    Sample(const Matrix<double>& data, int label) : m_data(data), m_label(label), m_type(Labeled)
    {}

    /**
     * Construct a unlabeled sample.
     * @param data Sample data.
     */
    Sample(const Matrix<double>& data) : m_data(data), m_type(Unlabeled)
    {}

    Matrix<double> m_data;
    int m_label = 0;
    Type m_type;
};

typedef std::shared_ptr<Sample> SamplePtr;

class SampleGenerator
{
public:
    SampleGenerator(int label) : m_label(label), m_type(Sample::Labeled) {}
    SampleGenerator() : m_type(Sample::Unlabeled) {}
    virtual ~SampleGenerator() {}

    /**
     * Generates one sample.
     * @return A sample
     */
    virtual SamplePtr get() = 0;

    /**
     * Generates a vector of samples
     * @param nbr Number of items.
     * @return Sample vector.
     */
    std::vector<SamplePtr> getN(size_t nbr)
    {
        std::vector<SamplePtr> sampleVect(nbr);
        std::generate(sampleVect.begin(), sampleVect.end(), [=] () mutable {
            return get(); });
        return sampleVect;
    }

protected:

    SamplePtr makeSample(const Matrix<double>& data) const
    {
        if( m_type == Sample::Unlabeled )
            return std::make_shared<Sample>(data);
        else
            return std::make_shared<Sample>(data, m_label);
    }

private:
    int m_label;
    Sample::Type m_type;
};


class DistributionSampleGenerator: public SampleGenerator
{
public:

    /**
     * Constructor of distribution sample generator.
     * @param mean Sample mean.
     * @param distAxes Distribution axes. One axes consists of a normalized direction vector and its corresponding range scale.
     * @param label Sample label.
     */
    DistributionSampleGenerator(const Matrix<double>& mean, const std::vector<std::pair<double,Matrix<double>>>& distAxes, int label)
            : SampleGenerator::SampleGenerator(label), m_mean(mean)
    {
        // check dist vectors
        for(auto d : distAxes)
        {
            if( Matrix<double>::equalDimension(d.second, m_mean) && d.second.cols() == 1 )
                m_distri.push_back(d);
            else
                throw InvalidInputException();
        }

        m_randomEngine = std::default_random_engine();
    }

    virtual ~DistributionSampleGenerator() {}

    /**
     * Generates a sample by taking the mean value and adding all given random uniform distribution vectors.
     * @return A sample.
     */
    SamplePtr get() override
    {
        Matrix<double> sample(m_mean); // starting from mean
        for(auto d: m_distri)
        {
            double randomValue = getRandomValue();
            sample = sample + ( d.second * (d.first*randomValue) );
        }

        return makeSample(sample);
    }

protected:

    /**
     * Get a random number. This function has to be implemented in order to follow
     * a certain distribution.
     * @return Random value
     */
    virtual double getRandomValue() = 0;

    Matrix<double> m_mean;
    std::vector<std::pair<double,Matrix<double>>> m_distri;
    std::default_random_engine m_randomEngine;
};


/**
 * Uniform sample generator.
 */
class UniformSampleGenerator: public DistributionSampleGenerator
{
public:

    /**
     * Constructor of uniform sample generator.
     * @param mean Sample mean.
     * @param distAxes Uniform distribution axes. One axes consists of a normalized direction vector and its corresponding range scale.
     * @param label Sample label.
     */
    UniformSampleGenerator(const Matrix<double>& mean, const std::vector<std::pair<double,Matrix<double>>>& distAxes, int label)
    : DistributionSampleGenerator::DistributionSampleGenerator(mean, distAxes, label)
    {
        m_uniform = std::uniform_real_distribution<double>(-1.0, 1.0);
    }


protected:

    double getRandomValue() override
    {
        return m_uniform(m_randomEngine);
    }

    std::uniform_real_distribution<double> m_uniform;
};


/**
 * Normal distribution sample generator.
 */
class NormalSampleGenerator: public DistributionSampleGenerator
{
public:

    /**
     * Constructor of normal sample generator.
     * @param mean Sample mean.
     * @param distAxes Normal distribution axes. One axes consists of a normalized direction vector and its corresponding range scale.
     * @param label Sample label.
     */
    NormalSampleGenerator(const Matrix<double>& mean, const std::vector<std::pair<double,Matrix<double>>>& distAxes, int label)
            : DistributionSampleGenerator::DistributionSampleGenerator(mean, distAxes, label)
    {
        m_normal = std::normal_distribution<double>(0.0, 1.0);
    }

protected:

    double getRandomValue() override
    {
        return m_normal(m_randomEngine);
    }

    std::normal_distribution<double> m_normal;
};


#endif //MY_SAMPLEGENERATION_H
