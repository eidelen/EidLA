/****************************************************************************
** Copyright (c) 2019 Adrian Schneider
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

#ifndef MY_SAMPLESTATISTIC_H
#define MY_SAMPLESTATISTIC_H

#include "matrix.hpp"
#include "exceptions.hpp"

#include <vector>
#include <algorithm>


/**
 * @brief This class offers common statistical functions on a sample set.
 */
class SampleStatistic
{
public:

    /**
     * Compute the sample mean.
     * @param samples Vector of samples.
     * @return Sample mean.
     */
    static Matrix<double> mean(const std::vector<SamplePtr>& samples)
    {
        size_t n = samples.size();
        if(n == 0)
            throw EmptyContainerException();

        Matrix<double> accum = samples.at(0)->m_data;
        std::for_each(samples.begin()+1, samples.end(), [&accum](const SamplePtr& s){
            accum.add(s->m_data);
        });

        return accum * (1.0 / n);
    }


};

#endif // MY_SAMPLESTATISTIC_H