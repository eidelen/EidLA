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

#ifndef MY_SAMPLE_H
#define MY_SAMPLE_H

#include "matrix.hpp"
#include <memory>


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

#endif //MY_SAMPLE_H
