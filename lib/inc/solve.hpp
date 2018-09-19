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

#ifndef MY_SOLVE_H
#define MY_SOLVE_H

#include "matrix.hpp"
#include "transformation.hpp"

class Solve
{
public:
    /**
    * Solves linear system of equations.
    * @param mat Matrix of coefficients.
    * @param b Condition vector.
    * @return Solution vector.
    */
    template <class T>
    static Matrix<double> solve_lseq(const Matrix<T>& mat, const Matrix<T>& b);
};

template <class T>
Matrix<double> Solve::solve_lseq(const Matrix<T>& mat, const Matrix<T>& b)
{
    // check input
    if (b.cols() != 1)
    {
        std::cout << "Error: Condition vector wrong dimension";
        std::exit(-1);
    }
    else if (b.rows() != mat.cols())
    {
        std::cout << "Error: Mismatching condition vector and coefficient matrix";
        std::exit(-1);
    }
    else if (!mat.isSquare())
    {
        std::cout << "Error: Coefficient matrix needs to be square matrix";
        std::exit(-1);
    }

    // make augmented matrix
    Matrix<double> a(mat.rows(), mat.cols() + 1);
    a.setSubMatrix(0, 0, mat);
    a.setSubMatrix(0, mat.cols(), b);

    // solve by using full-pivoting reduced echelon transformation
    std::vector<Matrix<double>> rowOps; // not used
    Matrix<double>              redEch = Transformation::reduced_echelon(a, rowOps, true);

    // solution vector is in last column
    return redEch.column(mat.cols());
}

#endif //MY_SOLVE_H
