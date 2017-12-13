#ifndef MY_MULTIPLIER_H
#define MY_MULTIPLIER_H

#include "matrix.hpp"

class Multiplier
{
public:

    template <class T>
    static Matrix<T> swapRow(const Matrix<T>& mat, size_t r0, size_t r1);
};


template <class T>
Matrix<T> Multiplier::swapRow(const Matrix<T>& mat, size_t r0, size_t r1)
{
    if (std::max(r0, r1) >= mat.rows())
    {
        std::cout << "row index exceeds matrix size";
        std::exit(-1);
    }

    Matrix<T> swapOp(mat.rows(), mat.rows());
    swapOp.setToIdentity();

    Matrix<T> zeroRow(1,mat.rows());
    zeroRow.fill(0);

    // r0
    swapOp.setRow(r0,zeroRow);
    swapOp(r0,r1) = 1;

    // r1
    swapOp.setRow(r1,zeroRow);
    swapOp(r1,r0) = 1;

    return swapOp;
}

#endif //MY_MULTIPLIER_H

