#ifndef MY_DECOMPOSITION_H
#define MY_DECOMPOSITION_H

#include "matrix.hpp"

class Decomposition
{
public:

    template <class T>
    static std::pair<Matrix<double>, Matrix<double>> luDecomposition(const Matrix<T>& mat);
};


template <class T>
std::pair<Matrix<double>, Matrix<double>> Decomposition::luDecomposition(const Matrix<T>& mat)
{

}

#endif //MY_DECOMPOSITION_H

