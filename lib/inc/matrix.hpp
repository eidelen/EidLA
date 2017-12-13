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

#ifndef MY_MATRIX_H
#define MY_MATRIX_H

#include <memory>
#include <iostream>
#include <cmath>
#include <tuple>
#include <limits>

#include <smmintrin.h>  // SSE4

template <class T>
class Matrix
{
public:

    /**
     * Constructor
     * @param rows number of rows
     * @param cols number of columns
     */
    Matrix( size_t rows, size_t cols );

    /**
     * Copy constructor.
     * @param mat Matrix from which to copy
     */
    template<class R>
    Matrix(const Matrix<R>& mat);

    Matrix(const Matrix<T>& mat);

    /**
     * Constructs a m x n matrix and assigns
     * the values taken from the data array
     * in row direction.
     * @param rows number of rows
     * @param cols number of columns
     * @param data Pointer to data
     */
    Matrix(size_t rows, size_t cols, const T* data);

    /**
     * Asignment operator: overwrite content of this matrix.
     * @param other
     * @return
     */
    Matrix<T>& operator=(const Matrix<T>& other);

    virtual ~Matrix();

    size_t rows() const;
    size_t cols() const;
    size_t getNbrOfElements() const;
    inline T*     data();
    inline const T* data() const;

    /**
     * Sets each element to the value val.
     * @param val Value
     */
    void fill(T val);

    /**
     * Set the value at position m, n.
     * @param m Row
     * @param n Column
     * @param val Value
     */
    inline void setValue( size_t m, size_t n, T val);

    /**
     * Gets the value at position m, n.
     * @param m Row
     * @param n Column
     * @return Value
     */
    inline T getValue( size_t m, size_t n ) const;

    /**
     * Returns the row at position m
     * @param m
     * @return row-vector of size 1xn
     */
    Matrix<T> row(size_t m) const;

    /**
     * Returns the column at position n
     * @param n
     * @return column-vector of size mx1
     */
    Matrix<T> column(size_t n) const;

    /**
     * Elementwise comparisson with the passed matrix mat.
     * @param epsilon The allowed tolerance.
     * @return True if all elements are within the tolerance. Otherwise false.
     */
    bool compare( const Matrix<T>& mat, double epsilon = std::numeric_limits<double>::min() ) const;

    /**
     * Gets the value at position m, n.
     * @param m
     * @param n
     * @return
     */
    const T operator() (size_t m, size_t n) const;
    T& operator() (size_t m, size_t n);

    /**
     * Multiplication of a matrix with a scalar.
     * @param scale Scalar
     * @return Matrix
     */
    Matrix<T> operator* (T scale) const;

    /**
     * Matrix multiplication
     * @param mat
     * @return Product of two matrix
     */
    Matrix<T> operator* (const Matrix<T>& mat) const;

    /**
     * Elementwise addition of two matrix.
     * @param mat
     * @return Resulting matrix.
     */
    Matrix<T> operator+ (const Matrix<T>& mat) const;

    /**
     * Elementwise subtraction of two matrix.
     * @param mat
     * @return Resulting matrix.
     */
    Matrix<T> operator- (const Matrix<T>& mat) const;

    /**
     * Creates a m x m identity matrix
     * @param m Matrix size.
     * @return Identity matrix.
     */
    static Matrix<T> identity(size_t m);

    /**
     * Set the elements of the matrix
     * to be identity.
     */
    void setToIdentity();

    /**
     * Returns the transpose of this matrix.
     * @return Matrix transpose.
     */
    Matrix<T> transpose() const;

    /**
     * Returns all elements on the central
     * diagonal as a column vector.
     * @return Diagonal as column vector.
     */
    Matrix<T> diagonal() const;

    /**
     * Sum of all elements in the matrix.
     * @return sum
     */
    T sum() const;

    /**
     * Get the position and value of the maximum element. If there
     * are two or more maximums, the first one is returned.
     * @return Position and value of maximum
     */
    std::tuple<size_t, size_t, T> max() const;

    /**
     * Get the position and value of the minimum element. If there
     * are two or more minimums, the first one is returned.
     * @return Position and value of minimum
     */
    std::tuple<size_t, size_t, T> min() const;

    /**
     * Swap rows m1 and m2
     * @param m1
     * @param m2
     */
    void swapRows( size_t m1, size_t m2 );

    /**
     * Sets the row rowIdx to values in row.
     * @param rowIdx
     * @param row
     */
    void setRow( size_t rowIdx, const Matrix<T>& row);

    /**
     * Compute the rank of the matrix, the number of
     * linearly independant rows.
     * http://stattrek.com/matrix-algebra/matrix-rank.aspx
     * @return Rank.
     */
    size_t getRank() const;

    /**
     * Returns the inverse of this matrix.
     * @param invertable On return, this is true if successfull
     * @return Inverse of this matrix
     */
    Matrix<double> inverted(bool* invertable) const;

protected:
    /**
     * Check if the dimensions of the two passed matrix are equal.
     * @param m1 Mat 1
     * @param m2 Mat 2
     * @return True if equal dimension.
     */
    static bool equalDimension( const Matrix<T> m1, const Matrix<T> m2);

    T elementwiseMultiplyAndSum(const T *arr1, const T *arr2, size_t length) const;

    T* getRowPtr(size_t row);

protected:
    size_t m_rows;
    size_t m_cols;

    std::shared_ptr<T> m_data;
    size_t m_nbrOfElements;
};


template <>
inline int Matrix<int>::elementwiseMultiplyAndSum(const int* arr1, const int* arr2, size_t length) const;

/*
template <>
inline double Matrix<double>::elementwiseMultiplyAndSum(const double* arr1, const double* arr2, size_t length) const;
 */


template <class T>
Matrix<T>::Matrix(size_t rows, size_t cols)
        : m_rows(rows), m_cols(cols), m_nbrOfElements(rows*cols)
{
    m_data.reset( new T[m_nbrOfElements] );
}

//** Specific copy transformations - to be extended
template <class T>
void copyMatData(const Matrix<T>& src, Matrix<T>& dst )
{
    const T* srcPtr = src.data();
    T* dstPtr = dst.data();
    std::copy(srcPtr, srcPtr+src.getNbrOfElements(), dstPtr);
}

inline void copyMatData(const Matrix<int>& src, Matrix<double>& dst )
{
    const int* srcPtr = src.data();
    double* dstPtr = dst.data();

    for( size_t i = 0; i < src.getNbrOfElements(); i++ )
        dstPtr[i] = double(srcPtr[i]);
}

//**

template <class T>
template <class R>
Matrix<T>::Matrix(const Matrix<R>& mat)
        : Matrix<T>(mat.rows(), mat.cols())
{
    // if you have a compile error here, you have to extend the
    // copyMatData function with the particular types.
    copyMatData(mat, *this);
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& mat)
        : Matrix<T>(mat.rows(), mat.cols())
{
    // if you have a compile error here, you have to extend the
    // copyMatData function with the particular types.
    copyMatData(mat, *this);
}

template <class T>
Matrix<T>::Matrix(size_t rows, size_t cols, const T* data)
        : Matrix<T>(rows, cols)
{
    T* dst = this->data();
    std::copy(data, data+this->getNbrOfElements(), this->data());
}

template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other)
{
    if (this != &other) // self-assignment check expected
    {
        // reinit matrix
        m_rows = other.rows();
        m_cols = other.cols();
        m_nbrOfElements = m_rows * m_cols;
        m_data.reset( new T[m_nbrOfElements] );

        // copy data
        copyMatData(other, *this);
    }

    return *this;
}

template <class T>
Matrix<T>::~Matrix()
{
}

template <class T>
size_t Matrix<T>::rows() const
{
    return m_rows;
}

template <class T>
inline T* Matrix<T>::data()
{
    return m_data.get();
}

template <class T>
inline const T* Matrix<T>::data() const
{
    return m_data.get();
}

template <class T>
size_t Matrix<T>::cols() const
{
    return m_cols;
}

template <class T>
size_t Matrix<T>::getNbrOfElements() const
{
    return m_nbrOfElements;
}

template <class T>
void Matrix<T>::fill(T val)
{
    T* dataPtr = data();
    for( size_t i = 0; i < m_nbrOfElements; i++ )
        dataPtr[i] = val;
}

template <class T>
inline void Matrix<T>::setValue( size_t m, size_t n, T val)
{
    data()[m*cols() + n] = val;
}

template <class T>
inline T Matrix<T>::getValue( size_t m, size_t n ) const
{
    return data()[m*cols() + n];
}

template <class T>
Matrix<T> Matrix<T>::row(size_t m) const
{
    if( m >= rows() )
    {
        std::cout << "Row access not possible. " << m << ", max " << rows();
        std::exit(-1);
    }

    Matrix<T> ret(1,cols());
    size_t offset = m * cols();
    std::copy(data()+offset, data()+offset+cols(), ret.data());

    return ret;
}

template <class T>
Matrix<T> Matrix<T>::column(size_t n) const
{
    if( n >= cols() )
    {
        std::cout << "Column access not possible. " << n << ", max " << cols();
        std::exit(-1);
    }

    Matrix<T> ret(rows(),1);
    for(size_t i = 0; i < rows(); i++)
        ret(i,0) = this->getValue(i,n);

    return ret;
}


template <class T>
const T Matrix<T>::operator() (size_t m, size_t n) const
{
    return data()[m*cols() + n];
}

template <class T>
T& Matrix<T>::operator() (size_t m, size_t n)
{
    return data()[m*cols() + n];
}

template <class T>
Matrix<T> Matrix<T>::identity(size_t m)
{
    Matrix<T> ident = Matrix<T>(m,m);
    ident.setToIdentity();
    return ident;
}

template <class T>
void Matrix<T>::setToIdentity()
{
    if( rows() != cols() )
    {
        std::cout << "Not a square matrix";
        std::exit(-1);
    }

    this->fill(0);
    for(size_t i = 0; i < rows(); i++)
    {
        this->setValue(i,i, 1);
    }
}

template <class T>
Matrix<T> Matrix<T>::transpose() const
{
    Matrix<T> ret( cols(), rows() );
    for( size_t m = 0; m < rows(); m++ )
    {
        for( size_t n = 0; n < cols(); n++ )
        {
            ret(n,m) = getValue(m,n);
        }
    }

    return ret;
}

template <class T>
Matrix<T> Matrix<T>::diagonal() const
{
    size_t diagLength = std::min(rows(), cols());
    Matrix<T> ret(diagLength, 1);

    for( size_t i = 0; i < diagLength; i++ )
    {
        ret(i,0) = getValue(i,i);
    }

    return ret;
}

template <class T>
T Matrix<T>::sum() const
{
    T ret = 0;

    const T* ptr = data();
    size_t nElem = getNbrOfElements();
    for( size_t i = 0; i < nElem; i++)
        ret += ptr[i];

    return ret;
}

template <class T>
bool Matrix<T>::equalDimension( const Matrix<T> m1, const Matrix<T> m2)
{
    return m1.cols()==m2.cols() && m1.rows()==m2.rows();
}

template <class T>
Matrix<T> Matrix<T>::operator* (T scale) const
{
    Matrix<T> res(this->rows(), this->cols());
    T* resD = res.data();
    const T* dataD = this->data();
    for(size_t i = 0; i < this->getNbrOfElements(); i++)
        resD[i] = dataD[i] * scale;

    return res;
}

template <class T>
Matrix<T> operator*(T scale, const Matrix<T>& mat)
{
    return mat * scale;
}

template <class T>
Matrix<T> Matrix<T>::operator* (const Matrix<T>& mat) const
{
    if( this->cols() != mat.rows() )
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    Matrix<T> res(this->rows(), mat.cols());

    //Create a lookup table of the right matrix split to column vectors
    // -> this is helpful because this way the column vector is in a
    //    continous memory block.
    std::vector<Matrix<T>> lookup;
    for(size_t k = 0; k < mat.cols(); k++)
        lookup.push_back(mat.column(k));

    for( size_t n = 0; n < mat.cols(); n++)
    {
        for( size_t m = 0; m < this->rows(); m++)
        {
            // get column mat -> continous memory block
            const T* matColNPtr = lookup.at(n).data();
            const T* thisRowPtr = data()+m*cols(); // this is fast

            res(m,n) = elementwiseMultiplyAndSum(thisRowPtr, matColNPtr, this->cols() );
        }
    }

    return res;
}

template <class T>
Matrix<T> Matrix<T>::operator+ (const Matrix<T>& mat) const
{
    // check dimension
    if(!equalDimension(*this, mat))
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    Matrix<T> res = Matrix<T>(rows(),cols());
    T* resD = res.data(); const T* thisData = data(); const T* matD = mat.data();
    for(size_t i = 0; i < getNbrOfElements(); i++)
        resD[i] = thisData[i] + matD[i];

    return res;
}

template <class T>
Matrix<T> Matrix<T>::operator- (const Matrix<T>& mat) const
{
    // check dimension
    if(!equalDimension(*this, mat))
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    Matrix<T> res = Matrix<T>(rows(),cols());
    T* resD = res.data(); const T* thisData = data(); const T* matD = mat.data();
    for(size_t i = 0; i < getNbrOfElements(); i++)
        resD[i] = thisData[i] - matD[i];

    return res;
}


template <class T>
std::ostream& operator<<(std::ostream& os, Matrix<T> const& mat)
{
    for( size_t m = 0; m < mat.rows(); m++ )
    {
        for( size_t n = 0; n < mat.cols(); n++ )
        {
            os << mat(m,n);
            if( n+1 < mat.cols() )
                os << ", ";
        }

        os << std::endl;
    }

    return os;
}

template <class T>
bool Matrix<T>::compare(const Matrix<T>& mat, double epsilon) const
{
    if(!equalDimension(*this, mat))
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    const T* thisData = data(); const T* matD = mat.data();
    for(size_t i = 0; i < getNbrOfElements(); i++)
    {
        double diff = thisData[i] - matD[i];
        if( std::fabs(diff) > epsilon )
            return false;
    }

    return true;

}

template <class T>
T Matrix<T>::elementwiseMultiplyAndSum(const T *arr1, const T *arr2, size_t length) const
{
    T accum = 0;
    for( size_t a = 0; a < length; a++ )
        accum += arr1[a] * arr2[a];

    return accum;
}

template <>
inline int Matrix<int>::elementwiseMultiplyAndSum(const int* arr1, const int* arr2, size_t length) const
{
    __m128i* arr_1_ptr = (__m128i*) arr1;
    __m128i* arr_2_ptr = (__m128i*) arr2;

    int accum = 0;
    size_t simdLoops = length / 4;
    size_t singleLoops = length - 4*simdLoops;

    for( size_t simd = 0; simd < simdLoops; simd++ )
    {
        __m128i reg_1_SSE = _mm_load_si128(arr_1_ptr);
        __m128i reg_2_SSE = _mm_load_si128(arr_2_ptr);
        __m128i mulRes = _mm_mullo_epi32 (reg_1_SSE, reg_2_SSE);
        __m128i sumRes = _mm_hadd_epi32(mulRes,mulRes);
                sumRes = _mm_hadd_epi32(sumRes,sumRes);

        accum += _mm_extract_epi32(sumRes, 0);

        arr_1_ptr++;
        arr_2_ptr++;
    }

    const int* offset1 = arr1 + 4*simdLoops;
    const int* offset2 = arr2 + 4*simdLoops;
    for(size_t a = 0; a < singleLoops; a++)
    {
        accum += offset1[a]* offset2[a];
    }

    return accum;
}

/* Makes the multiplication slower
template <>
inline double Matrix<double>::elementwiseMultiplyAndSum(const double* arr1, const double* arr2, size_t length) const
{
    const double* arr_1_ptr = arr1;
    const double* arr_2_ptr = arr2;

    int accum = 0;
    size_t simdLoops = length / 2;
    size_t singleLoops = length - 2*simdLoops;

    const int mask = 0x31;
    for( size_t simd = 0; simd < simdLoops; simd++ )
    {
        __m128d reg_1_SSE = _mm_load_pd(arr_1_ptr);
        __m128d reg_2_SSE = _mm_load_pd(arr_2_ptr);
        __m128d res = _mm_dp_pd(reg_1_SSE, reg_2_SSE, mask);

        accum += _mm_cvtsd_f64(res);

        arr_1_ptr+=2;
        arr_2_ptr+=2;
    }

    const double* offset1 = arr1 + 2*simdLoops;
    const double* offset2 = arr2 + 2*simdLoops;
    for(size_t a = 0; a < singleLoops; a++)
    {
        accum += offset1[a]* offset2[a];
    }

    return accum;
}
*/


template <class T>
std::tuple<size_t, size_t, T> Matrix<T>::max() const
{
    T maxVal = std::numeric_limits<T>::min();
    size_t maxPos = 0;
    const T* matData = data();

    for( size_t i = 0; i < getNbrOfElements(); i++ )
    {
        if( maxVal <  matData[i] )
        {
            maxPos = i;
            maxVal = matData[i];
        }
    }

    size_t m = maxPos / cols();
    size_t n = maxPos - (m*cols());

    return std::make_tuple(m,n,maxVal);
}

template <class T>
std::tuple<size_t, size_t, T> Matrix<T>::min() const
{
    T minVal = std::numeric_limits<T>::max();
    size_t minPos = 0;
    const T* matData = data();

    for( size_t i = 0; i < getNbrOfElements(); i++ )
    {
        if( minVal >  matData[i] )
        {
            minPos = i;
            minVal = matData[i];
        }
    }

    size_t m = minPos / cols();
    size_t n = minPos - (m*cols());

    return std::make_tuple(m,n,minVal);
}

template <class T>
T* Matrix<T>::getRowPtr(size_t row)
{
    return data() + row*cols();
}

template <class T>
void Matrix<T>::swapRows(size_t m1, size_t m2)
{
    if (std::max(m1, m2) >= rows())
    {
        std::cout << "row index exceeds matrix size";
        std::exit(-1);
    }

    T* row1 = getRowPtr(m1);
    T* row2 = getRowPtr(m2);

    for (size_t i = 0; i < cols(); i++)
    {
        std::swap(row1[i],row2[i]);
    }
}

template <class T>
void Matrix<T>::setRow(size_t rowIdx, const Matrix<T>& row)
{
    if (rowIdx >= rows())
    {
        std::cout << "row index exceeds matrix size";
        std::exit(-1);
    }

    if( row.cols() != cols() )
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    T* dstPtr = getRowPtr(rowIdx);
    const T* srcPtr = row.data();
    for( size_t i = 0; i < row.cols(); i++ )
        dstPtr[i] = srcPtr[i];

}

// Predefined Matrix Types
typedef std::shared_ptr<Matrix<int>> MatrixISP;
typedef std::shared_ptr<Matrix<double>> MatrixDSP;


#include "transformation.hpp"

template <class T>
size_t Matrix<T>::getRank() const
{
    Matrix<double> echelonMat = Transformation::echelon(*this);

    // count numbers of non zero rows
    Matrix<double> zeroRow(1, cols());
    zeroRow.fill(0.0);

    size_t rank = 0;
    for(size_t i = 0; i < rows(); i++)
    {
        Matrix<double> cRow = echelonMat.row(i);
        if( cRow.compare(zeroRow, std::numeric_limits<double>::min()) )
            break;
        else
            rank++;
    }

    return rank;
}

template <class T>
Matrix<double> Matrix<T>::inverted(bool* invertable) const
{
    // needs to be square matrix
    if( m_rows != m_cols )
    {
        std::cout << "Needs to be square matrix" << std::endl;
        *invertable = false;
        return *this;
    }

    // todo: check that determinant not = 0

    *invertable = true;

    std::vector<Matrix<double>> ops;
    Transformation::reduced_echelon(*this, ops);

    Matrix<double> inv = Matrix<double>(m_rows, m_cols);
    inv.setToIdentity();

    // multiply Echelon operators in reverse order
    // http://stattrek.com/matrix-algebra/how-to-find-inverse.aspx
    while( !ops.empty() )
    {
        inv = inv*ops.back();
        ops.pop_back();
    }

    return inv;
}



#endif //MY_MATRIX_H