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

#include "exceptions.hpp"

#include <memory>
#include <iostream>
#include <cmath>
#include <tuple>
#include <limits>
#include <random>
#include <fstream>
#include <algorithm>
#include <cstring>

#include <smmintrin.h> // SSE4

#ifdef OPENCVEIDLA
#include <opencv2/core/core.hpp>
#endif // OPENCVEIDLA

template <class T>
class Matrix
{
public:

    Matrix();

    /**
     * Constructs a matrix of size rows x cols and sets
     * all values to 0.
     * @param rows number of rows
     * @param cols number of columns
     */
    Matrix(size_t rows, size_t cols);

    /**
     * Copy constructor.
     * @param mat Matrix from which to copy
     */
    template <class R>
    Matrix(const Matrix<R>& mat);

    Matrix(const Matrix<T>& mat);

    /**
     * Move constructor
     * @param other
     */
    Matrix(Matrix&& other);

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
     * Constructs a matrix based on the specified
     * file at filepath. See the functions load
     * and save.
     * @param filepath
     */
    Matrix(const std::string& filepath);


#ifdef OPENCVEIDLA
    /**
     * Constructs a matrix from an OpenCV matrix.
     * @param mat OpenCV matrix
     */
    Matrix( const cv::Mat& mat);
#endif // OPENCVEIDLA

    /**
     * Asignment operator: overwrite content of this matrix.
     * @param other
     * @return
     */
    Matrix<T>& operator=(const Matrix<T>& other);

    Matrix<T>& operator=(Matrix<T>&& other);

    virtual ~Matrix();

    size_t          rows() const;
    size_t          cols() const;
    size_t          getNbrOfElements() const;
    inline T*       data();
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
    inline void setValue(size_t m, size_t n, T val);

    /**
     * Gets the value at position m, n.
     * @param m Row
     * @param n Column
     * @return Value
     */
    inline T getValue(size_t m, size_t n) const;

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
     * Removes the whole row at index m
     * @param m Index of row to remove.
     */
    void removeRow(size_t m);

    /**
     * Removes the whole column at index n
     * @param n Index of column to remove.
     */
    void removeColumn(size_t n);

    /**
     * Create a nbrOfRows x nbrOfColumns submatrix starting at rowStart/ columnStart.
     * @param rowStart Submatrix starts at rowStart.
     * @param columnStart Submatrix starts at columnStart.
     * @param nbrOfRows
     * @param nbrOfColumns
     * @return
     */
    Matrix<T> subMatrix(size_t rowStart, size_t columnStart, size_t nbrOfRows, size_t nbrOfColumns) const;

    /**
     * Writes the elements of the submatrix to this matrix at the given position.
     * @param rowStart Submatrix starts at rowStart.
     * @param columnStart Submatrix starts at columnStart.
     * @param subMat Matrix to write
     */
    template <class R>
    void setSubMatrix(size_t rowStart, size_t columnStart, const Matrix<R>& subMat);

    /**
     * Elementwise comparisson with the passed matrix mat.
     * @param mat Matrix to compare to.
     * @param useCustomTolerance Set if custom tolerance is used or not. Default false.
     * @param customTolerance Custom tolerance.
     * @return True if all elements are within the tolerance. Otherwise false.
     */
    bool compare(const Matrix<T>& mat, bool useCustomTolerance = false, T customTolerance = 0) const;

    /**
     * Gets the value at position m, n.
     * @param m
     * @param n
     * @return
     */
    const T operator()(size_t m, size_t n) const;
    T& operator()(size_t m, size_t n);

    /**
     * Multiplication of a matrix with a scalar.
     * @param scale Scalar
     * @return Matrix
     */
    Matrix<T> operator*(T scale) const;

    /**
     * Matrix multiplication
     * @param mat
     * @return Product of two matrix "this * mat"
     */
    Matrix<T> operator*(const Matrix<T>& mat) const;

    /**
     * Right side matrix multiplication.
     * @param mat
     * @return Product of two matrix "this * mat"
     */
    Matrix<T> matMulR(const Matrix<T>& mat) const;

    /**
     * Elementwise division
     * @param mat
     * @return Resulting matrix
     */
    Matrix<T> operator/(const Matrix<T>& mat) const;

    /**
     * Elementwise addition of two matrix.
     * @param mat
     * @return Resulting matrix.
     */
    Matrix<T> operator+(const Matrix<T>& mat) const;

    /**
     * Elementwise subtraction of two matrix.
     * @param mat
     * @return Resulting matrix.
     */
    Matrix<T> operator-(const Matrix<T>& mat) const;

    /**
     * Creates a matrix of the size m x n filled with random
     * values in the range of lower to upper.
     * @param m
     * @param
     * @param lower
     * @param upper
     * @return Matrix filled with random values
     */
    static Matrix<T> random(size_t m, size_t n, T lower, T upper);

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
     * Sum up elements in column direction.
     * @return rows x 1 matrix
     */
    Matrix<T> sumC() const;

    /**
     * Sum up elements in row direction.
     * @return 1 x columns matrix
     */
    Matrix<T> sumR() const;

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
     * Computes the mean of all matrix elements.
     * @return Mean value.
     */
    double mean() const;

    /**
     * Swap rows m1 and m2
     * @param m1
     * @param m2
     */
    void swapRows(size_t m1, size_t m2);

    /**
     * Sets the row rowIdx to values from row.
     * @param rowIdx
     * @param row
     */
    void setRow(size_t rowIdx, const Matrix<T>& row);

    /**
     * Sets the column colIdx to values from col.
     * @param colIdx
     * @param col
     */
    void setColumn(size_t colIdx, const Matrix<T>& col);

    /**
     * Compute the rank of the matrix, the number of
     * linearly independant rows.
     * http://stattrek.com/matrix-algebra/matrix-rank.aspx
     * @return Rank.
     */
    size_t getRank() const;

    /**
     * Returns the inverse of this matrix.
     * @param invertable On return, this is true if successful
     * @return Inverse of this matrix
     */
    Matrix<double> inverted(bool* invertable) const;

    /**
     * Returns the determinant of this matrix.
     * @param successful On return, this is true if successful
     * @return Determinant of this matrix
     */
    double determinant(bool* successful) const;

    /**
     * Returns a matrix consisting of normalized column vectors.
     * @return
     */
    Matrix<double> normalizeColumns() const;

    /**
     * Normalizes the matrix so that the elements are in
     * the range of -1 to +1.
     * @param mean On return, the applied mean.
     * @param scale On return, the applied scale.
     * @return Normalized Matrix
     */
    Matrix<double> normalize(double& mean, double& scale) const;

    /**
     * Denormalizes the matrix.
     * @param mean The mean
     * @param scale The scale
     * @return Denormalized matrix
     */
    Matrix<double> denormalize(double mean, double scale) const;


    /**
     * Returns the first minors of this matrix.
     * @return First minors matrix.
     */
    Matrix<double> firstMinors() const;

    /**
     * Return the cofactors of this matrix.
     * @return Cofactors matrix.
     */
    Matrix<double> cofactors() const;

    /**
     * Return the adjugate (adjoint) of this matrix.
     * @return Adjugate matrix.
     */
    Matrix<double> adjugate() const;

    /**
     * Computes the Euclidean norm of a mx1 or 1xn vector.
     * @return Euclidean norm (2-norm)
     */
    double norm() const;

    /**
     * Computes the square of the Euclidean vector norm
     * of a mx1 or 1xn vector. The difference to the
     * function norm() is that normSquare() is faster.
     * @return Euclidean norm (2-norm)
     */
    T normSquare() const;

    /**
     * Checks if this matris is symmetric: A = A'
     * @return True if symmetric. Otherwise false.
     */
    bool isSymmetric() const;

    /**
     * Checks if this matrix is a square matrix: NbrRows = NbrColumns
     * @return True if square. Otherwise false.
     */
    bool isSquare() const;

    /**
     * Checks if this matrix is an orthogonal matrix.
     * @param customTolerance Allowed tolerance of norm.
     * @return True if orthogonal. Otherwise false.
     */
    bool isOrthogonal(T customTolerance = 0) const;

    /**
     * Save matrix at the passed path.
     * @param path Path and filename.
     * @return True if successful. Otherwise false.
     */
    bool save(const std::string& path) const;

    /**
     * Load matrix from passed path.
     * @param path Path anf filename.
     * @return True if successful. Otherwise false.
     */
    bool load(const std::string& path);

    /**
     * Serialize the matrix as a data
     * array.
     * @return Data array in form of a string.
     */
    std::string serialize() const;

    /**
     * Deserialize the passed data into
     * this matrix.
     * @param data Data array in form of a string.
     */
    void deserialize(const std::string& data);

    /**
     * Sets matrix elements to zero, which are smaller
     * than the passed threshold:
     * |value| < threshold   ->  value = 0;
     * @param threshold
     */
    void setToZero( T threshold );

    /**
     * Computes the L1 norm of the matrix or vector.
     * @return L1 norm
     */
    T normL1() const;

    /**
     * Computes the infinity norm of the matrix or vector.
     * @return Infinity norm
     */
    T normInf() const;

    /**
     * Computes the Euclidean length of a vector. If this
     * is a matrix, it computes the maximum of the matrix.
     * @return Infinity norm
     */
    double normL2() const;

    /**
     * Computes the L1 norm condition number.
     * @return L1 condition number
     */
    double conditionNumberL1() const;

    /**
     * Computes the infinity norm condition number.
     * @return Infinity condition number
     */
    double conditionNumberInf() const;

    /**
     * Returns a matrix with absolute
     * values (all values are positive).
     * @return Absolute matrix.
     */
    Matrix<T> absolute() const;

    /**
     * Creates a repetition of this matrix.
     * @param rowDirection 1 - m repetitions in row direction.
     * @param columnDirection 1 - n repetitions in column direction.
     * @return Repeated matrix
     */
    Matrix<T> repMat(size_t rowDirection, size_t columnDirection) const;

#ifdef OPENCVEIDLA
    /**
     * Returns an OpenCV matrix of this
     * matrix.
     * @return OpendCV matrix.
     */
    cv::Mat toOpenCVMatrix( ) const;
#endif // OPENCVEIDLA

protected:

    /**
     * Check if the dimensions of the two passed matrix are equal.
     * @param m1 Mat 1
     * @param m2 Mat 2
     * @return True if equal dimension.
     */
    static bool equalDimension(const Matrix<T> m1, const Matrix<T> m2);

    T elementwiseMultiplyAndSum(const T* arr1, const T* arr2, size_t length) const;

    T* getRowPtr(size_t row);
    const T* getRowPtr(size_t row) const;

#ifdef OPENCVEIDLA
    cv::Mat createOpenCVMat( ) const;
#endif // OPENCVEIDLA

protected:
    size_t m_rows;
    size_t m_cols;

    std::shared_ptr<T> m_data;
    size_t             m_nbrOfElements;
};


#ifdef OPENCVEIDLA
template <>
inline cv::Mat Matrix<double>::createOpenCVMat( ) const;
template <>
inline cv::Mat Matrix<float>::createOpenCVMat( ) const;
template <>
inline cv::Mat Matrix<int>::createOpenCVMat( ) const;
#endif // OPENCVEIDLA

template <>
inline int Matrix<int>::elementwiseMultiplyAndSum(const int* arr1, const int* arr2, size_t length) const;

/*
template <>
inline double Matrix<double>::elementwiseMultiplyAndSum(const double* arr1, const double* arr2, size_t length) const;
 */

template <class T>
Matrix<T>::Matrix()
: m_rows(0), m_cols(0), m_nbrOfElements(0)
{
    m_data.reset();
}

template <class T>
Matrix<T>::Matrix(size_t rows, size_t cols)
: m_rows(rows), m_cols(cols), m_nbrOfElements(rows * cols)
{
    m_data.reset(new T[m_nbrOfElements]);
}


template <class T, class R>
void copyMatData(const Matrix<T>& src, Matrix<R>& dst)
{
    const T* srcPtr = src.data();
    R*       dstPtr = dst.data();
    std::copy(srcPtr, srcPtr + src.getNbrOfElements(), dstPtr);
}

template <class T>
void copyMatData(const Matrix<T>& src, Matrix<T>& dst)
{
    const T* srcPtr = src.data();
    T*       dstPtr = dst.data();
    std::memcpy(dstPtr, srcPtr, dst.getNbrOfElements() * sizeof(T));
}

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
Matrix<T>::Matrix(Matrix&& other)
{
    this->m_data          = std::move( other.m_data );
    this->m_rows          = other.rows();
    this->m_cols          = other.cols();
    this->m_nbrOfElements = m_rows * m_cols;
}

template <class T>
Matrix<T>::Matrix(size_t rows, size_t cols, const T* data)
: Matrix<T>(rows, cols)
{
    T* dst = this->data();
    std::copy(data, data + this->getNbrOfElements(), dst);
}

template <class T>
Matrix<T>::Matrix(const std::string& filepath)
{
    // allocates and overwrites
    load(filepath);
}

#ifdef OPENCVEIDLA
template <class T>
Matrix<T>::Matrix(const cv::Mat& mat)
: Matrix<T>(mat.rows, mat.cols)
{
    for(size_t m = 0; m < mat.rows; m++)
    {
        const T* src = mat.ptr<T>(m); // pointer to row m
        std::copy(src, src + mat.cols, this->getRowPtr(m));
    }
}
#endif // OPENCVEIDLA

template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other)
{
    if (this != &other) // self-assignment check expected
    {
        if (m_nbrOfElements != other.getNbrOfElements())
        {
            // reallocate data array
            m_data.reset(new T[other.getNbrOfElements()]);
        }

        m_rows          = other.rows();
        m_cols          = other.cols();
        m_nbrOfElements = m_rows * m_cols;

        // copy data
        copyMatData(other, *this);
    }

    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& other)
{
    if (this != &other) // self-assignment check expected
    {
        this->m_data = std::move( other.m_data );

        m_rows          = other.rows();
        m_cols          = other.cols();
        m_nbrOfElements = m_rows * m_cols;
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
    for (size_t i  = 0; i < m_nbrOfElements; i++)
        dataPtr[i] = val;
}

template <class T>
inline void Matrix<T>::setValue(size_t m, size_t n, T val)
{
    data()[m * cols() + n] = val;
}

template <class T>
inline T Matrix<T>::getValue(size_t m, size_t n) const
{
    return data()[m * cols() + n];
}

template <class T>
Matrix<T> Matrix<T>::row(size_t m) const
{
    if (m >= rows())
    {
        std::cout << "Row access not possible. " << m << ", max " << rows();
        std::exit(-1);
    }

    Matrix<T> ret(1, cols());
    size_t    offset = m * cols();
    std::copy(data() + offset, data() + offset + cols(), ret.data());

    return ret;
}

template <class T>
Matrix<T> Matrix<T>::column(size_t n) const
{
    if (n >= cols())
    {
        std::cout << "Column access not possible. " << n << ", max " << cols();
        std::exit(-1);
    }

    Matrix<T> ret(rows(), 1);
    for (size_t i = 0; i < rows(); i++)
        ret(i, 0) = this->getValue(i, n);

    return ret;
}

template <class T>
const T Matrix<T>::operator()(size_t m, size_t n) const
{
    return data()[m * cols() + n];
}

template <class T>
T& Matrix<T>::operator()(size_t m, size_t n)
{
    return data()[m * cols() + n];
}

template <class T>
Matrix<T> Matrix<T>::random(size_t m, size_t n, T lower, T upper)
{
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_real_distribution<> dis(lower, upper);

    Matrix<T> rand = Matrix<T>(m, n);
    for (size_t i = 0; i < rand.getNbrOfElements(); i++)
    {
        rand.data()[i] = static_cast<T>(dis(gen));
    }

    return rand;
}

template <>
inline Matrix<int> Matrix<int>::random(size_t m, size_t n, int lower, int upper)
{
    std::random_device              rd;
    std::mt19937                    gen(rd());
    std::uniform_int_distribution<> dis(lower, upper);

    Matrix<int> rand = Matrix<int>(m, n);
    for (size_t i = 0; i < rand.getNbrOfElements(); i++)
    {
        rand.data()[i] = dis(gen);
    }

    return rand;
}

template <class T>
Matrix<T> Matrix<T>::identity(size_t m)
{
    Matrix<T> ident = Matrix<T>(m, m);
    ident.setToIdentity();
    return ident;
}

template <class T>
void Matrix<T>::setToIdentity()
{
    if (rows() != cols())
    {
        std::cout << "Not a square matrix";
        std::exit(-1);
    }

    this->fill(0);
    for (size_t i = 0; i < rows(); i++)
    {
        this->setValue(i, i, 1);
    }
}

template <class T>
Matrix<T> Matrix<T>::transpose() const
{
    Matrix<T> ret(cols(), rows());
    for (size_t m = 0; m < rows(); m++)
    {
        for (size_t n = 0; n < cols(); n++)
        {
            ret(n, m) = getValue(m, n);
        }
    }

    return ret;
}

template <class T>
Matrix<T> Matrix<T>::diagonal() const
{
    size_t    diagLength = std::min(rows(), cols());
    Matrix<T> ret(diagLength, 1);

    for (size_t i = 0; i < diagLength; i++)
    {
        ret(i, 0) = getValue(i, i);
    }

    return ret;
}

template <class T>
T Matrix<T>::sum() const
{
    T ret = 0;

    const T* ptr   = data();
    size_t   nElem = getNbrOfElements();
    for (size_t i = 0; i < nElem; i++)
        ret += ptr[i];

    return ret;
}

template <class T>
Matrix<T> Matrix<T>::sumC() const
{
    Matrix<T> res( rows(), 1 );
    res.fill(0);

    for( size_t k = 0; k < cols(); k++ )
        res = res + column(k);

    return res;
}

template <class T>
Matrix<T> Matrix<T>::sumR() const
{
    Matrix<T> res( 1, cols() );
    res.fill(0);

    for( size_t k = 0; k < rows(); k++ )
        res = res + row(k);

    return res;
}

template <class T>
bool Matrix<T>::equalDimension(const Matrix<T> m1, const Matrix<T> m2)
{
    return m1.cols() == m2.cols() && m1.rows() == m2.rows();
}

template <class T>
Matrix<T> Matrix<T>::operator*(T scale) const
{
    Matrix<T> res(this->rows(), this->cols());
    T*        resD  = res.data();
    const T*  dataD = this->data();
    for (size_t i = 0; i < this->getNbrOfElements(); i++)
        resD[i]   = dataD[i] * scale;

    return res;
}

template <class T>
Matrix<T> operator*(T scale, const Matrix<T>& mat)
{
    return mat * scale;
}

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& mat) const
{
    return matMulR(mat);
}

template <class T>
Matrix<T> Matrix<T>::matMulR(const Matrix<T>& mat) const
{
    if (this->cols() != mat.rows())
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    Matrix<T> res(this->rows(), mat.cols());

    // Direct implementation for common square matrices
    if( this->isSquare() && Matrix<T>::equalDimension(*this,mat) && this->rows() < 5)
    {
        if (this->rows() == 1)
        {
            res(0,0) = mat(0,0) * getValue(0,0);
        }
        else if( this->rows() == 2)
        {
            res(0,0) = getValue(0,0) * mat(0,0) + getValue(0,1) * mat(1,0);
            res(0,1) = getValue(0,0) * mat(0,1) + getValue(0,1) * mat(1,1);
            res(1,0) = getValue(1,0) * mat(0,0) + getValue(1,1) * mat(1,0);
            res(1,1) = getValue(1,0) * mat(0,1) + getValue(1,1) * mat(1,1);
        }
        else if( this->rows() == 3)
        {
            res(0,0) = getValue(0,0) * mat(0,0) + getValue(0,1) * mat(1,0) + getValue(0,2) * mat(2,0);
            res(0,1) = getValue(0,0) * mat(0,1) + getValue(0,1) * mat(1,1) + getValue(0,2) * mat(2,1);
            res(0,2) = getValue(0,0) * mat(0,2) + getValue(0,1) * mat(1,2) + getValue(0,2) * mat(2,2);

            res(1,0) = getValue(1,0) * mat(0,0) + getValue(1,1) * mat(1,0) + getValue(1,2) * mat(2,0);
            res(1,1) = getValue(1,0) * mat(0,1) + getValue(1,1) * mat(1,1) + getValue(1,2) * mat(2,1);
            res(1,2) = getValue(1,0) * mat(0,2) + getValue(1,1) * mat(1,2) + getValue(1,2) * mat(2,2);

            res(2,0) = getValue(2,0) * mat(0,0) + getValue(2,1) * mat(1,0) + getValue(2,2) * mat(2,0);
            res(2,1) = getValue(2,0) * mat(0,1) + getValue(2,1) * mat(1,1) + getValue(2,2) * mat(2,1);
            res(2,2) = getValue(2,0) * mat(0,2) + getValue(2,1) * mat(1,2) + getValue(2,2) * mat(2,2);
        }
        else if( this->rows() == 4)
        {
            res(0,0) = getValue(0,0)*mat(0,0) + getValue(0,1)*mat(1,0) +  getValue(0,2)*mat(2,0) + getValue(0,3)*mat(3,0);
            res(0,1) = getValue(0,0)*mat(0,1) + getValue(0,1)*mat(1,1) +  getValue(0,2)*mat(2,1) + getValue(0,3)*mat(3,1);
            res(0,2) = getValue(0,0)*mat(0,2) + getValue(0,1)*mat(1,2) +  getValue(0,2)*mat(2,2) + getValue(0,3)*mat(3,2);
            res(0,3) = getValue(0,0)*mat(0,3) + getValue(0,1)*mat(1,3) +  getValue(0,2)*mat(2,3) + getValue(0,3)*mat(3,3);

            res(1,0) = getValue(1,0)*mat(0,0) + getValue(1,1)*mat(1,0) +  getValue(1,2)*mat(2,0) + getValue(1,3)*mat(3,0);
            res(1,1) = getValue(1,0)*mat(0,1) + getValue(1,1)*mat(1,1) +  getValue(1,2)*mat(2,1) + getValue(1,3)*mat(3,1);
            res(1,2) = getValue(1,0)*mat(0,2) + getValue(1,1)*mat(1,2) +  getValue(1,2)*mat(2,2) + getValue(1,3)*mat(3,2);
            res(1,3) = getValue(1,0)*mat(0,3) + getValue(1,1)*mat(1,3) +  getValue(1,2)*mat(2,3) + getValue(1,3)*mat(3,3);

            res(2,0) = getValue(2,0)*mat(0,0) + getValue(2,1)*mat(1,0) +  getValue(2,2)*mat(2,0) + getValue(2,3)*mat(3,0);
            res(2,1) = getValue(2,0)*mat(0,1) + getValue(2,1)*mat(1,1) +  getValue(2,2)*mat(2,1) + getValue(2,3)*mat(3,1);
            res(2,2) = getValue(2,0)*mat(0,2) + getValue(2,1)*mat(1,2) +  getValue(2,2)*mat(2,2) + getValue(2,3)*mat(3,2);
            res(2,3) = getValue(2,0)*mat(0,3) + getValue(2,1)*mat(1,3) +  getValue(2,2)*mat(2,3) + getValue(2,3)*mat(3,3);

            res(3,0) = getValue(3,0)*mat(0,0) + getValue(3,1)*mat(1,0) +  getValue(3,2)*mat(2,0) + getValue(3,3)*mat(3,0);
            res(3,1) = getValue(3,0)*mat(0,1) + getValue(3,1)*mat(1,1) +  getValue(3,2)*mat(2,1) + getValue(3,3)*mat(3,1);
            res(3,2) = getValue(3,0)*mat(0,2) + getValue(3,1)*mat(1,2) +  getValue(3,2)*mat(2,2) + getValue(3,3)*mat(3,2);
            res(3,3) = getValue(3,0)*mat(0,3) + getValue(3,1)*mat(1,3) +  getValue(3,2)*mat(2,3) + getValue(3,3)*mat(3,3);
        }
    }
    else
    {
        //Create a lookup table of the right matrix split to column vectors
        // -> this is helpful because this way the column vector is in a
        //    continous memory block.

        std::vector<Matrix<T>> lookup;
        lookup.reserve(mat.cols());
        for (size_t k = 0; k < mat.cols(); k++)
            lookup.emplace_back(mat.column(k)); //todo: the intermediate matrix column() is copied again

        for (size_t n = 0; n < mat.cols(); n++)
        {
            for (size_t m = 0; m < this->rows(); m++)
            {
                // get column mat -> continous memory block
                const T *matColNPtr = lookup.at(n).data();
                const T *thisRowPtr = data() + m * cols(); // this is fast

                res(m, n) = elementwiseMultiplyAndSum(thisRowPtr, matColNPtr, this->cols());
            }
        }
    }

    return std::move(res);
}

template <class T>
Matrix<T> Matrix<T>::operator/(const Matrix<T>& mat) const
{
    // check dimension
    if (!equalDimension(*this, mat))
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    Matrix<T> res      = Matrix<T>(rows(), cols());
    T*        resD     = res.data();
    const T*  thisData = data();
    const T*  matD     = mat.data();
    for (size_t i = 0; i < getNbrOfElements(); i++)
        resD[i]   = thisData[i] / matD[i];

    return res;
}

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& mat) const
{
    // check dimension
    if (!equalDimension(*this, mat))
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    Matrix<T> res      = Matrix<T>(rows(), cols());
    T*        resD     = res.data();
    const T*  thisData = data();
    const T*  matD     = mat.data();
    for (size_t i = 0; i < getNbrOfElements(); i++)
        resD[i]   = thisData[i] + matD[i];

    return res;
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& mat) const
{
    // check dimension
    if (!equalDimension(*this, mat))
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    Matrix<T> res      = Matrix<T>(rows(), cols());
    T*        resD     = res.data();
    const T*  thisData = data();
    const T*  matD     = mat.data();
    for (size_t i = 0; i < getNbrOfElements(); i++)
        resD[i]   = thisData[i] - matD[i];

    return res;
}

template <class T>
std::ostream& operator<<(std::ostream& os, Matrix<T> const& mat)
{
    for (size_t m = 0; m < mat.rows(); m++)
    {
        for (size_t n = 0; n < mat.cols(); n++)
        {
            os << mat(m, n);
            if (n + 1 < mat.cols())
                os << ", ";
        }

        os << std::endl;
    }

    return os;
}

// from http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
inline bool almost_equal(double x, double y)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x - y) <= std::numeric_limits<double>::epsilon() * std::abs(x + y) * 2 || std::abs(x - y) < std::numeric_limits<double>::min();
}

template <class T>
bool Matrix<T>::compare(const Matrix<T>& mat, bool useCustomTolerance, T customTolerance ) const
{
    if (!equalDimension(*this, mat))
    {
        std::cout << "mismatching matrix size";
        return false;
    }

    const T* thisData = data();
    const T* matD     = mat.data();
    for (size_t i = 0; i < getNbrOfElements(); i++)
    {
        T x = thisData[i];
        T y = matD[i];

        T absDiff = std::abs(x - y);

        if( useCustomTolerance )
        {
            if( absDiff > customTolerance )
                return false;
        }
        else
        {
            // automatic tolerance
            if (absDiff > std::numeric_limits<double>::epsilon() * std::abs(x + y) * 2)
                return false;
        }

    }

    return true;
}

template <class T>
T Matrix<T>::elementwiseMultiplyAndSum(const T* arr1, const T* arr2, size_t length) const
{
    T accum = 0;
    for (size_t a = 0; a < length; a++)
        accum += arr1[a] * arr2[a];

    return accum;
}

template <>
inline int Matrix<int>::elementwiseMultiplyAndSum(const int* arr1, const int* arr2, size_t length) const
{
    __m128i* arr_1_ptr = (__m128i*)arr1;
    __m128i* arr_2_ptr = (__m128i*)arr2;

    int    accum       = 0;
    size_t simdLoops   = length / 4;
    size_t singleLoops = length - 4 * simdLoops;

    for (size_t simd = 0; simd < simdLoops; simd++)
    {
        __m128i reg_1_SSE = _mm_load_si128(arr_1_ptr);
        __m128i reg_2_SSE = _mm_load_si128(arr_2_ptr);
        __m128i mulRes    = _mm_mullo_epi32(reg_1_SSE, reg_2_SSE);
        __m128i sumRes    = _mm_hadd_epi32(mulRes, mulRes);
        sumRes            = _mm_hadd_epi32(sumRes, sumRes);

        accum += _mm_extract_epi32(sumRes, 0);

        arr_1_ptr++;
        arr_2_ptr++;
    }

    const int* offset1 = arr1 + 4 * simdLoops;
    const int* offset2 = arr2 + 4 * simdLoops;
    for (size_t a = 0; a < singleLoops; a++)
    {
        accum += offset1[a] * offset2[a];
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
    T        maxVal  = std::numeric_limits<T>::min();
    size_t   maxPos  = 0;
    const T* matData = data();

    for (size_t i = 0; i < getNbrOfElements(); i++)
    {
        if (maxVal < matData[i])
        {
            maxPos = i;
            maxVal = matData[i];
        }
    }

    size_t m = maxPos / cols();
    size_t n = maxPos - (m * cols());

    return std::make_tuple(m, n, maxVal);
}

template <class T>
std::tuple<size_t, size_t, T> Matrix<T>::min() const
{
    T        minVal  = std::numeric_limits<T>::max();
    size_t   minPos  = 0;
    const T* matData = data();

    for (size_t i = 0; i < getNbrOfElements(); i++)
    {
        if (minVal > matData[i])
        {
            minPos = i;
            minVal = matData[i];
        }
    }

    size_t m = minPos / cols();
    size_t n = minPos - (m * cols());

    return std::make_tuple(m, n, minVal);
}

template <class T>
double Matrix<T>::mean() const
{
    return static_cast<double>(sum()) / static_cast<double>(getNbrOfElements());
}

template <class T>
T* Matrix<T>::getRowPtr(size_t row)
{
    return data() + row * cols();
}

template <class T>
const T* Matrix<T>::getRowPtr(size_t row) const
{
    return data() + row * cols();
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
        std::swap(row1[i], row2[i]);
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

    if (row.cols() != cols())
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    T*       dstPtr = getRowPtr(rowIdx);
    const T* srcPtr = row.data();
    for (size_t i = 0; i < row.cols(); i++)
        dstPtr[i] = srcPtr[i];
}

template <class T>
void Matrix<T>::setColumn(size_t colIdx, const Matrix<T>& col)
{
    if (colIdx >= cols())
    {
        std::cout << "col index exceeds matrix size";
        std::exit(-1);
    }

    if (col.rows() != rows())
    {
        std::cout << "mismatching matrix size";
        std::exit(-1);
    }

    for (size_t i = 0; i < col.rows(); i++)
        setValue(i, colIdx, col(i, 0));
}

// Predefined Matrix Types
typedef std::shared_ptr<Matrix<int>>    MatrixISP;
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
    for (size_t i = 0; i < rows(); i++)
    {
        Matrix<double> cRow = echelonMat.row(i);
        if (cRow.compare(zeroRow))
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
    if (m_rows != m_cols)
    {
        std::cout << "Needs to be square matrix" << std::endl;
        *invertable = false;
        return *this;
    }

    bool   detOk;
    double det = determinant(&detOk);

    if (detOk && std::abs(det) < std::numeric_limits<double>::min())
    {
        std::cout << "Inverse failed due to determinant equal zero:"

        << std::endl << *this << std::endl;
        *invertable = false;
        return *this;
    }

    *invertable = true;

    std::vector<Matrix<double>> ops;
    Transformation::reduced_echelon(*this, ops);

    Matrix<double> inv = Matrix<double>(m_rows, m_cols);
    inv.setToIdentity();

    // multiply Echelon operators in reverse order
    // http://stattrek.com/matrix-algebra/how-to-find-inverse.aspx
    while (!ops.empty())
    {
        inv = inv * ops.back();
        ops.pop_back();
    }

    return inv;
}

#include "decomposition.hpp"

template <class T>
double Matrix<T>::determinant(bool* successful) const
{
    // compute determinant by using LU decomposition
    // infos: https://s-mat-pcs.oulu.fi/~mpa/matreng/eem3_4-3.htm
    //        -> the example is wrong!
    // https://math.stackexchange.com/questions/831823/finding-determinant-of-44-matrix-via-lu-decomposition

    double det = 0.0;

    // needs to be square matrix
    if (m_rows != m_cols)
    {
        std::cout << "Needs to be square matrix" << std::endl;
        *successful = false;
        return det;
    }

    if( m_rows == 2 ) // m_cols eq. m_rows
    {
        // special case 2x2 matrix
        det = getValue(0,0) * getValue(1,1) - getValue(0,1) * getValue(1,0);
    }
    else
    {
        // sum of product of l and u diag elements.
        // by definition, product of l diag is 1!
        Decomposition::LUResult lu = Decomposition::luDecomposition(*this);
        Matrix<double> uDiag = lu.U.diagonal();

        // build diagonal products to get determinant of u
        double detU = 1.0;

        for (size_t i = 0; i < uDiag.rows(); i++)
        {
            detU = detU * uDiag(i, 0);
        }

        // determinant of permutation matrix: 1.0 or -1.0
        double detP;
        if( lu.NbrRowSwaps % 2 == 0 ) // check if even
            detP = 1.0;
        else
            detP = -1.0;

        det =  detP*detU;
    }

    *successful = true;
    return det;
}

template <class T>
Matrix<double> Matrix<T>::normalizeColumns() const
{
    Matrix<double> retMat = Matrix<double>(m_rows, m_cols);

    for (size_t n = 0; n < m_cols; n++)
    {
        Matrix<double> c = column(n);

        // root(square and add) -> length
        double length = std::sqrt(static_cast<double>((c.transpose() * c)(0, 0)));

        // scale column and copy it to the return matrix
        retMat.setColumn(n, c * (1.0 / length));
    }

    return retMat;
}

template <class T>
Matrix<double> Matrix<T>::normalize(double& mean, double& scale) const
{
    double minP = std::get<2>(this->min());
    double maxP = std::get<2>(this->max());
    mean = this->mean();

    double diff = maxP - minP;


    if( diff < std::numeric_limits<double>::min())
    {
        // all value are identical -> do not apply scale
        scale = 1.0;
    }
    else
    {
        scale = 2.0 / diff;
    }


    Matrix<double> ret(this->rows(), this->cols());

    const T* src = this->data();
    double* dst = ret.data();

    for( size_t i = 0; i < getNbrOfElements(); i++ )
    {
        dst[i] = (src[i] - mean) * scale;
    }

    return ret;
}

template <class T>
Matrix<double> Matrix<T>::denormalize(double mean, double scale) const
{
    Matrix<double> ret(this->rows(), this->cols());

    const T *src = this->data();
    double *dst = ret.data();

    for (size_t i = 0; i < getNbrOfElements(); i++)
    {
        dst[i] = src[i] / scale + mean;
    }

    return ret;
}

template <class T>
Matrix<double> Matrix<T>::firstMinors() const
{
    if (m_rows != m_cols)
    {
        std::cout << "First minors: Square matrix required" << std::endl;
        std::exit(-1);
    }

    // https://en.wikipedia.org/wiki/Minor_(linear_algebra)

    Matrix<double> fm = Matrix<double>(rows(),cols());
    Matrix<double> tDouble(*this);

    for( size_t m = 0; m < rows(); m++ )
    {
        for( size_t n = 0; n < cols(); n++ )
        {
            Matrix<double> subMat(tDouble);
            subMat.removeRow(m);
            subMat.removeColumn(n);

            bool detCalc;
            double det = subMat.determinant(&detCalc);

            if( detCalc )
            {
                fm(m,n) = det;
            }
            else
            {
                std::cout << "Compute first minors failed" << std::endl;
                std::exit(-1);
            }
        }
    }

    return fm;
}

template <class T>
Matrix<double> Matrix<T>::cofactors() const
{
    if (m_rows != m_cols)
    {
        std::cout << "Cofactors: Square matrix required" << std::endl;
        std::exit(-1);
    }

    Matrix<double> fm = firstMinors();

    // modify signs
    for(size_t m = 0; m < m_rows; m++)
    {
        for(size_t n = 0; n < m_cols; n++)
        {
            if( (m+n) % 2 != 0 ) // check if uneven
            {
                fm(m,n) = -fm(m,n);
            }
        }
    }

    return fm;
}

template <class T>
Matrix<double> Matrix<T>::adjugate() const
{
    if (m_rows != m_cols)
    {
        std::cout << "Adjugate: Square matrix required" << std::endl;
        std::exit(-1);
    }

    return cofactors().transpose();
}

template <class T>
void Matrix<T>::removeRow(size_t m)
{
    if (m >= rows())
    {
        std::cout << "Row removal not possible. " << m << ", max " << rows();
        std::exit(-1);
    }

    // swap the row to remove step by step to the very end
    // (swapping a row is relatively fast)
    for(size_t i = m; (i+1) < rows(); i++)
    {
        swapRows(i,i+1);
    }

    // adapt matrix shape - but not allocated memory
    m_rows = rows()-1;
    m_nbrOfElements = m_cols*m_rows;
}

template <class T>
void Matrix<T>::removeColumn(size_t n)
{
    if (n >= cols())
    {
        std::cout << "Column removal not possible. " << n << ", max " << cols();
        std::exit(-1);
    }

    size_t new_memLoc = 0;

    // Different to rows, columns are not mapped
    // to the memory. Therefore modifying a column
    // is relatively time consuming.
    Matrix<T> cpy = Matrix<T>( *this );
    for( size_t m = 0; m < rows(); m++ )
    {
        for( size_t c = 0; c < cols(); c++ )
        {
            if( c != n )
            {
                data()[new_memLoc++] = cpy(m,c);
            }
        }
    }

    // adapt matrix shape - but not allocated memory
    m_cols = cols()-1;
    m_nbrOfElements = m_cols*m_rows;
}

template <class T>
Matrix<T> Matrix<T>::subMatrix(size_t rowStart, size_t columnStart, size_t nbrOfRows, size_t nbrOfColumns) const
{
    if( rowStart + nbrOfRows > rows() || columnStart + nbrOfColumns > cols() )
    {
        std::cout << "subMatrix exceeds actual matrix size";
        std::exit(-1);
    }

    Matrix<T> res( nbrOfRows, nbrOfColumns );
    T* dst = res.data();
    const T* src = data() + rowStart*cols() + columnStart;

    for(size_t m = 0; m < nbrOfRows; m++ )
    {
        // copy partial lines at once
        const T* c_Src = src + m*cols();
        T* c_Dst = dst + m*nbrOfColumns;

        std::copy( c_Src, c_Src + nbrOfColumns, c_Dst);
    }

    return res;
}

template <class T>
template <class R>
void Matrix<T>::setSubMatrix(size_t rowStart, size_t columnStart, const Matrix<R>& subMat)
{
    if( rowStart + subMat.rows() > rows() || columnStart + subMat.cols() > cols() )
    {
        std::cout << "writing subMatrix exceeds actual matrix size";
        std::exit(-1);
    }

    for( size_t m = 0; m < subMat.rows(); m++ )
    {
        for( size_t n = 0; n < subMat.cols(); n++ )
        {
            setValue(m+rowStart, n+columnStart, subMat(m,n));
        }
    }
}

template <class T>
double Matrix<T>::norm() const
{
    return std::sqrt( normSquare() );
}

template <class T>
T Matrix<T>::normSquare() const
{
    T normSq = 0;

    if (cols() == 1 || rows() == 1)
    {
        const T *matData = data();
        for (size_t i = 0; i < getNbrOfElements(); i++)
        {
            T val = matData[i];
            normSq += val * val;
        }
    }
    else
    {
        std::cout << "NormSquare: needs to be a row or column vector." << std::endl;
    }

    return normSq;
}

template <class T>
bool Matrix<T>::isSymmetric() const
{
    return compare( transpose() );
}

template <class T>
bool Matrix<T>::isSquare() const
{
    return rows() == cols();
}

template <class T>
bool Matrix<T>::isOrthogonal(T customTolerance) const
{
    if( !isSquare() )
    {
        return false;
    }

    for( size_t i = 0; i < rows(); i++)
    {
        if(  std::abs(1.0 - row(i).norm() )  > customTolerance )
        {
            return false;
        }

        if(  std::abs(1.0 - column(i).norm())  > customTolerance )
        {
            return false;
        }
    }

    return true;
}

template <class T>
bool Matrix<T>::save(const std::string& path) const
{
    std::ofstream f;
    f.open(path, std::ofstream::trunc);
    if( ! f.is_open() )
    {
        std::cout << "Cannot save file " << path << std::endl;
        return false;
    }

    f << serialize();
    f.close();

    return true;
}

template <class T>
bool Matrix<T>::load(const std::string& path)
{
    std::ifstream f( path );
    if( ! f.is_open() )
    {
        std::cout << "Cannot open file " << path << std::endl;
        return false;
    }

    // read the whole file
    std::string buffer((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    f.close();

    deserialize(buffer);
    return true;
}

template <class T>
std::string Matrix<T>::serialize() const
{
    std::string data;

    size_t* dim = new size_t[2];
    dim[0] = m_rows;
    dim[1] = m_cols;

    data.append( std::string( (char*)dim, 2*sizeof(size_t) ) );
    data.append( std::string( (char*)this->data(), m_nbrOfElements* sizeof(T) ));

    delete[] dim;
    return data;
}

template <class T>
void Matrix<T>::deserialize(const std::string& data)
{
    const char* buffer = data.c_str();

    // init members
    m_rows = ((size_t*)(buffer))[0];
    m_cols = ((size_t*)(buffer))[1];

    m_nbrOfElements = m_rows * m_cols;
    m_data.reset(new T[m_nbrOfElements]);

    // copy data
    const T* src = (const T*)(buffer + 2* sizeof(size_t));
    std::copy(src, src + m_nbrOfElements, this->data());
}

template <class T>
void Matrix<T>::setToZero( T threshold )
{
    std::transform(data(), data()+m_nbrOfElements, data(),
                   [&threshold](T val) -> T
                   {
                       if( std::abs(val) < threshold)
                           return 0;
                       else
                           return val;
                   });
}

template <class T>
T Matrix<T>::normL1() const
{
    T norm = 0;
    Matrix<T> thisAbs = this->absolute();

    for( size_t c = 0; c < thisAbs.cols(); c++ )
    {
        T cSum = thisAbs.column(c).sum();
        if( cSum > norm )
            norm = cSum;
    }

    return norm;
}

template <class T>
T Matrix<T>::normInf() const
{
    T norm = 0;
    Matrix<T> thisAbs = this->absolute();

    for( size_t c = 0; c < thisAbs.rows(); c++ )
    {
        T cSum = thisAbs.row(c).sum();
        if( cSum > norm )
            norm = cSum;
    }

    return norm;
}

template <class T>
double Matrix<T>::normL2() const
{
    double normRet = 0;

    if (cols() == 1 || rows() == 1)
    {
        // if vector
        normRet = norm();
    }
    else
    {
        // if matrix -> max singular value
        Decomposition::SVDResult svdR = Decomposition::svd( *this );

        auto max = svdR.S.max();
        normRet = std::get<2>(max);
    }

    return normRet;
}

template <class T>
Matrix<T> Matrix<T>::absolute() const
{
    Matrix<T> ret( this->rows(), this->cols() );
    const T* src = this->data();
    T* dst = ret.data();
    for( size_t k = 0; k < ret.getNbrOfElements(); k++ )
    {
        dst[k] = std::abs(src[k]);
    }

    return ret;
}

template <class T>
Matrix<T> Matrix<T>::repMat(size_t rowDirection, size_t columnDirection) const
{
    if( rowDirection < 1 || columnDirection < 1 )
        throw InvalidInputException();

    Matrix<T> ret = Matrix<T>( rows() * rowDirection, cols() * columnDirection );
    for( size_t x = 0; x < columnDirection; x++ )
        for( size_t y = 0; y < rowDirection; y++ )
            ret.setSubMatrix(y * rows(), x * cols(), *this);

    return ret;
}

template <class T>
double Matrix<T>::conditionNumberL1() const
{
    bool invertable;
    Matrix<T> invMat = this->inverted(&invertable);

    if( !invertable )
    {
        std::cerr << "conditionNumberL1: matrix not invertable";
        return 0.0;
    }

    double n = normL1();
    double ni = invMat.normL1();
    return n * ni;
}

/**
 * Computes the infinity norm condition number.
 * @return Infinity condition number
 */
template <class T>
double Matrix<T>::conditionNumberInf() const
{
    bool invertable;
    Matrix<T> invMat = this->inverted(&invertable);

    if( !invertable )
    {
        std::cerr << "conditionNumberInf: matrix not invertable";
        return 0.0;
    }

    double n = normInf();
    double ni = invMat.normInf();
    return n * ni;
}

#ifdef OPENCVEIDLA


// general case -> see template specializations below

template <class T>
cv::Mat Matrix<T>::toOpenCVMatrix( ) const
{
    return createOpenCVMat();
}

// general case, return double matrix -> see template specializations below
template <class T>
inline cv::Mat Matrix<T>::createOpenCVMat( ) const
{
    cv::Mat oMat(rows(),cols(),CV_64F);

    for(size_t m = 0; m < rows(); m++)
    {
        for(size_t n = 0; n < cols(); n++)
        {
            oMat.at<double>(m,n) = (*this)(m,n);
        }
    }

    return oMat;
}

template <>
inline cv::Mat Matrix<double>::createOpenCVMat( ) const
{
    cv::Mat oMat(rows(),cols(),CV_64F);

    for(size_t m = 0; m < rows(); m++)
    {
        double* dst = oMat.ptr<double>(m);
        std::copy(getRowPtr(m), getRowPtr(m) + cols(), dst);
    }

    return oMat;
}

template <>
inline cv::Mat Matrix<float>::createOpenCVMat( ) const
{
    cv::Mat oMat(rows(),cols(),CV_32F);

    for(size_t m = 0; m < rows(); m++)
    {
        float* dst = oMat.ptr<float>(m);
        std::copy(getRowPtr(m), getRowPtr(m) + cols(), dst);
    }

    return oMat;
}


template <>
inline cv::Mat Matrix<int>::createOpenCVMat( ) const
{
    cv::Mat oMat(rows(),cols(),CV_32S);

    for(size_t m = 0; m < rows(); m++)
    {
        int* dst = oMat.ptr<int>(m);
        std::copy(getRowPtr(m), getRowPtr(m) + cols(), dst);
    }

    return oMat;
}

template <>
inline cv::Mat Matrix<uchar>::createOpenCVMat( ) const
{
    cv::Mat oMat(rows(),cols(),CV_8UC1);

    for(size_t m = 0; m < rows(); m++)
    {
        uchar* dst = oMat.ptr<uchar>(m);
        std::copy(getRowPtr(m), getRowPtr(m) + cols(), dst);
    }

    return oMat;
}


#endif // OPENCVEIDLA

#endif //MY_MATRIX_H