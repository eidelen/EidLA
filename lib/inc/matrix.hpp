
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

    /**
     * Constructs a m x n matrix and assigns
     * the values taken from the data array
     * in row direction.
     * @param rows number of rows
     * @param cols number of columns
     * @param data Pointer to data
     */
    Matrix(size_t rows, size_t cols, const T* data);

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
    bool compare( const Matrix<T>& mat, double epsilon = 0.0 ) const;

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

protected:
    /**
     * Check if the dimensions of the two passed matrix are equal.
     * @param m1 Mat 1
     * @param m2 Mat 2
     * @return True if equal dimension.
     */
    static bool equalDimension( const Matrix<T> m1, const Matrix<T> m2);

    T elementwiseMultiplyAndSum(const T *arr1, const T *arr2, size_t length) const;

protected:
    const size_t m_rows;
    const size_t m_cols;

    std::shared_ptr<T> m_data;
    const size_t m_nbrOfElements;
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
void copyMatData(const Matrix<T>& src, Matrix<T> dst )
{
    const T* srcPtr = src.data();
    T* dstPtr = dst.data();
    std::copy(srcPtr, srcPtr+src.getNbrOfElements(), dstPtr);
}

inline void copyMatData(const Matrix<int>& src, Matrix<double> dst )
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
Matrix<T>::Matrix(size_t rows, size_t cols, const T* data)
        : Matrix<T>(rows, cols)
{
    T* dst = this->data();
    std::copy(data, data+this->getNbrOfElements(), this->data());
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

/*
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

// Predefined Matrix Types
typedef std::shared_ptr<Matrix<int>> MatrixISP;
typedef std::shared_ptr<Matrix<double>> MatrixDSP;

#endif //MY_MATRIX_H