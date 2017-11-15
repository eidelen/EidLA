
#ifndef MY_MATRIX_H
#define MY_MATRIX_H

#include <memory>
#include <iostream>
#include <cmath>

template <class T>
class Matrix
{
public:
    Matrix( size_t rows, size_t cols );
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
    void setValue( T val );

    /**
     * Set the value at position m, n.
     * @param m Row
     * @param n Column
     * @param val Value
     */
    void setValue( size_t m, size_t n, T val);

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
    static Matrix<T> eye(size_t m);

protected:
    /**
     * Check if the dimensions of the two passed matrix are equal.
     * @param m1 Mat 1
     * @param m2 Mat 2
     * @return True if equal dimension.
     */
    static bool equalDimension( const Matrix<T> m1, const Matrix<T> m2);

protected:
    const size_t m_rows;
    const size_t m_cols;

    std::shared_ptr<T> m_data;
    const size_t m_nbrOfElements;
};


template <class T>
Matrix<T>::Matrix(size_t rows, size_t cols)
        : m_rows(rows), m_cols(cols), m_nbrOfElements(rows*cols)
{
    m_data.reset( new T[m_nbrOfElements] );
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
void Matrix<T>::setValue( T val )
{
    T* dataPtr = data();
    for( size_t i = 0; i < m_nbrOfElements; i++ )
        dataPtr[i] = val;
}

template <class T>
void Matrix<T>::setValue( size_t m, size_t n, T val)
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
    for(size_t i = 0; i < cols(); i++)
        ret.data()[i] = this->data()[offset+i];

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
Matrix<T> Matrix<T>::eye(size_t m)
{
    Matrix<T> ident = Matrix<T>(m,m);
    ident.setValue(0);
    for(size_t i = 0; i < m; i++)
    {
        ident(i,i) = 1;
    }

    return ident;
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

    for( size_t n = 0; n < mat.cols(); n++)
    {
        for( size_t m = 0; m < this->rows(); m++)
        {
            T accum = 0;
            for(size_t a = 0; a < this->cols(); a++)
            {
                accum = accum + (this->getValue(m,a) * mat(a,m));
            }
            res(m,n) = accum;
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

// Predefined Matrix Types
typedef std::shared_ptr<Matrix<int>> MatrixISP;
typedef std::shared_ptr<Matrix<double>> MatrixDSP;

#endif //MY_MATRIX_H