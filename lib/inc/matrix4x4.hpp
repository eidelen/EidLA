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

#ifndef MY_AFFINE_H
#define MY_AFFINE_H

#include "matrix.hpp"

/**
 * This class represents a 4x4 double matrix. This kind
 * of matrix is often used to define rigid transformations.
 */
class Matrix4x4 : public Matrix<double>
{
public:

    /**
     * Constructs a 4x4 identity double matrix.
     */
    Matrix4x4();

    /**
     * Constructs a 4x4 double matrix with content.
     * @param m(M,N) Matrix entry m of row M and column N.
     */
    Matrix4x4(  double m00, double m01, double m02, double m03,
                double m10, double m11, double m12, double m13,
                double m20, double m21, double m22, double m23,
                double m30, double m31, double m32, double m33 );

    /**
     * Copy constructor.
     * @param mat Matrix from which to copy
     */
    template <class R>
    Matrix4x4(const Matrix<R>& mat);

    Matrix4x4(const Matrix4x4& mat);

    /**
     * Asignment operator: overwrite content of this matrix.
     * @param other
     * @return
     */
    Matrix4x4& operator=(const Matrix4x4& other);

    /**
     * Asignment operator: overwrite content of this matrix.
     * @param other
     * @return
     */
    Matrix4x4& operator=(const Matrix<double>& other);

    /**
     * This function returns the upper left 3x3
     * matrix. In case of rigid transformations,
     * this is the rotation matrix.
     * @return Rotation matrix.
     */
    Matrix<double> getRotation() const;

    /**
     * This function returns the most right column.
     * In case of rigid transformations, this is
     * the translation vector.
     * @return Translation vector.
     */
    Matrix<double> getTranslation() const;

    /**
     * Sets the 3x3 upper left rotation matrix.
     * @param rot 3x3 rotation matrix.
     */
    void setRotation( const Matrix<double>& rot );

    /**
     * Sets the most right column, which is the
     * translation vector.
     * @param trans Translation vector.
     */
    void setTranslation( const Matrix<double>& trans);

    /**
     * Sets the most right column, which is the
     * translation vector.
     * @param x
     * @param y
     * @param z
     */
    void setTranslation( double x, double y, double z);

    /**
     * Matrix multiplication
     * @param mat
     * @return Product of two 4x4 matrix
     */
    Matrix4x4 operator*(const Matrix4x4& mat) const;
    template <class R>
    Matrix<double> operator*(const Matrix<R>& mat) const;

    /**
     * Rotate around the z-axis.
     * @param radian Rotation angle.
     */
    void rotZ(double radian);

    /**
     * Rotate around the y-axis.
     * @param radian Rotation angle.
     */
    void rotY(double radian);

    /**
     * Rotate around the x-axis.
     * @param radian Rotation angle.
     */
    void rotX(double radian);

    /**
     * Computes the inverse of this 4x4 matrix,
     * by assuming it is a rigid transformation.
     * In that way, the invert of the 3x3 rot-mat
     * can be computed fast.
     * @return Inverse of the rigid transformation.
     */
    Matrix4x4 inverted_rg() const;
};


template <class R>
Matrix4x4::Matrix4x4(const Matrix<R>& mat)
        : Matrix4x4()
{
    if( Matrix::equalDimension(*this, mat) )
    {
        // if you have a compile error here, you have to extend the
        // copyMatData function with the particular types.
        copyMatData(mat, *this);
    }
    else
    {
        std::cout << "Cannot create Matrix4x4 from matrix with unequal dimension" << std::endl;
        std::exit(-1);
    }
}

template <class R>
Matrix<double> Matrix4x4::operator*(const Matrix<R>& mat) const
{
    return matMulR(mat);
}

#endif //MY_AFFINE_H
