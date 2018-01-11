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

#include "matrix4x4.hpp"

Matrix4x4::Matrix4x4()
: Matrix<double>(4, 4)
{
    setToIdentity();
}

Matrix4x4::Matrix4x4(   double m00, double m01, double m02, double m03,
                        double m10, double m11, double m12, double m13,
                        double m20, double m21, double m22, double m23,
                        double m30, double m31, double m32, double m33 )
: Matrix<double>(4,4)
{
    setValue(0,0, m00); setValue(0,1, m01); setValue(0,2, m02); setValue(0,3, m03);
    setValue(1,0, m10); setValue(1,1, m11); setValue(1,2, m12); setValue(1,3, m13);
    setValue(2,0, m20); setValue(2,1, m21); setValue(2,2, m22); setValue(2,3, m23);
    setValue(3,0, m30); setValue(3,1, m31); setValue(3,2, m32); setValue(3,3, m33);
}


Matrix4x4& Matrix4x4::operator=(const Matrix4x4& other)
{
    if (this != &other) // self-assignment check expected
    {
        // copy data
        copyMatData(other, *this);
    }

    return *this;
}

Matrix4x4& Matrix4x4::operator=(const Matrix<double>& other)
{
    if (this != &other) // self-assignment check expected
    {
        if( Matrix::equalDimension(*this, other) )
        {
            copyMatData(other, *this);
        }
        else
        {
            std::cout << "Cannot assign Matrix4x4 from matrix with unequal dimension" << std::endl;
            std::exit(-1);
        }
    }

    return *this;
}


Matrix4x4::Matrix4x4(const Matrix4x4& mat)
        : Matrix4x4()
{
    copyMatData(mat, *this);
}


// https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations

void Matrix4x4::rotZ(double radian)
{
    Matrix4x4 rotZMat( std::cos(radian), -std::sin(radian), 0.0, 0.0,
                       std::sin(radian), std::cos(radian), 0.0, 0.0,
                       0.0, 0.0, 1.0, 0.0,
                       0.0, 0.0, 0.0, 1.0);

    *this = rotZMat * (*this);
}

void Matrix4x4::rotY(double radian)
{
    Matrix4x4 rotYMat( std::cos(radian), 0.0, std::sin(radian), 0.0,
                       0.0, 1.0, 0.0, 0.0,
                       -std::sin(radian), 0.0, std::cos(radian), 0.0,
                       0.0, 0.0, 0.0, 1.0 );

    *this = rotYMat * (*this);
}

void Matrix4x4::rotX(double radian)
{
    Matrix4x4 rotXMat( 1.0, 0.0, 0.0, 0.0,
                       0.0, std::cos(radian), -std::sin(radian), 0.0,
                       0.0, std::sin(radian),  std::cos(radian), 0.0,
                       0.0, 0.0, 0.0, 1.0 );

    *this = rotXMat * (*this);
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4& m) const
{
    Matrix4x4 res;
    const Matrix4x4& t = *this;

    res(0,0) = t(0,0)*m(0,0) + t(0,1)*m(1,0) +  t(0,2)*m(2,0) + t(0,3)*m(3,0);
    res(0,1) = t(0,0)*m(0,1) + t(0,1)*m(1,1) +  t(0,2)*m(2,1) + t(0,3)*m(3,1);
    res(0,2) = t(0,0)*m(0,2) + t(0,1)*m(1,2) +  t(0,2)*m(2,2) + t(0,3)*m(3,2);
    res(0,3) = t(0,0)*m(0,3) + t(0,1)*m(1,3) +  t(0,2)*m(2,3) + t(0,3)*m(3,3);

    res(1,0) = t(1,0)*m(0,0) + t(1,1)*m(1,0) +  t(1,2)*m(2,0) + t(1,3)*m(3,0);
    res(1,1) = t(1,0)*m(0,1) + t(1,1)*m(1,1) +  t(1,2)*m(2,1) + t(1,3)*m(3,1);
    res(1,2) = t(1,0)*m(0,2) + t(1,1)*m(1,2) +  t(1,2)*m(2,2) + t(1,3)*m(3,2);
    res(1,3) = t(1,0)*m(0,3) + t(1,1)*m(1,3) +  t(1,2)*m(2,3) + t(1,3)*m(3,3);

    res(2,0) = t(2,0)*m(0,0) + t(2,1)*m(1,0) +  t(2,2)*m(2,0) + t(2,3)*m(3,0);
    res(2,1) = t(2,0)*m(0,1) + t(2,1)*m(1,1) +  t(2,2)*m(2,1) + t(2,3)*m(3,1);
    res(2,2) = t(2,0)*m(0,2) + t(2,1)*m(1,2) +  t(2,2)*m(2,2) + t(2,3)*m(3,2);
    res(2,3) = t(2,0)*m(0,3) + t(2,1)*m(1,3) +  t(2,2)*m(2,3) + t(2,3)*m(3,3);

    res(3,0) = t(3,0)*m(0,0) + t(3,1)*m(1,0) +  t(3,2)*m(2,0) + t(3,3)*m(3,0);
    res(3,1) = t(3,0)*m(0,1) + t(3,1)*m(1,1) +  t(3,2)*m(2,1) + t(3,3)*m(3,1);
    res(3,2) = t(3,0)*m(0,2) + t(3,1)*m(1,2) +  t(3,2)*m(2,2) + t(3,3)*m(3,2);
    res(3,3) = t(3,0)*m(0,3) + t(3,1)*m(1,3) +  t(3,2)*m(2,3) + t(3,3)*m(3,3);

    return res;
}

