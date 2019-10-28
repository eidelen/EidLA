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
#include "exceptions.hpp"

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
    Matrix4x4()
    : Matrix<double>(4, 4)
    {
        setToIdentity();
    }

    /**
     * Constructs a 4x4 double matrix with content.
     * @param m(M,N) Matrix entry m of row M and column N.
     */
    Matrix4x4(double m00, double m01, double m02, double m03,
              double m10, double m11, double m12, double m13,
              double m20, double m21, double m22, double m23,
              double m30, double m31, double m32, double m33)
    : Matrix<double>(4, 4)
    {
        setValue(0, 0, m00);
        setValue(0, 1, m01);
        setValue(0, 2, m02);
        setValue(0, 3, m03);
        setValue(1, 0, m10);
        setValue(1, 1, m11);
        setValue(1, 2, m12);
        setValue(1, 3, m13);
        setValue(2, 0, m20);
        setValue(2, 1, m21);
        setValue(2, 2, m22);
        setValue(2, 3, m23);
        setValue(3, 0, m30);
        setValue(3, 1, m31);
        setValue(3, 2, m32);
        setValue(3, 3, m33);
    }

    /**
     * Copy constructor.
     * @param mat Matrix from which to copy
     */
    template <class R>
    Matrix4x4(const Matrix<R>& mat);

    Matrix4x4(const Matrix4x4& mat)
    : Matrix<double>(4, 4)
    {
        copyMatData(mat, *this);
    }

    virtual ~Matrix4x4()
    {
    }

    /**
     * Asignment operator: overwrite content of this matrix.
     * @param other
     * @return
     */
    Matrix4x4& operator=(const Matrix4x4& other)
    {
        if (this != &other) // self-assignment check expected
        {
            // copy data
            copyMatData(other, *this);
        }

        return *this;
    }

    /**
     * Asignment operator: overwrite content of this matrix.
     * @param other
     * @return
     */
    Matrix4x4& operator=(const Matrix<double>& other)
    {
        if (this != &other) // self-assignment check expected
        {
            if (Matrix::equalDimension(*this, other))
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

    /**
     * This function returns the upper left 3x3
     * matrix. In case of rigid transformations,
     * this is the rotation matrix.
     * @return Rotation matrix.
     */
    Matrix<double> getRotation() const
    {
        return subMatrix(0, 0, 3, 3);
    }

    /**
     * This function returns the most right column.
     * In case of rigid transformations, this is
     * the translation vector.
     * @return Translation vector.
     */
    Matrix<double> getTranslation() const
    {
        return subMatrix(0, 3, 3, 1);
    }

    /**
     * Sets the 3x3 upper left rotation matrix.
     * @param rot 3x3 rotation matrix.
     */
    void setRotation(const Matrix<double>& rot)
    {
        if (rot.cols() != 3 || rot.rows() != 3)
        {
            std::cout << "Rotation matrix has wrong dimension";
            std::exit(-1);
        }

        const double* src = rot.data();
        double*       dst = data();

        std::copy(src, src + 3, dst);
        std::copy(src + 3, src + 6, dst + 4);
        std::copy(src + 6, src + 9, dst + 8);
    }

    /**
     * Sets the most right column, which is the
     * translation vector.
     * @param trans Translation vector.
     */
    void setTranslation(const Matrix<double>& trans)
    {
        if (trans.cols() != 1 || (trans.rows() != 3 && trans.rows() != 4))
            throw NoVectorException();

        setTranslation(trans(0, 0), trans(1, 0), trans(2, 0));
    }

    /**
     * Sets the most right column, which is the
     * translation vector.
     * @param x
     * @param y
     * @param z
     */
    void setTranslation(double x, double y, double z)
    {
        setValue(0, 3, x);
        setValue(1, 3, y);
        setValue(2, 3, z);
        setValue(3, 3, 1.0);
    }

    /**
     * Rotate around the z-axis.
     * @param radian Rotation angle.
     * info: https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
     */
    void rotZ(double radian)
    {
        Matrix4x4 rotZMat(std::cos(radian), -std::sin(radian), 0.0, 0.0,
                          std::sin(radian), std::cos(radian), 0.0, 0.0,
                          0.0, 0.0, 1.0, 0.0,
                          0.0, 0.0, 0.0, 1.0);

        *this = rotZMat * (*this);
    }

    /**
     * Rotate around the y-axis.
     * @param radian Rotation angle.
     * https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
     */
    void rotY(double radian)
    {
        Matrix4x4 rotYMat(std::cos(radian), 0.0, std::sin(radian), 0.0,
                          0.0, 1.0, 0.0, 0.0,
                          -std::sin(radian), 0.0, std::cos(radian), 0.0,
                          0.0, 0.0, 0.0, 1.0);

        *this = rotYMat * (*this);
    }

    /**
     * Rotate around the x-axis.
     * @param radian Rotation angle.
     * https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
     */
    void rotX(double radian)
    {
        Matrix4x4 rotXMat(1.0, 0.0, 0.0, 0.0,
                          0.0, std::cos(radian), -std::sin(radian), 0.0,
                          0.0, std::sin(radian), std::cos(radian), 0.0,
                          0.0, 0.0, 0.0, 1.0);

        *this = rotXMat * (*this);
    }

    /**
     * Computes the inverse of this 4x4 matrix,
     * by assuming it is a rigid transformation.
     * In that way, the invert of the 3x3 rot-mat
     * can be computed fast.
     * @return Inverse of the rigid transformation.
     * note: https://math.stackexchange.com/questions/1234948/inverse-of-a-rigid-transformation
     */
    Matrix4x4 inverted_rg() const
    {
        Matrix<double> rotMat    = subMatrix(0, 0, 3, 3).transpose();
        Matrix<double> transVect = (rotMat * (-1)) * subMatrix(0, 3, 3, 1);

        Matrix4x4 ret;
        ret.setRotation(rotMat);
        ret.setTranslation(transVect);
        return ret;
    }

    /**
     * Finds the rigid transformation between two corresponding 3D
     * point sets. The algorithm is described in the paper
     * "Least-Squares Fitting of Two 3-D Point Sets".
     * @param setA Point set A (3 x n matrix)
     * @param setB Point set B (3 x n matrix)
     * @param error On return, this holds the registration error.
     * @return Rigid transformation
     */
    template <typename R, typename Q>
    static Matrix4x4 findRigidTransformation(const Matrix<R>& setA, const Matrix<Q>& setB, double& error);

    /**
     * Computes the average Euclidean (L2) error of the transformation t, which
     * maps the setA to the setB: Norm( setB - t * setA )
     * @param setA Set A
     * @param setB Set B
     * @param transformation Transformation
     * @return Average Euclidean error
     */
    template <typename R, typename Q>
    static double computeTransformationError(const Matrix<R>& setA, const Matrix<Q>& setB, const Matrix4x4& t);
};

template <class R>
Matrix4x4::Matrix4x4(const Matrix<R>& mat)
: Matrix4x4()
{
    if (Matrix::equalDimension(*this, mat))
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

template <typename R, typename Q>
Matrix4x4 Matrix4x4::findRigidTransformation(const Matrix<R>& setA, const Matrix<Q>& setB, double& error)
{
    // input to double
    Matrix<double> dA(setA);
    Matrix<double> dB(setB);

    // at least 3 point correspondences
    if (dA.cols() != dB.cols() || dA.cols() < 3 || dA.rows() != 3 || dB.rows() != 3)
        throw InvalidInputException();

    size_t n = dA.cols();

    // compute centers
    Matrix<double> cA = dA.sumC() * (1.0 / n);
    Matrix<double> cB = dB.sumC() * (1.0 / n);

    // centered point sets
    Matrix<double> qA = dA - cA.repMat(1, n);
    Matrix<double> qB = dB - cB.repMat(1, n);

    // compute rotation
    Matrix<double> h = Matrix<double>(3, 3);
    h.fill(0.0);

    for (size_t k = 0; k < n; k++)
        h = h + (qA.column(k) * qB.column(k).transpose());

    Decomposition::SVDResult dec = Decomposition::svdGolubKahan(h);

    Matrix<double> rotation = dec.V * dec.U.transpose();

    // If necessary, fix rotation matrix -> reflection
    double det = rotation.determinant();

    if (std::abs(1.0 - std::abs(det)) > 0.001)
    {
        h.save("failingsvd.mat");
        throw NoRotationMatrixException();
    }

    if (det < 0.0)
    {
        //std::cout << "Rotation is a reflection (determinant = " << det << "). Let's fix it!" << std::endl;

        Matrix<double> vModified = dec.V;
        vModified.setColumn(2, vModified.column(2) * (-1.0));

        rotation = vModified * dec.U.transpose();

        // Modified rotation determinant should be positive... actually 1
        double detFixed = rotation.determinant();
        if (detFixed < 0.0)
            throw NoRotationMatrixException();
    }

    // compose resulting transformation
    Matrix4x4 res = Matrix4x4();
    res.setRotation(rotation);
    Matrix<double> translation = cB - (rotation * cA);
    res.setTranslation(translation);

    error = computeTransformationError(dA, dB, res);

    return res;
}

template <typename R, typename Q>
double Matrix4x4::computeTransformationError(const Matrix<R>& setA, const Matrix<Q>& setB, const Matrix4x4& t)
{
    // check input
    if (!Matrix::equalDimension(setA, setB))
        throw InvalidInputException();

    if (setA.rows() < 3)
        throw InvalidInputException();

    size_t n = setA.cols();

    // works for homogene and non-homogene coordinates
    Matrix<double> setA_h(4, n);
    setA_h.fill(1.0);
    setA_h.setSubMatrix(0, 0, setA);

    // transform setA and compute difference to setB
    Matrix<double> tDiff = setB - (t * setA_h).subMatrix(0, 0, 3, n);

    // compute average L2 error
    double error = 0.0;
    for (size_t k = 0; k < n; k++)
        error += tDiff.column(k).norm();

    error = error / n;

    return error;
}

#endif //MY_AFFINE_H
