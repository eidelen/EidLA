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

#ifndef MY_DECOMPOSITION_H
#define MY_DECOMPOSITION_H

#include "matrix.hpp"

class Decomposition
{
public:
    struct LUResult
    {
        LUResult(Matrix<double> l, Matrix<double> u, Matrix<double> p, size_t n)
        : L(l), U(u), P(p), NbrRowSwaps(n)
        {
        }
        Matrix<double> L;           // Lower triangle matrix
        Matrix<double> U;           // Upper triangle matrix
        Matrix<double> P;           // Row Swaps
        size_t         NbrRowSwaps; // Number of row swaps
    };

    struct EigenPair
    {
        EigenPair(Matrix<double> v, double l, bool valid)
        : V(v), L(l), Valid(valid)
        {
        }
        Matrix<double> V;     // Eigen vector
        double         L;     // Eigen value
        bool           Valid; // Is eigen pair valid? It is if precision was reached.
    };

    static void sortDescending(std::vector<EigenPair>& pairs)
    {
        std::sort(pairs.begin(), pairs.end(), [](const EigenPair& a, const EigenPair& b) {
            return a.L > b.L;
        });
    }

    static void sortAscending(std::vector<EigenPair>& pairs)
    {
        std::sort(pairs.begin(), pairs.end(), [](const EigenPair& a, const EigenPair& b) {
            return a.L < b.L;
        });
    }

    static void filter(std::vector<EigenPair>& pairs, double eigenValueThreshold)
    {
        std::vector<EigenPair> copyEP = pairs;
        pairs.clear();
        std::copy_if(copyEP.begin(), copyEP.end(), std::back_inserter(pairs), [eigenValueThreshold](EigenPair i) { return i.L >= eigenValueThreshold; });
    }

    struct QRResult
    {
        QRResult(const Matrix<double>& q, const Matrix<double>& r)
        : Q(q), R(r)
        {
        }

        QRResult(Matrix<double>&& q, Matrix<double>&& r)
        {
            Q = std::move(q);
            R = std::move(r);
        }

        QRResult(const QRResult& qr)
        : Q(qr.Q), R(qr.R)
        {
        }
        QRResult()
        : Q(Matrix<double>::identity(1)), R(Matrix<double>::identity(1))
        {
        }

        QRResult& operator=(const QRResult& other)
        {
            // check for self-assignment
            if (&other == this)
                return *this;

            this->Q = other.Q;
            this->R = other.R;
            return *this;
        }

        QRResult& operator=(QRResult&& other)
        {
            // check for self-assignment
            if (&other == this)
                return *this;

            this->Q = std::move(other.Q);
            this->R = std::move(other.R);
            return *this;
        }

        Matrix<double> Q; // Orthogonal matrix
        Matrix<double> R; // Upper triangle matrix
    };

    struct SVDResult
    {
        SVDResult(Matrix<double> u, Matrix<double> s, Matrix<double> v)
        : U(u), S(s), V(v)
        {
        }
        Matrix<double> U; // Left singular vectors
        Matrix<double> S; // Diagonal matrix of singular values
        Matrix<double> V; // Right singular vectors
    };

    struct HouseholderResult
    {
        HouseholderResult(double b, Matrix<double> v)
        : B(b), V(v)
        {
        }
        double         B; // scalar
        Matrix<double> V; // householder vector
    };

    struct DiagonalizationResult
    {
        DiagonalizationResult(Matrix<double> u, Matrix<double> d, Matrix<double> v)
        : U(u), D(d), V(v)
        {
        }
        Matrix<double> U; // Left orthogonal matrix
        Matrix<double> D; // Bidiagonal matrix
        Matrix<double> V; // Right orthogonal matrix
    };

    struct GivensRotation
    {
        GivensRotation()
        : S(0), C(0), R(0)
        {
        }

        GivensRotation(double s, double c, double r)
        : S(s), C(c), R(r)
        {
        }

        double S;
        double C;
        double R;
    };

public:
    /**
     * LU decomposition of the matrix mat.
     * @param mat
     * @param pivoting Enable or disable row pivoting. If disabled, the returning P is the identity matrix.
     * @return LUResult
     */
    template <class T>
    static LUResult luDecomposition(const Matrix<T>& mat, bool pivoting = true);

    enum EigenMethod
    {
        PowerIterationAndHotellingsDeflation, //Power iteration and hotelling's deflation
        QRAlgorithm,                          // QR algorithm
    };

    /**
     * Eigen decomposition of the matrix mat. Finds all Eigen pairs
     * in symmetric real matrices. In non-symmetric matrices, it finds
     * at least the most significant Eigen pair.
     * @param mat
     * @param method Eigen computation method to used.
     * @return Vector of Eigen pairs.
     */
    template <class T>
    static std::vector<EigenPair> eigen(const Matrix<T>& mat, EigenMethod method = QRAlgorithm);

    /**
     * Eigen decomposition of 2x2 matrix.
     * @param mat 2x2 matrix
     * @return  Vector of Eigen pairs.
     */
    template <class T>
    static std::vector<EigenPair> eigen2x2(const Matrix<T>& mat);

    /**
     * Converges to Eigen pair depending on initial Eigen vector and Eigen value.
     * @param mat Matrix of which to perform Eigen decomposition.
     * @param initialEigenVector Initial Eigen vector
     * @param initialEigenValue Initial Eigen value
     * @param maxIteration Maximum number of rayleigh iterations. If precision not reached, no Eigen pair was found.
     * @param precision Rayleigh iteration stops when the Eigen value change is below the passed precision value
     * @return Eigen pair.
     */
    template <class T>
    static EigenPair rayleighIteration(const Matrix<T>& mat, const Matrix<double>& initialEigenVector, double initialEigenValue, size_t maxIteration, double precision);

    /**
     * Converges to Eigen pair with most significant Eigen value.
     * @param mat  Matrix of which to perform Eigen decomposition.
     * @param maxIteration Maximum number of power iterations. If precision not reached, no Eigen pair was found.
     * @param precision Power iteration stops when the Eigen value change is below the passed precision value
     * @return Eigen pair
     */
    template <class T>
    static EigenPair powerIteration(const Matrix<T>& mat, size_t maxIteration, double precision);

    /**
     * The QR algorithm iteratively finds all Eigen values and Eigen vectors of a matrix
     * by applying iteratively the QR decomposition.
     * @param mat  Matrix of which to perform Eigen decomposition.
     * @param maxIteration Maximum number of power iterations. If precision not reached, no Eigen pair was found.
     * @param precision Power iteration stops when the Eigen value change is below the passed precision value
     * @param showProgress Prints the progress of the algorithm
     * @return Eigen pairs
     */
    template <class T>
    static std::vector<EigenPair> qrAlgorithm(const Matrix<T>& mat, size_t maxIteration, double precision, bool showProgress = false);

    /**
     * Compute Rayleigh quotient of a matrix and a vector. This can be
     * used to find the Eigenvalue to a corresponding Eigenvector and
     * matrix.
     * @param mat
     * @param vec
     * @return Rayleigh quotient
     */
    template <class T>
    static double rayleighQuotient(const Matrix<T>& mat, const Matrix<T>& vec);

    /**
     * Computes the Householder vector and the corresponding scalar, which forms
     * together the Householder matrix P. P reflects the vector x.
     * @param x Vector.
     * @return Householder result.
     */
    template <class T>
    static HouseholderResult householder(const Matrix<T>& x);

    /**
     * Generate Householder matrix P from Householder vector and
     * corresponding scalar value, such that
     * P = I - b*v*v'
     * @param v Householder vector.
     * @param b Scalar value.
     * @return Householder matrix.
     */
    template <class T>
    static Matrix<double> householderMatrix(const Matrix<T>& v, double b);

    /**
     * Bidiagonalization of a Matrix A, so that
     * U'*A*V = B,
     * where B is bidiagonal.
     * @param a Input matrix.
     * @return Diagonalization result.
     */
    template <class T>
    static DiagonalizationResult bidiagonalization(const Matrix<T>& a);

    enum QRMethod
    {
        Householder, /* Householder reflection */
        Givens       /* Givens rotation */
    };

    /**
     * QR decomposition, where the passed matrix A is decomposed into an
     * orthogonal matrix Q and an upper triangle matrix R, such that
     * A = Q * R.
     * @param mat Matrix A
     * @param positive If true, the diagonal elements of R are positive. This needs additional computation.
     * @param method Applied method to compute the QR decomposition.
     * @return QR decomposition
     */
    template <class T>
    static QRResult qr(const Matrix<T>& mat, bool positive = true, QRMethod method = QRMethod::Householder);

    template <class T>
    static QRResult qr_householder(const Matrix<T>& mat, bool positive = true);

    template <class T>
    static QRResult qr_givens(const Matrix<T>& mat, bool positive = true);

    /**
     * RQ decomposition, where the passed matrix A is decomposed into an
     * upper triangle matrix R and an orthogonal matrix Q, such that
     * A = R * Q.
     * @param mat Matrix A
     * @param positive If true, the diagonal elements of R are positive. This needs additional computation.
     * @param method Applied method to compute the QR decomposition.
     * @return QR decomposition
     */
    template <class T>
    static QRResult rq(const Matrix<T>& mat, QRMethod method = QRMethod::Householder);

    /**
     * The QR decomposition is not unique. Changing the sign of a row and
     * its corresponding column leads to another valid decomposition. This
     * function modifies an existing QR decomposition and returns another
     * solution.
     * @param q Orthogonal matrix
     * @param r Upper triangle matrix
     * @param row Which row in r is modified (change sign)
     * @return  Another QR decomposition
     */
    template <class T>
    static QRResult qrSignModifier(const Matrix<T>& q, const Matrix<T>& r, size_t row);

    /**
     * Performs a singular value decomposition of the
     * passed matrix mat.
     * @param mat Passed matrix.
     * @return SVD decomposition
     */
    template <class T>
    static SVDResult svd(const Matrix<T>& mat);

    /**
     * Performs a singular value decomposition of the
     * passed matrix mat. The method applied is the
     * Golub and Kahan algorithm.
     * @param mat Passed matrix.
     * @return SVD decomposition
     */
    template <class T>
    static SVDResult svdGolubKahan(const Matrix<T>& mat);

    struct SvdStepResult
    {
        Matrix<double> u;
        Matrix<double> v;
    };

    static SvdStepResult svdStepGolubKahan(Matrix<double>& b)
    {
        // b is an upper bidiagonal matrix

        size_t m = b.rows();
        size_t n = b.cols();

        // lower right 2x2 submatrix - eigen decomposition
        Matrix<double> t = b.transpose() * b; // todo: speed up by forming t_minor directly
        Matrix<double> t_minor = t.subMatrix(t.rows() - 2, t.cols() - 2, 2, 2);

        std::vector<Decomposition::EigenPair> eig = Decomposition::eigen(t_minor);

        // choose eigen value closer to tnn
        double l;
        double l1 = eig.at(0).L;
        double l2 = eig.at(1).L;
        double tnn = t_minor(1, 1);
        if (std::abs(tnn - l1) < std::abs(tnn - l2))
            l = l1;
        else
            l = l2;

        // use eigenvalue to perform the first Givens rotation
        double y = t(0, 0) - l;
        double z = t(0, 1);

        Matrix<double> vR = Matrix<double>::identity(n);
        Matrix<double> uL = Matrix<double>::identity(m);

        for (size_t k = 0; k < n - 1; k++)
        {
            Matrix<double> b_in = b;

            Matrix<double> gr = givensRotationRowDirection(y, z, n, k, k, k + 1);
            b = b * gr;
            vR = vR * gr;

            y = b(k, k);
            z = b(k + 1, k);
            Matrix<double> gl = givensRotatioColumnDirection(y, z, m, k, k, k + 1);
            b = gl * b;
            uL = gl * uL;

            if (k < n - 1)
            {
                // set the values for the next column
                y = b(k, k + 1);
                z = b(k, k + 2);
            }
        }

        SvdStepResult res;
        res.v = vR;
        res.u = uL;

        // print out diagonal elements and their upper diagonal element
        /*
        for(size_t k = 0; k < n; k++ )
        {
            if( k < n-1 )
                std::cout << b(k,k) << "  ( " << b(k,k+1) << " )" << std::endl;
            else
                std::cout << b(k,k) << "  ( - )" << std::endl;
        }

        std::cout << std::endl;*/

        return res;
    }

    /*
     *  Example matrix n = 5
     *
     *      d   0   0   0   0
     *      0   d   f   0   0
     *      0   0   d   0   0
     *      0   0   0   d   0
     *      0   0   0   0   d
     *
     *    |  p |       |  q   |
     */
    static void svdCheckMatrixGolubKahan(const Matrix<double>& mat, size_t& p, size_t& q)
    {
        double eps = std::numeric_limits<double>::epsilon() * 10;
        size_t n = mat.cols();

        // find q - diagonality from back
        q = 0;
        for( size_t i = 0; i < n; i++ )
        {
            size_t x = (n - 1) - i;
            size_t y = (n - 1) - i - 1;

            if( x == 0 )
            {
                // first diagonal entry  does not have an upper element
                q++;
                break;
            }
            else if( std::abs(mat(y, x)) < eps )
            {
                // still diagonal
                q++;
            }
            else
            {
                // stop - not anymore diagonal, as upper diagonal element not zero.
                break;
            }
        }

        // find q - diagonality from back
        p = n - q;
        for( size_t i = q; i < n; i++ )
        {
            size_t x = (n - 1) - i;
            size_t y = (n - 1) - i - 1;

            if( x == 0 )
            {
                // first diagonal entry  does not have an upper element
                p--;
                break;
            }
            else if( std::abs(mat(y, x)) >= eps )
            {
                // still not diagonal again
                p--;
            }
            else
            {
                // stop - upper element is 0. However, this belongs to B22 too
                p--;
                break;
            }
        }

    }

    /*
    *  Example matrix n = 5 (p=2,q=1)
    *
    *      1   0   0   0   0
    *      0   1   0   0   0
    *      0   0   m   m   0
    *      0   0   m   m   0
    *      0   0   0   0   1
    *
    *    |   p   |       | q |
    */
    static Matrix<double> svdPaddingRotation(const Matrix<double>& mat, size_t p, size_t q)
    {
        Matrix<double> padded = Matrix<double>::identity( mat.cols() + p + q );
        padded.setSubMatrix(p,p,mat);
        return padded;
    }



    // Chapter 8, p. 490 -> If dk = 0 for some k < n, then premultiplication by a sequence of
    //                      Givens transformations can zero fk.
    static Matrix<double> svdZeroRow(Matrix<double>& b, size_t row)
    {
        size_t n = b.cols();

        Matrix<double> u_left = Matrix<double>::identity(b.rows());
        for( size_t i = row ; i < (n-1); i++ )
        {
            Matrix<double> uS = Decomposition::givensRotatioColumnDirection(b, i+1, i+1, row);
            b = uS * b;
            u_left = uS * u_left;
        }

        return u_left;
    }

    static SvdStepResult svdHandleZeroDiagonalEntries( Matrix<double>& b, size_t p, size_t q, bool& modified)
    {
        double eps = std::numeric_limits<double>::epsilon() * 10;
        size_t n = b.cols();

        for( size_t r = p; r < n - q && !modified; r++ )
        {
            // if any diagonal entry in B22 is zero, then zero the
            // superdiagonal entry in the same row.
            if( std::abs(b(r,r)) < eps )
            {
                if( r < (n-1) )
                {
                    Matrix<double> extraULeft = svdZeroRow(b, r);
                }
                else
                {
                    // special case: last diagonal element is zero
                    // "If dn = 0, then the last column can be zeroed with a series of column rotations in planes
                    //(n - 1 , n) , (n - 2, n) , . . . , ( 1 , n) . Thus, we can decouple if Ji · · · fn-1 = 0 or di · · · dn =
                    //o."

                }

                modified = true;
            }
        }

    };

    /**
     * Computes the given rotation based on a and b as input. See
     * the article https://en.wikipedia.org/wiki/Givens_rotation#Stable_calculation
     * @param a First number
     * @param b Second number -> Givens Rotation will zero b!
     * @return Givens Rotation parameters
     */
    static GivensRotation givensRotation(double a, double b)
    {
        // based on https://en.wikipedia.org/wiki/Givens_rotation#Stable_calculation

        GivensRotation res;
        if (std::abs(b) > std::numeric_limits<double>::epsilon())
        {
            double h = std::hypot(a, b); // no precession problem
            double d = 1.0 / h;

            res.C = std::abs(a) * d;
            res.S = -(std::copysign(d, a) * b);
            res.R = std::copysign(1.0, a) * h;
        }
        else
        {
            res.S = 0.0;
            res.C = 1.0;
            res.R = a;
        }

        return res;
    }

    /**
     * Generates a Givens rotation matrix for the passed values a and b
     * and the two row indices.
     * @param a First number
     * @param b Second number -> Givens Rotation will zero b!
     * @param m Matrix size
     * @param col Plane index
     * @param a_row Plane a index
     * @param b_row Plane b index, which will be set to zero.
     * @return Givens rotation.
     */
    static Matrix<double> givensRotatioColumnDirection(double a, double b, size_t m, size_t col, size_t a_row, size_t b_row)
    {
        Decomposition::GivensRotation gr = givensRotation(a, b);

        Matrix<double> gMat = Matrix<double>::identity(m);
        gMat(a_row, a_row) = gr.C;
        gMat(b_row, b_row) = gr.C;
        gMat(a_row, b_row) = -gr.S;
        gMat(b_row, a_row) = gr.S;

        return gMat;
    }

    /**
     * Generates a Givens rotation matrix for the passed matrix
     * mat and the two row indices. The Givens rotation rotates the
     * element mat(b_row,col) to zero, by G*mat. Note: a_row < b_row || a_row >= col
     * @param mat Input matrix
     * @param col Plane index
     * @param a_row Plane a index
     * @param b_row Plane b index, which will be set to zero.
     * @return
     */
    template <class T>
    static Matrix<double> givensRotatioColumnDirection(const Matrix<T> &mat, size_t col, size_t a_row, size_t b_row);

    /**
     * Generates a Givens rotation matrix for the passed values a and b
     * and the two column indices.
     * @param a First number
     * @param b Second number -> Givens Rotation will zero b!
     * @param m Matrix size
     * @param row Plane index
     * @param a_col Plane a index
     * @param b_col Plane b index, which will be set to zero.
     * @return Givens rotation.
     */
    static Matrix<double> givensRotationRowDirection(double a, double b, size_t m, size_t row, size_t a_col, size_t b_col)
    {
        Decomposition::GivensRotation gr = givensRotation(a, b);

        Matrix<double> gMat = Matrix<double>::identity(m);
        gMat(a_col, a_col) = gr.C;
        gMat(b_col, b_col) = gr.C;
        gMat(a_col, b_col) = gr.S;
        gMat(b_col, a_col) = -gr.S;

        return gMat;
    }

    /**
     * Generates a Givens rotation matrix for the passed matrix
     * mat, the row  and the two column indices. The Givens rotation rotates the
     * element mat(row,b_col) to zero, by mat*G. Note: a_col < b_col || a_col >= row
     * @param mat Input matrix
     * @param row Plane index
     * @param a_col Plane a index
     * @param b_col Plane b index, which will be set to zero.
     * @return
     */
    template <class T>
    static Matrix<double> givensRotationRowDirection(const Matrix<T> &mat, size_t row, size_t a_col, size_t b_col);

private:
    template <class T>
    static LUResult doolittle(const Matrix<T>& a, bool pivoting);

};

// Infos from:
//  http://mathonline.wikidot.com/the-algorithm-for-doolittle-s-method-for-lu-decompositions

template <class T>
Decomposition::LUResult Decomposition::luDecomposition(const Matrix<T>& mat, bool pivoting)
{
    if (mat.rows() != mat.cols())
    {
        std::cout << "Square matrix required";
        std::exit(-1);
    }

    return doolittle(mat, pivoting);
}

template <class T>
Decomposition::LUResult Decomposition::doolittle(const Matrix<T>& aIn, bool pivoting)
{
    Matrix<double> u = Matrix<double>(aIn);
    size_t         n = u.rows();
    Matrix<double> l = Matrix<double>(n, n);

    l.setToIdentity();

    Matrix<double> p           = Matrix<double>::identity(n);
    size_t         nbrRowSwaps = 0;

    for (size_t k = 0; k < n; k++)
    {
        if (pivoting)
        {
            // find pivot row and swap
            size_t pivotRow   = k;
            double pivotValue = 0.0;
            for (size_t searchPivotIdx = k; searchPivotIdx < n; searchPivotIdx++)
            {
                double cPivotElement = std::abs(u(searchPivotIdx, k));
                if (cPivotElement > pivotValue)
                {
                    pivotRow   = searchPivotIdx;
                    pivotValue = cPivotElement;
                }
            }

            // swap row if different from k
            if (pivotRow != k)
            {
                nbrRowSwaps++;

                // swap rows in u
                u.swapRows(pivotRow, k);
                p = Multiplier::swapRow(u, pivotRow, k) * p;

                //swap the subdiagonal entries of in l
                for (size_t lswapIdx = 0; lswapIdx < k; lswapIdx++)
                {
                    double tmpVal         = l(pivotRow, lswapIdx);
                    l(pivotRow, lswapIdx) = l(k, lswapIdx);
                    l(k, lswapIdx)        = tmpVal;
                }
            }
        }

        // process beneath rows, so that first pivot column element of u is zero
        double pivot = u(k, k);

        if (std::abs(pivot) > std::numeric_limits<double>::min()) // check if pivot is not zero
        {
            Matrix<double> pivotRow = u.row(k);
            for (size_t i = k + 1; i < n; i++)
            {
                double cFactor = -(u(i, k) / pivot);

                // modify row in u
                u.setRow(i, cFactor * pivotRow + u.row(i));

                // modify corresponding entry in l
                l(i, k) = -cFactor;
            }
        }
        else
        {
            // This is a singular matrix, therefore l(i,k) can be freely chosen -> lu decomposition is not unique
            // U does not be to be modified in this step
            // info: https://math.stackexchange.com/questions/2010470/doolittle-transformation-is-non-unique-for-singular-matrices
            for (size_t i = k + 1; i < n; i++)
            {
                l(i, k) = 0;
            }
        }
    }

    LUResult ret(l, u, p, nbrRowSwaps);
    return ret;
}



template <class T>
std::vector<Decomposition::EigenPair> Decomposition::eigen(const Matrix<T>& mat, Decomposition::EigenMethod method)
{
    if (!mat.isSquare())
    {
        std::cout << "Eigen: Square matrix required";
        std::exit(-1);
    }

    std::vector<EigenPair> pairs;

    // special solution for 2x2 matrix
    if( mat.rows() == 2 )
    {
        pairs = eigen2x2(mat);
    }
    else
    {
        Matrix<double> cMat(mat);

        if (cMat.isSymmetric())
        {
            switch (method)
            {
                case QRAlgorithm:
                    // QR algorithm
                    pairs = qrAlgorithm(cMat, 200, std::numeric_limits<double>::epsilon(), false);
                    break;

                case PowerIterationAndHotellingsDeflation:
                    // Power iteration and hotelling's deflation
                    // http://www.robots.ox.ac.uk/~sjrob/Teaching/EngComp/ecl4.pdf
                    for (size_t i = 0; i < cMat.rows(); i++)
                    {
                        EigenPair ePair = powerIteration(cMat, 50, std::numeric_limits<double>::epsilon());
                        if (ePair.Valid)
                        {
                            pairs.push_back(ePair);

                            // Hotelling's deflation -> remove found Eigen pair from cMat
                            cMat = cMat - (ePair.L * ePair.V * ePair.V.transpose());
                        }
                    }
                    break;
            }
        }
        else
        {
            std::cout
                << "Warning: The Eigen decomposition of non-symmetric matrices does not work yet properly in this library!!"
                << std::endl;

            // compute the largest eigenvalue
            pairs.push_back(powerIteration(cMat, 30, std::numeric_limits<double>::epsilon()));
        }
    }

    sortDescending(pairs);

    return pairs;
}

// http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/
template <class T>
std::vector<Decomposition::EigenPair> Decomposition::eigen2x2(const Matrix<T>& mat)
{
    if( mat.rows() != 2 || mat.cols() != 2 )
        throw InvalidInputException();

    double a = mat(0,0);
    double b = mat(0,1);
    double c = mat(1,0);
    double d = mat(1,1);

    double t = a + d;
    double k = a*d - b*c;

    double ll = std::sqrt(std::pow(t,2.0)/4.0-k);
    double l1 = t/2.0 + ll;
    double l2 = t/2.0 - ll;
    Matrix<double> v1(2,1);
    Matrix<double> v2(2,1);

    if( c != 0.0 )
    {
        v1(0,0) = l1 - d;
        v1(1,0) = c;

        v2(0,0) = l2 - d;
        v2(1,0) = c;
    }
    else if( b != 0.0 )
    {
        v1(0,0) = b;
        v1(1,0) = l1-a;

        v2(0,0) = b;
        v2(1,0) = l2-a;
    }
    else
    {
        l1 = a;
        v1(0,0) = 1.0;
        v1(1,0) = 0.0;

        l2 = d;
        v2(0,0) = 0.0;
        v2(1,0) = 1.0;
    }

    std::vector<Decomposition::EigenPair> pairs;
    pairs.push_back( Decomposition::EigenPair(v1,l1,true));
    pairs.push_back( Decomposition::EigenPair(v2,l2,true));

    return pairs;
}

template <class T>
Decomposition::EigenPair Decomposition::rayleighIteration(const Matrix<T>& mat, const Matrix<double>& initialEigenVector, double initialEigenValue, size_t maxIteration, double precision)
{
    Matrix<double> matD = mat;

    // Rayleigh quotient iteration: https://en.wikipedia.org/wiki/Rayleigh_quotient_iteration

    Matrix<double> e_vec        = initialEigenVector;
    double         e_val        = initialEigenValue;
    double         e_val_before = e_val;

    // init used variables
    Matrix<double> ident           = Matrix<double>::identity(matD.cols());
    Matrix<double> e_vec_unscaled  = Matrix<double>(matD.cols(), 1);
    bool           go              = true;
    size_t         nbrOfIterations = 0;
    bool           validEigenPair  = false;

    while (go)
    {
        e_vec_unscaled = (matD - (ident * e_val)).adjugate() * e_vec;
        e_vec          = e_vec_unscaled.normalizeColumns();
        e_val          = rayleighQuotient(matD, e_vec);

        //std::cout << "Iteration: " << nbrOfIterations << std::endl << e_vec << std::endl << e_val << "--------------" << std::endl;

        // check stopping criteria of Rayleigh iteration
        if (std::abs(e_val - e_val_before) < precision * std::abs(e_val + e_val_before))
        {
            go             = false;
            validEigenPair = true;
        }
        else if (nbrOfIterations >= maxIteration)
        {
            go             = false;
            validEigenPair = false;
        }

        e_val_before = e_val;
        nbrOfIterations++;
    }

    return EigenPair(e_vec, e_val, validEigenPair);
}

template <class T>
std::vector<Decomposition::EigenPair> Decomposition::qrAlgorithm(const Matrix<T>& mat, size_t maxIteration, double precision, bool showProgress)
{
    // https://en.wikipedia.org/wiki/QR_algorithm
    std::vector<EigenPair> ret;

    Matrix<double> a = mat;

    size_t         nbrOfIterations = 0;
    bool           go              = true;
    bool           foundEig        = false;
    Matrix<double> q_before        = Matrix<double>(a.rows(), a.rows());
    q_before.fill(0);
    Matrix<double> qProd = Matrix<double>::identity(a.rows());

    while (go)
    {
        Decomposition::QRResult qr = Decomposition::qr(a, false); // note: not important to have positive elements on diagonal of R

        // check stopping criteria
        if (q_before.compare(qr.Q, true, precision))
        {
            // q changed less than required precission.
            // q is good choice since its normalized
            go       = false;
            foundEig = true;
        }
        else if (nbrOfIterations >= maxIteration)
        {
            go       = false;
            foundEig = false;
        }
        else
        {
            // continue -> prepare next loop

            a     = qr.R * qr.Q;  // iteratively converge to a - diag(a) are eigenvalues
            qProd = qProd * qr.Q; // accumlate q transformations to get eigenvectors

            if (showProgress)
            {
                std::cout << "qrAlgorithm progress = " << static_cast<double>(nbrOfIterations) / static_cast<double>(maxIteration) << std::endl;
            }

            q_before = qr.Q;
            nbrOfIterations++;
        }
    }

    // The eigenvalues are in a and eigenvectors in qProd
    for (size_t i = 0; i < qProd.cols(); i++)
    {
        EigenPair eP(qProd.column(i), a(i, i), foundEig);
        ret.push_back(eP);
    }

    return ret;
}

template <class T>
Decomposition::EigenPair Decomposition::powerIteration(const Matrix<T>& mat, size_t maxIteration, double precision)
{
    // Power iteration : https://en.wikipedia.org/wiki/Power_iteration
    // Finds most significant eigenvalue

    Matrix<double> matD = mat;
    Matrix<double> eVec(matD.rows(), 1);
    eVec.fill(1.0);
    double eVal       = 1;
    double eValBefore = eVal;

    bool   validEigenPair  = false;
    size_t nbrOfIterations = 0;
    bool   go              = true;

    while (go)
    {
        eVec = (matD * eVec).normalizeColumns();
        eVal = rayleighQuotient(matD, eVec);

        // check stopping criteria of power iteration
        if (std::abs(eVal - eValBefore) < precision * std::abs(eVal + eValBefore))
        {
            go             = false;
            validEigenPair = true;
        }
        else if (nbrOfIterations >= maxIteration)
        {
            go             = false;
            validEigenPair = false;
        }

        eValBefore = eVal;
        nbrOfIterations++;
    }

    return EigenPair(eVec, eVal, validEigenPair);
}

// info: https://www.mathematik.uni-wuerzburg.de/~borzi/RQGradient_Chapter_10.pdf
template <class T>
double Decomposition::rayleighQuotient(const Matrix<T>& m, const Matrix<T>& v)
{
    Matrix<T> vT = v.transpose();
    return static_cast<double>((vT * m * v)(0, 0)) / static_cast<double>((vT * v)(0, 0));
}

template <class T>
Decomposition::HouseholderResult Decomposition::householder(const Matrix<T>& x)
{
    // generate basis vector [1.0, 0, 0, ..]
    Matrix<double> e(x.rows(), 1);
    e.fill(0.0);
    e(0, 0) = 1.0;

    double sign = 1.0;
    if (x(0, 0) > 0.0)
        sign = -1.0;

    // vh -> Householder vector
    Matrix<double> v = x - (e * x.norm() * sign);

    // scalar
    double b = 2.0 / (v.transpose() * v)(0, 0);

    return HouseholderResult(b, v);
}

template <class T>
Matrix<double> Decomposition::householderMatrix(const Matrix<T>& v, double b)
{
    Matrix<double> p = Matrix<double>::identity(v.rows());

    if (std::isfinite(b))
        p = p - b * v * v.transpose();

    return p;
}

// Sources:
// Householder bidiagonalization, Matrix computation, 4th ed, Golub & Loan, p.284
// documents/bidiagonalization.pdf -> Martin Plesinger
template <class T>
Decomposition::DiagonalizationResult Decomposition::bidiagonalization(const Matrix<T>& a_m)
{
    size_t m = a_m.rows();
    size_t n = a_m.cols();

    if (m < n)
    {
        std::cout << "bidiagonalization: Invalid matrix dimension";
        std::exit(-1);
    }

    Matrix<double> a = a_m;
    Matrix<double> u = Matrix<double>::identity(m);
    Matrix<double> v = Matrix<double>::identity(n);

    for (size_t j = 0; j < n; j++)
    {
        // row direction
        Matrix<double>    x_r     = a.subMatrix(j, j, m - j, 1);
        HouseholderResult h_r     = householder(x_r);
        Matrix<double>    h_mat_r = householderMatrix(h_r.V, h_r.B);

        Matrix<double> a_sub_r = a.subMatrix(j, j, m - j, n - j);

        // transform a with householder matrix
        a_sub_r = h_mat_r * a_sub_r;
        a.setSubMatrix(j, j, a_sub_r); // place submatrix into a

        // concatenate householder matrix to u
        Matrix<double> H_MAT_R = Matrix<double>::identity(m);
        H_MAT_R.setSubMatrix(j, j, h_mat_r);
        u = u * H_MAT_R;

        // column direction
        if (j < n - 2)
        {
            Matrix<double>    x_c     = a.subMatrix(j, j + 1, 1, n - (j + 1));
            HouseholderResult h_c     = householder(x_c.transpose());
            Matrix<double>    h_mat_c = householderMatrix(h_c.V, h_c.B);

            Matrix<double> a_sub_c = a.subMatrix(j, j + 1, m - j, n - j - 1);
            a_sub_c                = a_sub_c * h_mat_c;
            a.setSubMatrix(j, j + 1, a_sub_c); // place submatrix into a

            // concatenate householder matrix to v
            Matrix<double> H_MAT_C = Matrix<double>::identity(n);
            H_MAT_C.setSubMatrix(j + 1, j + 1, h_mat_c);
            v = v * H_MAT_C;
        }
    }

    return DiagonalizationResult(u, a, v);
}

// QR decomposition by using Householder reflection -> see documents/qr_decomposition.pdf
template <class T>
Decomposition::QRResult Decomposition::qr(const Matrix<T>& mat, bool positive, QRMethod method)
{
    switch (method)
    {
        case Householder:
            return qr_householder(mat, positive);

        case Givens:
            return qr_givens(mat, positive);
    }
}

// QR decomposition by using Householder reflection -> see documents/qr_decomposition.pdf
template <class T>
Decomposition::QRResult Decomposition::qr_householder(const Matrix<T>& mat, bool positive)
{
    Matrix<double> r = mat;
    size_t         m = r.rows();
    size_t         n = r.cols();

    // initialize q as identity
    Matrix<double> q = Matrix<double>::identity(m);

    for (size_t i = 0; i < n; i++)
    {
        // copy current column
        Matrix<double>    x     = r.subMatrix(i, i, m - i, 1);
        HouseholderResult house = householder(x);
        Matrix<double>    h     = householderMatrix(house.V, house.B);

        // create the matrix for next loop.
        Matrix<double> rSub = h * r.subMatrix(i, i, m - i, n - i);
        r.setSubMatrix(i, i, rSub);

        // use the smaller Householder matrix h to
        // create one of the right size, H.
        // each H is used to get step by step to to
        // the matrix Q.
        Matrix<double> H = Matrix<double>::identity(m);
        H.setSubMatrix(i, i, h);
        q = q * H;
    }

    QRResult retResult(q, r);

    if (positive)
    {
        // There exist multiple qr solutions. To get a unique result,
        // the diagonal elements on r are chosen to be positive.
        for (size_t i = 0; i < m; i++)
        {
            if (retResult.R(i, i) < 0)
            {
                retResult = qrSignModifier(retResult.Q, retResult.R, i);
            }
        }
    }

    return retResult;
}

// QR decomposition by using Givens rotations -> see documents/qr_decomposition.pdf
template <class T>
Decomposition::QRResult Decomposition::qr_givens(const Matrix<T>& mat, bool positive)
{
    Matrix<double> r = mat;
    size_t         m = r.rows();
    size_t         n = r.cols();

    // initialize q as identity
    Matrix<double> q = Matrix<double>::identity(m);

    for (size_t j = 0; j < n; j++)
    {
        for (size_t i = m - 1; i > j; i--)
        {
            Matrix<double> gs = givensRotatioColumnDirection(r, j, i - 1, i);
            r                 = gs * r;
            q                 = q * gs.transpose();
        }
    }

    QRResult res(q, r);

    if (positive)
    {
        // There exist multiple qr solutions. To get a unique result,
        // the diagonal elements on r are chosen to be positive.
        for (size_t i = 0; i < m; i++)
        {
            if (res.R(i, i) < 0)
            {
                res = qrSignModifier(res.Q, res.R, i);
            }
        }
    }

    return res;
}

template <class T>
Decomposition::QRResult Decomposition::rq(const Matrix<T>& mat, QRMethod method)
{
    // info: https://math.stackexchange.com/questions/1640695/rq-decomposition
    size_t m = mat.rows();

    Matrix<double> p = Matrix<double>(m, m);
    p.fill(0.0);
    for (size_t i = 0; i < m; i++)
    {
        p(i, m - 1 - i) = 1.0;
    }

    Matrix<double> ad = p * mat; // reverses the rows in mat

    Decomposition::QRResult qrD = Decomposition::qr(ad.transpose(), method);

    Matrix<double> q = p * qrD.Q.transpose();
    Matrix<double> r = p * qrD.R.transpose() * p;

    return QRResult(q, r);
}

template <class T>
Decomposition::QRResult Decomposition::qrSignModifier(const Matrix<T>& q, const Matrix<T>& r, size_t row)
{
    // https://math.stackexchange.com/questions/2237262/is-there-a-correct-qr-factorization-result

    Matrix<double> Q = q;
    Matrix<double> R = r;

    for (size_t i = 0; i < Q.rows(); i++)
        Q(i, row) = -Q(i, row);

    for (size_t i = 0; i < R.cols(); i++)
        R(row, i) = -R(row, i);

    return QRResult(std::move(Q), std::move(R));
}

#include "solve.hpp"

template <class T>
Decomposition::SVDResult Decomposition::svd(const Matrix<T>& mat)
{
    // U*S*V
    Matrix<double> aTa = mat.transpose() * mat;

    // aTa are symmetric
    std::vector<EigenPair> ep = eigen(aTa, QRAlgorithm);

    // compute singular values -> S
    // compute left orthogonal matrix  > U
    Matrix<double> s_diag_mat = Matrix<double>(mat.rows(), mat.cols());
    s_diag_mat.fill(0.0);

    Matrix<double> u_left = Matrix<double>(s_diag_mat.rows(), s_diag_mat.rows());
    u_left.fill(0.0);
    Matrix<double> v_right = Matrix<double>(s_diag_mat.cols(), s_diag_mat.cols());
    v_right.fill(0.0);

    bool forceOrthogonalization = false;

    for (size_t p = 0; p < ep.size(); p++)
    {
        double         tEigenValue = ep.at(p).L; // left and right eigenvalues should be equal
        Matrix<double> eVect       = ep.at(p).V;

        // eigenvalues of a symmetric matrix cannot be negative. if so,
        // this comes from rounding errors in eigen algorithm.
        if (tEigenValue <= 0.0)
        {
            if ((ep.size() - 1) == p)
            {
                // this is allowed to happen in the very last run
                s_diag_mat(p, p) = 0.0;
                v_right.setColumn(p, eVect);
                forceOrthogonalization = true;
                // the last column of u will be built up in the post-processing orthogonalization
            }
            else
            {
                throw SVDFailedException();
            }
        }
        else
        {
            s_diag_mat(p, p) = std::sqrt(tEigenValue);
            v_right.setColumn(p, eVect);
            u_left.setColumn(p, mat * eVect * (1.0 / s_diag_mat(p, p)));
        }
    }

    // in case of small singular values, the matrix u might be not orthogonal.
    // if so, the last column can be modified.
    if (forceOrthogonalization || !u_left.isOrthogonal(0.001))
    {
        if (u_left.rows() == 3)
        {
            // use crossproduct to make third column
            double ax = u_left(0, 0);
            double ay = u_left(1, 0);
            double az = u_left(2, 0);

            double bx = u_left(0, 1);
            double by = u_left(1, 1);
            double bz = u_left(2, 1);

            u_left(0, 2) = ay * bz - az * by;
            u_left(1, 2) = az * bx - ax * bz;
            u_left(2, 2) = ax * by - ay * bx;
        }
        else
        {
            size_t lastColumn = (ep.size() - 1);
            for (size_t u_row = 0; u_row < u_left.rows(); u_row++)
                u_left(u_row, lastColumn) = std::sqrt(1 - u_left.row(u_row).normSquare());
        }
    }

    return SVDResult(u_left, s_diag_mat, v_right);
}

template <class T>
Matrix<double> Decomposition::givensRotatioColumnDirection(const Matrix<T> &mat, size_t col, size_t a_row, size_t b_row)
{
    double b = mat(b_row, col);
    double a = mat(a_row, col);

    return givensRotatioColumnDirection(a, b, mat.rows(), col, a_row, b_row);
}

template <class T>
Matrix<double> Decomposition::givensRotationRowDirection(const Matrix<T> &mat, size_t row, size_t a_col, size_t b_col)
{
    double b = mat(row, b_col);
    double a = mat(row, a_col);

    return givensRotationRowDirection(a, b, mat.cols(), mat.rows(), a_col, b_col);;
}

// Described in Matrix Computations, 4th edition, Golub & van Loan, p.
template <class T>
Decomposition::SVDResult Decomposition::svdGolubKahan(const Matrix<T>& mat)
{
    // U*S*V
    double eps = std::numeric_limits<double>::epsilon() * 10;

    Decomposition::DiagonalizationResult diag = Decomposition::bidiagonalization(mat);
    Matrix<double> b = diag.D;

    size_t m = b.rows();
    size_t n = b.cols();

    std::vector<Matrix<double>> uLeftV;
    std::vector<Matrix<double>> vRightV;

    // svd step
    size_t q = 0;
    while( q != n )
    {
        // go through the upper band and zero values close to zero
        size_t zeroCnt = 0;
        for(size_t r = 0; r < n-1; r++)
        {
            if(std::abs(b(r,r+1)) <= eps * (std::abs(b(r,r)) + std::abs(b(r+1,r+1)) ) )
            {
                b(r,r+1) = 0.0;
                zeroCnt++;
            }
        }

        std::cout << "Current b:" << std::endl << b << std::endl << std::endl;

        size_t p;
        svdCheckMatrixGolubKahan(b,p,q);

        if( q < n )
        {
            bool modified = false;

            // check B22 (middle matrix) for zero diagonal element
            Decomposition::SvdStepResult extraRotations = Decomposition::svdHandleZeroDiagonalEntries(b,p,q, modified);

            if( !modified )
            {
                // size_t subMatSize = n-q-p;
                // Matrix<double> b22 = b.subMatrix(p,p,subMatSize,subMatSize);
                // std::cout << "B22:" << std::endl << b22 << std::endl << std::endl;

                //todo: continue here by adapting the size of b and padd returned V and U

                Decomposition::SvdStepResult step = Decomposition::svdStepGolubKahan(b);

                vRightV.push_back(step.v);
                uLeftV.push_back(step.u);
            }
        }

    }

    Matrix<double> u_left_accum = Matrix<double>::identity(b.rows());
    Matrix<double> v_right_accum = Matrix<double>::identity(b.cols());
    for( size_t k = 0; k < uLeftV.size(); k++ )
    {
        v_right_accum = v_right_accum * vRightV.at(k);
        u_left_accum = uLeftV.at(k) * u_left_accum;
    }

    Matrix<double> u = diag.U * u_left_accum.transpose();
    Matrix<double> v = diag.V * v_right_accum;

    // Make singular values positive
    Matrix<double> sing(b.rows(), b.cols());
    sing.fill(0.0);

    for( size_t k = 0; k < b.cols(); k++ )
    {
        double s = b(k,k);

        if( s < 0.0 )
        {
            // invert s and row k in u
            s = -s;
            for( size_t r = 0; r < u.rows(); r++ )
                u(r,k) = -u(r,k);
        }

        sing(k,k) = s;
    }

    return SVDResult(u,sing,v);
}

#endif //MY_DECOMPOSITION_H
