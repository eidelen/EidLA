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
        Matrix<double> L; // Lower triangle matrix
        Matrix<double> U; // Upper triangle matrix
        Matrix<double> P; // Row Swaps
        size_t NbrRowSwaps;     // Number of row swaps
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
        std::sort(pairs.begin(), pairs.end(), [](const EigenPair& a, const EigenPair& b)
        {
            return a.L > b.L;
        });
    }

    static void sortAscending(std::vector<EigenPair>& pairs)
    {
        std::sort(pairs.begin(), pairs.end(), [](const EigenPair& a, const EigenPair& b)
        {
            return a.L < b.L;
        });
    }

    static void filter(std::vector<EigenPair>& pairs, double eigenValueThreshold)
    {
        std::vector<EigenPair> copyEP = pairs;
        pairs.clear();
        std::copy_if (copyEP.begin(), copyEP.end(), std::back_inserter(pairs), [eigenValueThreshold](EigenPair i){return i.L >= eigenValueThreshold; } );
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
            if(&other == this)
                return *this;

            this->Q = other.Q;
            this->R = other.R;
            return *this;
        }

        QRResult& operator=(QRResult&& other)
        {
            // check for self-assignment
            if(&other == this)
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
        HouseholderResult(double b, Matrix<double> v) :
        B(b), V(v)
        {
        }
        double B;         // scalar
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
        QRAlgorithm, // QR algorithm
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


    /**
     * QR decomposition, where the passed matrix A is decomposed into an
     * orthogonal matrix Q and an upper triangle matrix R, such that
     * A = Q * R.
     * @param mat Matrix A
     * @param positive If true, the diagonal elements of R are positive. This needs additional computation.
     * @return QR decomposition
     */
    template <class T>
    static QRResult qr(const Matrix<T>& mat, bool positive = true);

    /**
     * RQ decomposition, where the passed matrix A is decomposed into an
     * upper triangle matrix R and an orthogonal matrix Q, such that
     * A = R * Q.
     * @param mat Matrix A
     * @param positive If true, the diagonal elements of R are positive. This needs additional computation.
     * @return QR decomposition
     */
    template <class T>
    static QRResult rq(const Matrix<T>& mat);


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
    static SVDResult svd( const Matrix<T> mat );



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

    Matrix<double> p = Matrix<double>::identity(n);
    size_t nbrRowSwaps = 0;

    for (size_t k = 0; k < n; k++)
    {
        if( pivoting )
        {
            // find pivot row and swap
            size_t pivotRow = k;
            double pivotValue = 0.0;
            for (size_t searchPivotIdx = k; searchPivotIdx < n; searchPivotIdx++)
            {
                double cPivotElement = std::abs(u(searchPivotIdx, k));
                if (cPivotElement > pivotValue)
                {
                    pivotRow = searchPivotIdx;
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
                for( size_t lswapIdx = 0; lswapIdx < k; lswapIdx++)
                {
                    double tmpVal = l(pivotRow, lswapIdx);
                    l(pivotRow, lswapIdx) = l(k,lswapIdx);
                    l(k,lswapIdx) = tmpVal;
                }
            }
        }


        // process beneath rows, so that first pivot column element of u is zero
        double pivot = u(k,k);

        if( std::abs(pivot) > std::numeric_limits<double>::min()) // check if pivot is not zero
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
            for(size_t i = k + 1; i < n; i++)
            {
                l(i, k) = 0;
            }
        }
    }

    LUResult ret(l, u, p, nbrRowSwaps);
    return ret;
}

template <class T>
std::vector<Decomposition::EigenPair> Decomposition::eigen(const Matrix<T>& mat, Decomposition::EigenMethod method )
{
    if(!mat.isSquare())
    {
        std::cout << "Eigen: Square matrix required";
        std::exit(-1);
    }

    Matrix<double> cMat(mat);
    std::vector<EigenPair> pairs;

    if( cMat.isSymmetric() )
    {
        switch( method )
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
                    EigenPair ePair = powerIteration(cMat, 50,std::numeric_limits<double>::epsilon());
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
        std::cout << "Warning: The Eigen decomposition of non-symmetric matrices does not work yet properly in this library!!" << std::endl;

        // compute the largest eigenvalue
        pairs.push_back(powerIteration(cMat, 30,std::numeric_limits<double>::epsilon()));
    }

    sortDescending(pairs);

    return pairs;
}

template <class T>
Decomposition::EigenPair Decomposition::rayleighIteration(const Matrix<T>& mat, const Matrix<double>& initialEigenVector, double initialEigenValue, size_t maxIteration, double precision)
{
    Matrix<double> matD      = mat;

    // Rayleigh quotient iteration: https://en.wikipedia.org/wiki/Rayleigh_quotient_iteration

    Matrix<double> e_vec = initialEigenVector;
    double e_val = initialEigenValue;
    double e_val_before = e_val;

    // init used variables
    Matrix<double> ident = Matrix<double>::identity(matD.cols());
    Matrix<double> e_vec_unscaled = Matrix<double>(matD.cols(),1);
    bool go = true;
    size_t nbrOfIterations = 0;
    bool validEigenPair = false;

    while (go)
    {
        e_vec_unscaled = (matD - (ident*e_val) ).adjugate() * e_vec;
        e_vec = e_vec_unscaled.normalizeColumns();
        e_val = rayleighQuotient(matD,e_vec);

        //std::cout << "Iteration: " << nbrOfIterations << std::endl << e_vec << std::endl << e_val << "--------------" << std::endl;

        // check stopping criteria of Rayleigh iteration
        if( std::abs(e_val - e_val_before) < precision*std::abs(e_val + e_val_before) )
        {
            go = false;
            validEigenPair = true;
        }
        else if( nbrOfIterations >= maxIteration )
        {
            go = false;
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

    size_t nbrOfIterations = 0;
    bool go = true;
    bool foundEig = false;
    Matrix<double> q_before = Matrix<double>(a.rows(), a.rows());
    q_before.fill(0);
    Matrix<double> qProd = Matrix<double>::identity(a.rows());

    while (go)
    {
        Decomposition::QRResult qr = Decomposition::qr(a, false); // note: not important to have positive elements on diagonal of R

        // check stopping criteria
        if( q_before.compare(qr.Q,true,precision))
        {
            // q changed less than required precission.
            // q is good choice since its normalized
            go = false;
            foundEig = true;
        }
        else if( nbrOfIterations >= maxIteration )
        {
            go = false;
            foundEig = false;
        }
        else
        {
            // continue -> prepare next loop

            a = qr.R * qr.Q; // iteratively converge to a - diag(a) are eigenvalues
            qProd = qProd * qr.Q; // accumlate q transformations to get eigenvectors

            if( showProgress )
            {
                std::cout << "qrAlgorithm progress = " << static_cast<double>(nbrOfIterations) / static_cast<double> (maxIteration) << std::endl;
            }

            q_before = qr.Q;
            nbrOfIterations++;
        }
    }

    // The eigenvalues are in a and eigenvectors in qProd
    for( size_t i = 0; i < qProd.cols(); i++ )
    {
        EigenPair eP( qProd.column(i), a(i,i), foundEig );
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
    Matrix<double> eVec(matD.rows(),1); eVec.fill(1.0);
    double eVal = 1;
    double eValBefore = eVal;

    bool validEigenPair = false;
    size_t nbrOfIterations = 0;
    bool go = true;

    while (go)
    {
        eVec = (matD * eVec).normalizeColumns();
        eVal = rayleighQuotient(matD,eVec);

        // check stopping criteria of power iteration
        if( std::abs(eVal - eValBefore) < precision*std::abs(eVal + eValBefore) )
        {
            go = false;
            validEigenPair = true;
        }
        else if( nbrOfIterations >= maxIteration )
        {
            go = false;
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
    Matrix<double> e(x.rows(),1);
    e.fill(0.0);
    e(0,0) = 1.0;

    double sign = 1.0;
    if( x(0,0) > 0.0 )
        sign = -1.0;

    // vh -> Householder vector
    Matrix<double> v = x - (e*x.norm()*sign);

    // scalar
    double b = 2.0 / (v.transpose() * v)(0,0);

    return HouseholderResult(b,v);
}

template <class T>
Matrix<double> Decomposition::householderMatrix(const Matrix<T>& v, double b)
{
    Matrix<double> p = Matrix<double>::identity(v.rows());

    if( std::isfinite(b) )
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

    if( m < n )
    {
        std::cout << "bidiagonalization: Invalid matrix dimension";
        std::exit(-1);
    }

    Matrix<double> a = a_m;
    Matrix<double> u = Matrix<double>::identity(m);
    Matrix<double> v = Matrix<double>::identity(n);

    for( size_t j = 0; j < n; j++)
    {
        // row direction
        Matrix<double> x_r = a.subMatrix(j,j, m-j, 1);
        HouseholderResult h_r = householder(x_r);
        Matrix<double> h_mat_r = householderMatrix(h_r.V, h_r.B);

        Matrix<double> a_sub_r = a.subMatrix(j,j, m-j, n-j);

        // transform a with householder matrix
        a_sub_r = h_mat_r * a_sub_r;
        a.setSubMatrix(j,j, a_sub_r); // place submatrix into a

        // concatenate householder matrix to u
        Matrix<double> H_MAT_R = Matrix<double>::identity(m);
        H_MAT_R.setSubMatrix(j,j,h_mat_r);
        u = u * H_MAT_R;

        // column direction
        if( j < n-2 )
        {
            Matrix<double> x_c = a.subMatrix(j,j+1, 1, n-(j+1));
            HouseholderResult h_c = householder(x_c.transpose());
            Matrix<double> h_mat_c = householderMatrix(h_c.V, h_c.B);

            Matrix<double> a_sub_c = a.subMatrix(j,j+1, m-j, n-j-1);
            a_sub_c = a_sub_c * h_mat_c;
            a.setSubMatrix(j, j+1, a_sub_c); // place submatrix into a

            // concatenate householder matrix to v
            Matrix<double> H_MAT_C = Matrix<double>::identity(n);
            H_MAT_C.setSubMatrix(j+1,j+1,h_mat_c);
            v =  v * H_MAT_C;
        }
    }

    return DiagonalizationResult(u,a,v);
}

// QR decomposition by using Householder reflection -> see documents/qr_decomposition.pdf
template <class T>
Decomposition::QRResult Decomposition::qr(const Matrix<T>& mat, bool positive)
{
    Matrix<double> r = mat;
    size_t m = r.rows();
    size_t n = r.cols();

    // initialize q as identity
    Matrix<double> q = Matrix<double>::identity(m);

    for( size_t i = 0; i < n; i++ )
    {
        // copy current column
        Matrix<double> x = r.subMatrix(i,i,m-i,1);
        HouseholderResult house = householder(x);
        Matrix<double> h = householderMatrix(house.V, house.B);

        // create the matrix for next loop.
        Matrix<double> rSub = h * r.subMatrix(i,i, m-i, n-i);
        r.setSubMatrix(i,i,rSub);

        // use the smaller Householder matrix h to
        // create one of the right size, H.
        // each H is used to get step by step to to
        // the matrix Q.
        Matrix<double> H = Matrix<double>::identity(m);
        H.setSubMatrix(i,i,h);
        q = q * H;
    }

    QRResult retResult(std::move(q),std::move(r));

    if( positive )
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

    return std::move(retResult);
}

template <class T>
Decomposition::QRResult Decomposition::rq(const Matrix<T>& mat)
{
    // info: https://math.stackexchange.com/questions/1640695/rq-decomposition
    size_t m = mat.rows();

    Matrix<double> p = Matrix<double>(m, m);
    p.fill(0.0);
    for( size_t i = 0; i < m; i++ )
    {
        p(i,m-1-i) = 1.0;
    }

    Matrix<double> ad = p * mat; // reverses the rows in mat

    Decomposition::QRResult qrD = Decomposition::qr(ad.transpose());

    Matrix<double> q = p * qrD.Q.transpose();
    Matrix<double> r = p * qrD.R.transpose() * p;

    return QRResult(q,r);
}

template <class T>
Decomposition::QRResult Decomposition::qrSignModifier(const Matrix<T>& q, const Matrix<T>& r, size_t row)
{
    // https://math.stackexchange.com/questions/2237262/is-there-a-correct-qr-factorization-result

    Matrix<double> Q = q;
    Matrix<double> R = r;

    for( size_t i = 0; i < Q.rows(); i++)
        Q(i, row) = -Q(i, row);

    for( size_t i = 0; i < R.cols(); i++)
        R(row,i) = -R(row,i);

    return QRResult(std::move(Q),std::move(R));
}

#include "solve.hpp"

template <class T>
Decomposition::SVDResult Decomposition::svd(const Matrix<T> mat)
{
    // U*S*V
    Matrix<double> aTa = mat.transpose()*mat;

    // aTa are symmetric
    std::vector<EigenPair> ep = eigen(aTa, QRAlgorithm);

    // compute singular values -> S
    // compute left orthogonal matrix  > U
    Matrix<double> s_diag_mat = Matrix<double>(mat.rows(), mat.cols());
    s_diag_mat.fill(0.0);

    Matrix<double> u_left = Matrix<double>(s_diag_mat.rows(), s_diag_mat.rows()); u_left.fill(0.0);
    Matrix<double> v_right = Matrix<double>(s_diag_mat.cols(), s_diag_mat.cols()); v_right.fill(0.0);

    bool forceOrthogonalization = false;

    for( size_t p = 0; p < ep.size(); p++ )
    {
        double tEigenValue = ep.at(p).L; // left and right eigenvalues should be equal
        Matrix<double> eVect = ep.at(p).V;

        // eigenvalues of a symmetric matrix cannot be negative. if so,
        // this comes from rounding errors in eigen algorithm.
        if( tEigenValue <= 0.0 )
        {
            if( (ep.size() - 1) == p )
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
            u_left.setColumn( p, mat * eVect * (1.0 / s_diag_mat(p,p)));
        }
    }

    // in case of small singular values, the matrix u might be not orthogonal.
    // if so, the last column can be modified.
    if( forceOrthogonalization || !u_left.isOrthogonal(0.001) )
    {
        if( u_left.rows() == 3 )
        {
            // use crossproduct to make third column
            double ax = u_left(0,0);
            double ay = u_left(1,0);
            double az = u_left(2,0);

            double bx = u_left(0,1);
            double by = u_left(1,1);
            double bz = u_left(2,1);

            u_left(0,2) = ay*bz - az*by;
            u_left(1,2) = az*bx - ax*bz;
            u_left(2,2) = ax*by - ay*bx;
        }
        else
        {
            size_t lastColumn = (ep.size() - 1);
            for(size_t u_row = 0; u_row < u_left.rows(); u_row++)
                u_left(u_row,lastColumn) = std::sqrt(1 - u_left.row( u_row ).normSquare());
        }

    }

    return SVDResult(u_left, s_diag_mat, v_right);
}


#endif //MY_DECOMPOSITION_H
