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
        const Matrix<double> L; // Lower triangle matrix
        const Matrix<double> U; // Upper triangle matrix
        const Matrix<double> P; // Row Swaps
        size_t NbrRowSwaps;     // Number of row swaps
    };

    struct EigenPair
    {
        EigenPair(Matrix<double> v, double l, bool valid)
        : V(v), L(l), Valid(valid)
        {
        }
        const Matrix<double> V;     // Eigen vector
        const double         L;     // Eigen value
        const bool           Valid; // Is eigen pair valid? It is if precision was reached.
    };

    struct QRResult
    {
        QRResult(Matrix<double> q, Matrix<double> r)
        : Q(q), R(r)
        {
        }
        const Matrix<double> Q; // Orthogonal matrix
        const Matrix<double> R; // Upper triangle matrix
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

    /**
     * Eigen decomposition of the matrix mat. Finds all Eigen pairs
     * in symmetric real matrices. In non-symmetric matrices, it finds
     * at least the most significant Eigen pair.
     * @param mat
     * @return Vector of Eigen pairs .
     */
    template <class T>
    static std::vector<EigenPair> eigen(const Matrix<T>& mat);

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
     * Compute Rayleigh quotient of a matrix and a vector. This can be
     * used to find the Eigenvalue to a corresponding Eigenvector and
     * matrix.
     * @param mat
     * @param vec
     * @return Rayleigh quotient
     */
    template <class T>
    static double rayleighQuotient(const Matrix<T>& mat, const Matrix<T> vec);

    template <class T>
    static QRResult qrDecomposition(const Matrix<T> mat);

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

/*
        std::cout << "step " << k << std::endl;
        std::cout << "l:" << l << std::endl;
        std::cout << "u:" << u << std::endl;
        std::cout << "p:" << p << std::endl;
        std::cout << "-------------------" << std::endl;*/
    }

    LUResult ret(l, u, p, nbrRowSwaps);
    return ret;
}

template <class T>
std::vector<Decomposition::EigenPair> Decomposition::eigen(const Matrix<T>& mat)
{
    if(!mat.isSquare())
    {
        std::cout << "Eigen: Square matrix required";
        std::exit(-1);
    }

    Matrix<double> cMat(mat);
    std::vector<EigenPair> ret;

    if( cMat.isSymmetric() )
    {
        // use power iteration and hotelling's deflation
        // http://www.robots.ox.ac.uk/~sjrob/Teaching/EngComp/ecl4.pdf

        for( size_t i = 0; i < cMat.rows(); i++ )
        {
            EigenPair ePair = powerIteration(cMat, 30, 10e-10);
            if (ePair.Valid)
            {
                ret.push_back(ePair);

                // Hotelling's deflation -> remove found Eigen pair from cMat
                cMat = cMat - (ePair.L * ePair.V * ePair.V.transpose());
            }
        }
    }
    else
    {
        std::cout << "Note: This Eigen decomposition on non-symmetric matrices does not work properly!!" << std::endl;

        /*
        for( int k = 0; k < 3; k++ )
        {
            // Step 1: Find most significant Eigen pair by using the power iteration
            EigenPair msPair = powerIteration(cMat, 20, std::numeric_limits<double>::epsilon());

            if (msPair.Valid)
            {

                std::cout << "Mat input: " << std::endl << cMat << std::endl;
                std::cout << "Eigenvalue = " << msPair.L << " Eigenvector = " << msPair.V.transpose() << std::endl;
                std::cout << "Diff deflated matrix: " << std::endl << (cMat*msPair.V - msPair.L*msPair.V).transpose() << std::endl;
                std::cout << "Diff original matrix: " << std::endl << (mat*msPair.V - msPair.L*msPair.V).transpose() << "----------------" << std::endl ;

                // modify based on former eigenvectors
                // http://zoro.ee.ncku.edu.tw/na/res/09-power_method.pdf

                ret.push_back(msPair);

                // Hotelling's deflation -> remove most significant Eigen pair
                cMat = cMat - (msPair.L * msPair.V * msPair.V.transpose());
            }
            else
            {
                std::cout << "Power iteration did not find an Eigen pair" << std::endl;
                return ret;
            }


            for( size_t i = 0; i < cMat.cols()*2; i++ )
    {
        // initial eigen value and eigen vector
        Matrix<double> initialEigenVector = Matrix<double>::random(cMat.cols(), 1, -1, 1);
        initialEigenVector.normalizeColumns();

        double initialEigenValue = Matrix<double>::random(1, 1, -10, +10)(0,0); // random value -/+ most significant eigen value

        EigenPair cPair = rayleighIteration(cMat, initialEigenVector, initialEigenValue, 100, std::numeric_limits<double>::epsilon()*1000);
        ret.push_back(cPair);

        // Hotelling deflation
        if(cPair.Valid)
        {
            cMat = cMat - (cPair.L * (cPair.V * cPair.V.transpose()));
        }
    }
        }*/
    }

    return ret;
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
double Decomposition::rayleighQuotient(const Matrix<T>& m, const Matrix<T> v)
{
    Matrix<T> vT = v.transpose();
    return static_cast<double>((vT * m * v)(0, 0)) / static_cast<double>((vT * v)(0, 0));
}

template <class T>
Decomposition::QRResult Decomposition::qrDecomposition(const Matrix<T> mat)
{
    Matrix<double> a = mat;
    size_t m = a.rows();
    size_t n = a.cols();

    Matrix<double> q = Matrix<double>::identity(m);

    for( size_t i = 0; i < n; i++ )
    {
        Matrix<double> x = a.subMatrix(i,i,m-i,1);
        Matrix<double> e(m-i,1);
        e(0,0) = 1.0;

        double sign = 1.0;
        if( x(0,0) > 0.0 )
            sign = -1;

        Matrix<double> vh = x - (e*x.norm()*sign);

        double c = 2.0 / (vh.transpose() * vh)(0,0);

        Matrix<double> h = Matrix<double>::identity(vh.rows()) - c * vh * vh.transpose();

        Matrix<double> aSub = h * a.subMatrix(i,i, m-i, n-i);

        // copy aSub to the region of a
        a.setSubMatrix(i,i,aSub);


        // compute part of q
        Matrix<double> H = Matrix<double>::identity(m);
        H.setSubMatrix(i,i,h);
        q = q * H;


        std::cout << i << std::endl <<  "x =" << std::endl << x << std::endl <<  "xnorm =" << x.norm() << std::endl <<  "vh =" << std::endl << vh << std::endl << "c=" << c << std::endl
                  << "h" << std::endl << h << "aSub" << std::endl << aSub << "a" << std::endl << a;

        std::cout << "-------------------" << std::endl;
    }

    std::cout << "******************************" << std::endl;
    std::cout << "Q = " << std::endl << q;
    std::cout << "R = " << std::endl << a;

}

#endif //MY_DECOMPOSITION_H
