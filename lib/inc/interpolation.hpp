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

#ifndef MY_INTERPOLATION_H
#define MY_INTERPOLATION_H

#include "matrix.hpp"
#include <limits>

class Interpolation1D
{
public:

    /**
     * Interpolates the function at position x.
     * @param x Position
     * @return Function value
     */
    virtual double interpolate(double x) const = 0;

    /**
     * Interpolates the function at each position in x.
     * @param xV Positions in form of n x 1 vector
     * @return Function values y and positions x in form of n x 2 matrix (x,y)
     */
    virtual Matrix<double> interpolate(const Matrix<double>& xV) const
    {
        Matrix<double> ip( xV.rows(), 2 );
        for( size_t m = 0; m < xV.rows(); m++ )
        {
            double x = xV(m,0);
            ip(m,0) = x;
            ip(m,1) = interpolate(x);
        }

        return ip;
    }

protected:

    struct range1D
    {
        double startR;
        double endR;
    };

    virtual bool isInRange(range1D range, double value) const
    {
        return range.startR <= value && range.endR >= value;
    }
};

class LinearInterpolation : public Interpolation1D
{
public:

    /**
     * Constructor of linear interpolation.
     * @param data n x 2 data matrix, with n >= 2, where each row is a data point (x,y)
     */
    template <typename T>
    LinearInterpolation(const Matrix<T>& data)
    {
        prepare(data);
    }

    virtual ~LinearInterpolation(){}

    using Interpolation1D::interpolate;

    double interpolate(double x) const override
    {
        // find last function which is covering the range of x
        linfunc f = m_linf.front();
        for( size_t k = 1; k < m_linf.size(); k++ )
        {
            linfunc lf = m_linf.at(k);
            if( isInRange(lf.r,x) )
                f = lf;
            else
                break;
        }

        // compute interpolation result
        double y = f.a * x + f.b;
        return y;
    }

private:

    template <typename T>
    void prepare( const Matrix<T>& interpdata )
    {
        // Sort following rising x-vals;
        Matrix<double> data( interpdata );
        data.sortRows(0,Matrix<double>::SortDirection::Ascending);

        if( data.rows() < 2 || data.cols() < 2 )
            throw new InvalidInputException();

        for( size_t i = 0; i < data.rows() - 1; i++ )
        {
            double x0 = data(i,0); double y0 = data(i,1);
            double x1 = data(i+1,0); double y1 = data(i+1,1);

            linfunc f;
            f.a = (y0-y1)/(x0-x1);
            f.b = y1 - f.a * x1;

            // special treatment for first function
            if( i == 0 )
                f.r.startR = std::numeric_limits<double>::lowest();
            else
                f.r.startR = x0;

            // not used in this linear interpolation
            f.r.endR = std::numeric_limits<double>::max();

            m_linf.push_back(f);
        }
    }

    struct linfunc
    {
        range1D r;
        double a = 0.0;
        double b = 0.0;
    };

    std::vector<linfunc> m_linf;
};

// from https://en.wikipedia.org/wiki/Cubic_Hermite_spline
class CubicSplineInterpolation : public Interpolation1D
{
public:

    /**
     * Constructor of cubic spline interpolation.
     * @param data n x 2 data matrix, with n >= 3, where each row is a data point (x,y).
     */
    template <typename T>
    CubicSplineInterpolation(const Matrix<T>& data)
    {
        prepare(data);
    }

    virtual ~CubicSplineInterpolation(){}

    using Interpolation1D::interpolate;

    double interpolate(double x) const override
    {
        // find spline
        bool foundSpline;
        spline sp;
        for( spline s : m_spline )
        {
            if( isInRange(s.r, x ) )
            {
                sp = s;
                foundSpline = true;
                break;
            }
        }

        if( !foundSpline )
            throw OutOfRangeException();

        // compute interpolation
        double h = sp.r.endR-sp.r.startR;
        double t = (x-sp.r.startR)/h;
        double res = h00(t)*sp.pk + h10(t)*h*sp.mk + h01(t)*sp.pkp1 + h11(t)*h*sp.mkp1;
        return res;
    }

private:

    template <typename T>
    void prepare( const Matrix<T>& interpdata )
    {
        // Sort following rising x-vals;
        Matrix<double> data( interpdata );
        data.sortRows(0,Matrix<double>::SortDirection::Ascending);

        if( data.rows() < 3 || data.cols() < 2 )
            throw new InvalidInputException();

        // compute tangents
        std::vector<double> mk;
        for( size_t k = 0; k < data.rows(); k++ )
            mk.push_back( finite_difference_tangent(data, k) );

        m_spline.clear();
        for(size_t k = 0; k < mk.size()-1; k++)
        {
            spline sp;
            sp.mk = mk.at(k);
            sp.mkp1 = mk.at(k+1);
            sp.pk = data(k,1);
            sp.pkp1 = data(k+1,1);
            sp.r.startR = data(k,0);
            sp.r.endR = data(k+1,0);
            m_spline.push_back(sp);
        }
    }

    double finite_difference_tangent( const Matrix<double>& data, size_t k ) const
    {
        if( k >= data.rows() )
            throw new InvalidInputException();

        double p_k = data(k, 1);
        double t_k = data(k, 0);

        if( k == 0 )
        {
            // first segment - one sided
            double p_kp1 = data(k+1, 1);
            double t_kp1 = data(k+1, 1);
            return (p_kp1 - p_k)/(t_kp1 - t_k);
        }
        else if( k > 0 && k < data.rows() - 1)
        {
            // in between
            double p_kp1 = data(k+1, 1);
            double p_km1 = data(k-1, 1);
            double t_kp1 = data(k+1, 0);
            double t_km1 = data(k-1, 0);

            return 0.5*( (p_kp1 - p_k)/(t_kp1 - t_k) + (p_k - p_km1)/(t_k - t_km1) );
        }
        else
        {
            // last segment - one sided
            double p_km1 = data(k-1, 1);
            double t_km1 = data(k-1, 1);
            return (p_k - p_km1)/(t_k - t_km1);
        }
    }

    // Hermite basis functions
    double h00(double t) const
    {
        return 2.0*std::pow(t,3.0) - 3.0*std::pow(t,2.0) + 1;
    }
    double h10(double t) const
    {
        return std::pow(t,3.0) - 2.0*std::pow(t,2.0) + t;
    }
    double h01(double t) const
    {
        return -2.0*std::pow(t,3.0) + 3.0*std::pow(t,2.0);
    }
    double h11(double t) const
    {
        return std::pow(t,3.0) - std::pow(t,2.0);
    }

    struct spline
    {
            double mk;
            double mkp1;
            double pk;
            double pkp1;
            range1D r;
    };

    std::vector<spline> m_spline;
};

#endif // MY_INTERPOLATION_H

