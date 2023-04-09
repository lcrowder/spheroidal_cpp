#ifndef GSL_UTILS_H_
#define GSL_UTILS_H_

#include <spheroidal/gsl_wrapper.h>
#include <iostream>
#include <functional>
#include <utility>

namespace gsl
{
    /*! \brief Store Gauss-Legendre quadrature nodes and weights */
    void leggauss(size_t n, gsl::vector &x, gsl::vector &w, double a = -1.0, double b = 1.0);

    /*! \brief Get a vector Gauss-Legendre quadrature nodes */
    vector leggauss(size_t n, double a = -1.0, double b = 1.0);

    /*! \brief Get a gsl::vector of evenly spaced points on the interval [a, b] (inclusive) */
    gsl::vector linspace(double a, double b, size_t N = 100);

    /*! \brief Store 2D grid coordinates based on 1D input gsl::vectors */
    void meshgrid(const gsl::vector &x, const gsl::vector &y, gsl::matrix &X, gsl::matrix &Y);

    /*! \brief Return the nxn identity matrix */
    gsl::matrix eye(size_t n);
    
    /*! \brief Compute arccos of each element in array */
    gsl::vector acos(const gsl::vector &x);

    //! \brief Move version of acos
    gsl::vector acos(gsl::vector &&x);

    /*! \brief Apply lambda function to each element of a gsl::vector, akin to MATLAB arrayfun
    * \param func Lambda function to apply to each element
    * \param x Input matrix
    * 
    * \note This function is not optimized for speed, but rather for convenience,
    *  as there is overhead in using templated Lambda functions. If speed is necessary,
    *  add the function to gsl_utils specifically. 
    * \return Matrix of same size as \p x with \p func applied to each element
    */
    template <typename Lambda>
    gsl::vector arrayfun(Lambda&& func, const gsl::vector &x)
    {
        gsl::vector y(x.size());
        for (int i = 0; i < x.size(); ++i)
            y(i) = func(x(i));
        return y;

    }

    //! \brief Move version of arrayfun for vectors
    template <typename Lambda>
    gsl::vector arrayfun(Lambda&& func, gsl::vector &&x)
    {
        for (int i = 0; i < x.size(); ++i)
            x(i) = func(x(i));
        return std::move(x);
    }

    /*! \brief Apply lambda function to each element of a gsl::matrix, akin to MATLAB arrayfun
    * \param func Lambda function to apply to each element
    * \param x Input matrix
    * 
    * \note This function is not optimized for speed, but rather for convenience,
    *  as there is overhead in using templated Lambda functions. If speed is necessary,
    *  add the function to gsl_utils specifically. 
    * \return Matrix of same size as \p x with \p func applied to each element
    */
    template <typename Lambda>
    gsl::matrix arrayfun(Lambda&& func, const gsl::matrix &x)
    {
        gsl::matrix y(x.nrows(), x.ncols());
        for (int i = 0; i < x.nrows(); ++i)
            for (int j = 0; j < x.ncols(); ++j)
                y(i, j) = func(x(i, j));
        return y;
    }

    //! \brief Move version of arrayfun
    template <typename Lambda>
    gsl::matrix arrayfun(Lambda&& func, gsl::matrix &&x)
    {
        for (int i = 0; i < x.nrows(); ++i)
            for (int j = 0; j < x.ncols(); ++j)
                x(i, j) = func(x(i, j));
        return std::move(x);
    }

} // namespace gsl

#endif // GSL_UTILS_H_