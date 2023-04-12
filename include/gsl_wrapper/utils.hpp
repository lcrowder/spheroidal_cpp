#ifndef GSL_WRAPPER_UTILS_HPP_
#define GSL_WRAPPER_UTILS_HPP_

#include <gsl_wrapper/core.h>
#include <utility>

namespace gsl
{
    /*! \brief Store Gauss-Legendre quadrature nodes and weights */
    void leggauss(size_t n, gsl::vector &x, gsl::vector &w, double a = -1.0, double b = 1.0);

    /*! \brief Get a vector Gauss-Legendre quadrature nodes */
    vector leggauss(size_t n, double a = -1.0, double b = 1.0);

    /*! \brief Get a gsl::vector of evenly spaced points on the interval [a, b] (inclusive) */
    gsl::vector linspace(double a, double b, size_t N = 100);

    //! \brief Complex version of linspace
    gsl::cvector linspace(gsl::complex a, gsl::complex b, size_t n);

    gsl::vector arange(double a, double b, double step = 1.0);

    /*! \brief Store 2D grid coordinates based on 1D input gsl::vectors */
    void meshgrid(const gsl::vector &x, const gsl::vector &y, gsl::matrix &X, gsl::matrix &Y);

    /*! Complex version of gsl::meshgrid */
    void meshgrid(const gsl::cvector &x, const gsl::cvector &y, gsl::cmatrix &X, gsl::cmatrix &Y);

    /*! \brief Return the nxn identity matrix */
    gsl::matrix eye(size_t n);

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
    gsl::vector arrayfun(Lambda &&func, const gsl::vector &x)
    {
        gsl::vector y(x.size());
        for (int i = 0; i < x.size(); ++i)
            y(i) = func(x(i));
        return y;
    }

    //! \brief Move version of arrayfun for vectors
    template <typename Lambda>
    gsl::vector arrayfun(Lambda &&func, gsl::vector &&x)
    {
        for (int i = 0; i < x.size(); ++i)
            x(i) = func(x(i));
        return std::move(x);
    }

    //! \brief Copy version of arrayfun for complex vectors
    template <typename Lambda>
    gsl::cvector arrayfun(Lambda &&func, const gsl::cvector &x)
    {
        gsl::cvector y(x.size());
        for (int i = 0; i < x.size(); ++i)
            //y(i) = func(x(i)); // Desired usage
            y.set( func(x.get(i)), i );
        return y;
    }

    //! \brief Move version of arrayfun for complex vectors
    template <typename Lambda>
    gsl::cvector arrayfun(Lambda &&func, gsl::cvector &&x)
    {
        for (int i = 0; i < x.size(); ++i)
            //x(i) = func(x(i)); // Desired usage
            x.set( func(x.get(i)), i );
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
    gsl::matrix arrayfun(Lambda &&func, const gsl::matrix &x)
    {
        gsl::matrix y(x.nrows(), x.ncols());
        for (int i = 0; i < x.nrows(); ++i)
            for (int j = 0; j < x.ncols(); ++j)
                y(i, j) = func(x(i, j));
        return y;
    }

    //! \brief Move version of arrayfun
    template <typename Lambda>
    gsl::matrix arrayfun(Lambda &&func, gsl::matrix &&x)
    {
        for (int i = 0; i < x.nrows(); ++i)
            for (int j = 0; j < x.ncols(); ++j)
                x(i, j) = func(x(i, j)); 
        return std::move(x);
    }

    //! \brief Copy version of arrayfun for complex matrices
    template <typename Lambda>
    gsl::cmatrix arrayfun(Lambda &&func, const gsl::cmatrix &x)
    {
        gsl::cmatrix y(x.nrows(), x.ncols());
        for (int i = 0; i < x.nrows(); ++i)
            for (int j = 0; j < x.ncols(); ++j)
                //y(i, j) = func(x(i, j)); // Desired usage
                y.set( func(x.get(i, j)), i, j );
        return y;
    }

    //! \brief Move version of arrayfun for complex matrices
    template <typename Lambda>
    gsl::cmatrix arrayfun(Lambda &&func, gsl::cmatrix &&x)
    {
        for (int i = 0; i < x.nrows(); ++i)
            for (int j = 0; j < x.ncols(); ++j)
                //x(i, j) = func(x(i, j));// Desired usage
                x.set( func(x.get(i, j)), i, j );
        return std::move(x);
    }

} // namespace gsl

#endif // GSL_WRAPPER_UTILS_HPP_