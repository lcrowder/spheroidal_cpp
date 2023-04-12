#include <gsl_wrapper/core.h>
#include <gsl_wrapper/utils.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <iostream>

/*!
 * \brief Computes \p n Gauss-Legendre quadrature nodes and weights for the interval [a, b] (default [-1, 1])
 * \param n Number of nodes
 * \param x Reference to a gsl::vector of nodes
 * \param w Reference to a gsl::vector of weights
 * \param a Lower bound of the interval (Default -1)
 * \param b Upper bound of the interval (Default +1)
 *
 * \note This could be made more efficient by caching the gsl_integration_glfixed_table,
 * rather than allocating and freeing it each time.
 */
void gsl::leggauss(size_t n, gsl::vector &x, gsl::vector &w, double a, double b)
{
    x.resize(n);
    w.resize(n);

    if (n != 0)
    {
        gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(n);
        for (int i = 0; i < n; ++i)
            gsl_integration_glfixed_point(a, b, i, &x(i), &w(i), t);

        gsl_integration_glfixed_table_free(t);
    }
}
/*!
 * \brief Computes \p n Gauss-Legendre quadrature nodes and weights for the interval [a, b] (default [-1, 1])
 * \param n Number of nodes
 * \param a Lower bound of the interval (Default -1)
 * \param b Upper bound of the interval (Default +1)

 * \return gsl::vector of nodes
 */
gsl::vector gsl::leggauss(size_t n, double a, double b)
{
    gsl::vector x(n);
    double trash_w;

    if (n != 0)
    {
        gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(n);
        for (int i = 0; i < n; ++i)
            gsl_integration_glfixed_point(a, b, i, &x(i), &trash_w, t);

        gsl_integration_glfixed_table_free(t);
    }

    return x;
}

/*!
 * \brief Computes \p n evenly spaced points on the interval [\p a, \p b] (inclusive)
 * \param a Lower bound of the interval
 * \param b Upper bound of the interval
 * \param n Number of points
 *
 * \return gsl::vector of points
 */
gsl::vector gsl::linspace(double a, double b, size_t n)
{
    if (n == 0)
        return gsl::vector();
    if (n == 1)
    {
        gsl::vector x(1);
        x(0) = b;
        return x;
    }

    gsl::vector x(n);
    double dx = (b - a) / (n - 1);
    for (int i = 0; i < n; ++i)
        x(i) = a + i * dx;
    return x;
}

//! \brief Complex version of linspace
gsl::cvector gsl::linspace(gsl::complex a, gsl::complex b, size_t n)
{
    if (n == 0)
        return gsl::cvector();
    if (n == 1)
    {
        gsl::cvector x(1);
        x.set(0, b);
        return x;
    }

    gsl::cvector x(n);
    double dx = (b - a).abs() / (n - 1);
    for (int i = 0; i < n; ++i)
        x.set( i, a );//+ dx * i );
    return x;
}

/*! Do matlab thing */
void gsl::meshgrid(const gsl::vector &x, const gsl::vector &y, gsl::matrix &X, gsl::matrix &Y)
{
    X.resize(x.size(), y.size());
    Y.resize(x.size(), y.size());

    for (int i = 0; i < x.size(); ++i)
        for (int j = 0; j < y.size(); ++j)
        {
            X(i, j) = x(i);
            Y(i, j) = y(j);
        }
}

/*! Complex version of gsl::meshgrid */
void gsl::meshgrid(const gsl::cvector &x, const gsl::cvector &y, gsl::cmatrix &X, gsl::cmatrix &Y)
{
    X.resize(x.size(), y.size());
    Y.resize(x.size(), y.size());

    for (int i = 0; i < x.size(); ++i)
        for (int j = 0; j < y.size(); ++j)
        {
            //X(i, j) = x(i);
            //Y(i, j) = y(j);
            X.set( i, j, x.get(i) );
            Y.set( i, j, y.get(j) );
        }
}

/*! Return the nxn identity matrix */
gsl::matrix gsl::eye(size_t n)
{
    gsl::matrix I(n, n);
    for (int i = 0; i < n; ++i)
        I(i, i) = 1;
    return I;
}

//! \brief Compute FFT for each column of a gsl::matrix

