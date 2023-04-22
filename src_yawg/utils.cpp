#include <yawg/core.h>
#include <yawg/utils.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>

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
        x.set(i, a + dx * i);
    return x;
}

/*! \brief Get a gsl::vector of \p step spaced points on the interval [ \p a, \p b )
 *  \param a Lower bound of the interval
 *  \param b Upper bound of the interval
 *  \param step Spacing between points
 *
 *  \note Will return an empty vector if \p b cannot be reached by steping from \p a by \p step.
 *
 *  \return gsl::vector of points
 */
gsl::vector gsl::arange(double start, double stop, double step)
{
    if (step == 0)
        return gsl::vector();

    if ((stop - start) / step < 0)
        return gsl::vector();

    size_t n = ceil((stop - start) / step);
    gsl::vector x(n);
    for (int i = 0; i < n; ++i)
        x(i) = start + i * step;
    return x;
}

/*! \brief Perform a circular shift of the elements of a gsl::vector
 *  \param x Input vector
 *  \param k Number of positions to shift
 *
 *  \return Shifted vector
 */
gsl::vector gsl::circshift(const gsl::vector &x, int k)
{
    gsl::vector y(x);

    // If k is zero, don't shift
    if (k == 0)
        return x;

    // Negative shifts are equivalent to positive shifts in the opposite direction
    if (k < 0)
        k = k * (1 - x.size());

    for (int i = 0; i < x.size(); ++i)
        y((i + k) % x.size()) = x(i);
    return y;
}

/*! \brief Perform a circular shift of the elements of a gsl::cvector
 *  \param x Input vector
 *  \param k Number of positions to shift
 *
 *  \return Shifted vector
 */
gsl::cvector gsl::circshift(const gsl::cvector &x, int k)
{
    gsl::cvector y(x);

    // If k is zero, don't shift
    if (k == 0)
        return x;

    // Negative shifts are equivalent to positive shifts in the opposite direction
    if (k < 0)
        k = k * (1 - x.size());

    for (int i = 0; i < x.size(); ++i)
        y((i + k) % x.size()) = x(i);
    return y;
}

/*! \brief Store 2D grid coordinates based on 1D input gsl::vectors
 * \param x 1D vector of x-coordinates
 * \param y 1D vector of y-coordinates
 * \param X 2D matrix of x-coordinates
 * \param Y 2D matrix of y-coordinates
 */
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

/*! \brief Store 2D grid coordinates based on 1D input gsl::cvectors
 * \param x 1D complex vector of x-coordinates
 * \param y 1D complex vector of y-coordinates
 * \param X 2D complex matrix of x-coordinates
 * \param Y 2D complex matrix of y-coordinates
 */
void gsl::meshgrid(const gsl::cvector &x, const gsl::cvector &y, gsl::cmatrix &X, gsl::cmatrix &Y)
{
    X.resize(x.size(), y.size());
    Y.resize(x.size(), y.size());

    for (int i = 0; i < x.size(); ++i)
        for (int j = 0; j < y.size(); ++j)
        {
            // X(i, j) = x(i);
            // Y(i, j) = y(j);
            X.set(i, j, x.get(i));
            Y.set(i, j, y.get(j));
        }
}

/*! \brief Return the nxn identity matrix */
gsl::matrix gsl::eye(size_t n)
{
    gsl::matrix I(n, n);
    for (int i = 0; i < n; ++i)
        I(i, i) = 1;
    return I;
}


gsl::vector gsl::diag(const gsl::matrix &A)
{
    gsl::vector d(std::min(A.nrows(), A.ncols()));
    for (int i = 0; i < d.size(); ++i)
        d(i) = A(i, i);
    return d;
}

gsl::vector gsl::pow(const gsl::vector &x, double p)
{
    gsl::vector y(x.size());
    for (int i = 0; i < x.size(); ++i)
        y(i) = std::pow(x(i), p);
    return y;
}

gsl::cvector gsl::pow(const gsl::cvector &x, gsl::complex p)
{
    gsl::cvector y(x.size());
    for (int i = 0; i < x.size(); ++i)
        y(i) = gsl_complex_pow(x(i), p);
    return y;
}

gsl::matrix gsl::pow(const gsl::matrix &x, double p)
{
    gsl::matrix y(x.nrows(),x.ncols());
    for (int i = 0; i < x.nrows(); ++i)
        for (int j = 0; j < x.ncols(); ++j)
            y(i,j) = std::pow(x(i,j), p);
    return y;
}

gsl::cmatrix gsl::pow(const gsl::cmatrix &x, gsl::complex p)
{
    gsl::cmatrix y(x.nrows(),x.ncols());
    for (int i = 0; i < x.nrows(); ++i)
        for (int j = 0; j < x.ncols(); ++j)
            y(i,j) = gsl_complex_pow(x(i,j), p);
    return y;
}


