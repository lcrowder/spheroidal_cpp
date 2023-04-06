#ifndef GSL_UTILS_H_
#define GSL_UTILS_H_

#include <spheroidal/gsl_wrapper.h>

/*! \brief Store Gauss-Legendre quadrature nodes and weights */
void leggauss( size_t n, gsl::vector& x, gsl::vector& w, double a=-1.0, double b=1.0 );

/*! \brief Get a vector Gauss-Legendre quadrature nodes */
gsl::vector leggauss( size_t n, double a=-1.0, double b=1.0 );

/*! \brief Get a vector of evenly spaced points on the interval [a, b] (inclusive) */
gsl::vector linspace( double a, double b, size_t N=100);

/*! \brief Store 2D grid coordinates based on 1D input vectors */
void meshgrid( const gsl::vector& x, const gsl::vector& y, gsl::matrix& X, gsl::matrix& Y );

/*! \brief Compute arccos of each element in array */
gsl::vector acos( const gsl::vector& x );

#endif // GSL_UTILS_H_