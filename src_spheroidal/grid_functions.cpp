#include <spheroidal/grid_functions.h>
#include <yawg/utils.hpp>
#include <yawg/core.h>
#include <gsl/gsl_math.h>

/*!
 * \brief Computes grid size of spheroidal harmonics grid given the order.
 * \param p Order of the spheroidal harmonics
 * \param nu Number of points in the theta direction
 * \param nv Number of points in the lambda direction
 *
 * \return The order of the spheroidal harmonics (legacy usage)
 */
int spharm_grid_size_ord(int p, int &nu, int &nv)
{
    nu = p + 1;
    nv = 2 * p;
    return p;
}

/*!
 * \brief Computes grid size of spheroidal harmonics grid given the total number of points.
 * \param ntot The total number of points in the grid
 * \param nu Number of points in the theta direction
 * \param nv Number of points in the lambda direction
 *
 * \note Will cause an error if the number of points is not valid for a spheroidal harmonics grid
 * \return The order of the spheroidal harmonics (legacy usage)
 */
int spharm_grid_size_tot(int ntot, int &nu, int &nv)
{
    int p = (round(sqrt(2 * ntot + 1)) - 1) / 2;
    nu = p + 1;
    nv = 2 * p;

    // assert( ntot == nu*nv );  // TODO: Include this or equivalent
    return p;
}

/*!
 * \brief Populates the matrices \p U and \p V with a spheroidal harmonics grid of order \p
 * \param p Order of the spheroidal harmonics
 * \param U Reference to gsl::matrix of Gaussian spaced rows
 * \param V Reference to gsl::matrix of Uniform spaced columns
 *
 * Mapped to a spheroid, \p u is spaced on [0, pi] and \p v is spaced on [0, 2pi].
 */
void gl_grid(size_t p, gsl::matrix &U, gsl::matrix &V)
{
    int nu, nv;
    spharm_grid_size_ord(p, nu, nv);

    gsl::vector lambda = gsl::linspace(0, 2 * M_PI, nv);
    gsl::vector theta = gsl::arrayfun([](double x)
                                      { return acos(x); },
                                      gsl::leggauss(nu));
    gsl::meshgrid(theta, lambda, U, V);
}