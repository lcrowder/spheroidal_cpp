#include <yawg/core.h>
#include <yawg/legendre.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>

/*! \brief Compute the normalized associated Legendre polynomial P_n^m(x)
 * \param n Degree of the polynomial
 * \param m Order of the polynomial
 * \param x Point at which to evaluate the polynomial
 *
 * \note Assumes n >= 0, |m| <= n, and -1 <= x <= 1
 * 
 * \return The value of the polynomial at \p x
 */
double gsl::spherical_harmonic(int n, int m, double x)
{
    if (m < 0)
        return (m % 2 == 0) ? gsl_sf_legendre_sphPlm(n, -m, x) : -gsl_sf_legendre_sphPlm(n, -m, x);
    else
        return gsl_sf_legendre_sphPlm(n, m, x);
}

/*! \brief Compute the normalized associated spherical harmonic Y_n^m(theta, phi)
 * \param n Degree of the polynomial
 * \param m Order of the polynomial
 * \param theta Azimuthal angle
 * \param phi Polar angle
 *
 * \note Assumes n >= 0, |m| <= n, and 0 <= theta <= pi, 0 <= phi <= 2pi
 * 
 * \return The value of the spherical harmonic at \p x
 */
gsl::complex gsl::spherical_harmonic(int n, int m, double theta, double phi)
{
    using namespace gsl::complex_literals;
    double Ynm;

    if (m < 0)
        Ynm = (m % 2 == 0) ? gsl_sf_legendre_sphPlm(n, -m, cos(theta)) : -gsl_sf_legendre_sphPlm(n, -m, cos(theta));
    else
        Ynm = gsl_sf_legendre_sphPlm(n, m, cos(theta));

    return Ynm * gsl::complex(gsl_complex_exp(1.0_i * m * phi));
}