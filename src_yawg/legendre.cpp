#include <yawg/core.h>
#include <yawg/legendre.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <stdio.h>

//! \brief Compute Legendre polynomials P_l(x) of degree \p l at points \p x
gsl::vector gsl::legendre_P(size_t nmax, const gsl::vector &x, legendre_norm mode)
{
    gsl_vector *P = gsl_vector_calloc(gsl_sf_legendre_nlm(nmax));
    gsl_sf_legendre_array((gsl_sf_legendre_t)mode, nmax, x(0), P->data);
    gsl::vector P_vec(P);

    return P_vec;
}

double gsl::spherical_harmonic(int n, int m, double x)
{
    if (m < 0)
        return (m % 2 == 0) ? gsl_sf_legendre_sphPlm(n, -m, x) : -gsl_sf_legendre_sphPlm( n, -m, x );
    else 
        return gsl_sf_legendre_sphPlm(n, m, x);
}

gsl::complex gsl::spherical_harmonic(int n, int m, double theta, double phi)
{
    using namespace gsl::complex_literals;
    double Ynm;

    if (m < 0)
        Ynm =(m % 2 == 0) ? gsl_sf_legendre_sphPlm(n, -m, cos(theta)) : -gsl_sf_legendre_sphPlm( n, -m, cos(theta) );
    else
        Ynm = gsl_sf_legendre_sphPlm(n, m, cos(theta));

    return Ynm * gsl::complex(gsl_complex_exp(1.0_i * m * phi));
}