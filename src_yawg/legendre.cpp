#include <yawg/core.h>
#include <yawg/legendre.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <stdio.h>

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