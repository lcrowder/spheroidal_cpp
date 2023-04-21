#ifndef YAWG_LEGENDRE_H_
#define YAWG_LEGENDRE_H_

#include <yawg/core.h>
#include <gsl/gsl_sf_legendre.h>

namespace gsl
{
    enum class legendre_norm : int
    {
        none = GSL_SF_LEGENDRE_NONE,
        schmidt = GSL_SF_LEGENDRE_SCHMIDT,
        spharm = GSL_SF_LEGENDRE_SPHARM,
        full = GSL_SF_LEGENDRE_FULL
    };
    
    vector legendre_P(size_t nmax, const vector& x, legendre_norm mode = legendre_norm::none );

    double spherical_harmonic( int n, int m, double x );

    complex spherical_harmonic(int n, int m, double theta, double phi);
    
}

#endif // YAWG_LEGENDRE_H_