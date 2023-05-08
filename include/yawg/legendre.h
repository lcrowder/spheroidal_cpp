#ifndef YAWG_LEGENDRE_H_
#define YAWG_LEGENDRE_H_

#include <yawg/core.h>
#include <gsl/gsl_sf_legendre.h>

namespace gsl
{
    /*! \brief Wrapper enum for the GSL gsl_sf_legendre_t */
    enum class legendre_norm : int
    {
        none = GSL_SF_LEGENDRE_NONE,
        schmidt = GSL_SF_LEGENDRE_SCHMIDT,
        spharm = GSL_SF_LEGENDRE_SPHARM,
        full = GSL_SF_LEGENDRE_FULL
    };

    //! \brief Compute the normalized associated Legendre polynomial P_n^m(x)
    double spherical_harmonic( int n, int m, double x );

    //! \brief Compute the normalized associated spherical harmonic Y_n^m(theta, phi)
    complex spherical_harmonic(int n, int m, double theta, double phi);
    
}

#endif // YAWG_LEGENDRE_H_