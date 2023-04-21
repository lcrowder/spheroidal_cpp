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
    
    vector legendre_P(size_t lmax, const vector& x, legendre_norm mode = legendre_norm::none );
}

#endif // YAWG_LEGENDRE_H_