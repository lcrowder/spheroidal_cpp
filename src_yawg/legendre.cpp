#include <yawg/core.h>
#include <yawg/legendre.h>
#include <stdio.h>

//! \brief Compute Legendre polynomials P_l(x) of degree \p l at points \p x
gsl::vector gsl::legendre_P( size_t lmax, const gsl::vector& x, legendre_norm mode)
{
    gsl_vector *P = gsl_vector_calloc( gsl_sf_legendre_array_n(lmax));
    gsl_sf_legendre_array( (gsl_sf_legendre_t) mode, lmax, x(0), P->data);
    gsl::vector P_vec( P );

    return P_vec;
}