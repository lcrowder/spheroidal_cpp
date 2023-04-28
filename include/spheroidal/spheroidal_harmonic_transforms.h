#ifndef SPHEROIDAL_HARMONIC_TRANSFORMS_H_
#define SPHEROIDAL_HARMONIC_TRANSFORMS_H_

#include <yawg/core.h>

gsl::matrix get_legendre_matrix(int p, int m);
gsl::matrix get_legendre_matrix_inv(int p, int m);
gsl::cmatrix spheroidal_analysis(gsl::cmatrix f);
gsl::cmatrix spheroidal_snythesis(gsl::cmatrix shc);

#endif // _SPHEROIDAL_HARMONIC_TRANSFORMS_H_