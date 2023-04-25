#ifndef SPHEROIDAL_HARMONIC_TRANSFORMS_H_
#define SPHEROIDAL_HARMONIC_TRANSFORMS_H_

#include <yawg/core.h>

gsl::matrix get_legendre_matrix(int p, int m);
gsl::matrix get_legendre_matrix_inv(int p, int m);
int geti(int n, int m);
gsl::cmatrix spheroidal_analysis(gsl::matrix f);
gsl::cmatrix spheroidal_snythesis(gsl::cmatrix shc);

#endif // _SPHEROIDAL_HARMONIC_TRANSFORMS_H_