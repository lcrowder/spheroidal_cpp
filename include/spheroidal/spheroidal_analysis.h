#ifndef SPHEROIDAL_ANALYSIS_H_
#define SPHEROIDAL_ANALYSIS_H_

#include <yawg/core.h>

gsl::matrix get_legendre_matrix(int p, int m);
gsl::cmatrix spheroidal_analysis(gsl::matrix f);

#endif // _SPHEROIDAL_ANALYSIS_H_