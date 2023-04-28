#ifndef SPHEROIDAL_DOUBLE_LAYER_H_
#define SPHEROIDAL_DOUBLE_LAYER_H_

#include <yawg/core.h>

gsl::matrix solid_harmonic(int p, gsl::vector u_x, int region=0);
void DLspectrum(int p, double u0, gsl::vector &lambda_int, gsl::vector &lambda_surf, gsl::vector &lambda_ext);
gsl::cmatrix Ynm_matrix(int p, gsl::vector v, gsl::vector phi);
gsl::cmatrix spheroidal_double_layer(gsl::cmatrix sigma, double u0, gsl::matrix X, int target_coords=0);
gsl::cmatrix spheroidal_double_layer(gsl::cmatrix sigma, double u0);
gsl::cvector spheroidal_double_layer(gsl::cvector sigma, double u0, gsl::matrix X, int target_coords=0);
gsl::cvector spheroidal_double_layer(gsl::cvector sigma, double u0);

#endif // _SPHEROIDAL_DOUBLE_LAYER_H_