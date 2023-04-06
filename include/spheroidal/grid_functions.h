#ifndef GRID_FUNCTIONS_H_
#define GRID_FUNCTIONS_H_

#include <spheroidal/gsl_wrapper.h>

int spharm_grid_size_ord( int p, int& nu, int& nv );
int spharm_grid_size_tot( int ntot, int& nu, int& nv );
void gl_grid( size_t p, gsl::matrix& U, gsl::matrix& V);

#endif // GRID_FUNCTIONS_H_