#ifndef GRID_FUNCTIONS_H_
#define GRID_FUNCTIONS_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int spharm_grid_size_ord( int p, int& nu, int& nv );
int spharm_grid_size_tot( int ntot, int& nu, int& nv );
void g_grid( int n, gsl_vector* x, gsl_vector* w );
void gl_grid( gsl_matrix * u, gsl_matrix * v );


#endif // GRID_FUNCTIONS_H_