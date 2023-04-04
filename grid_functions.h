#ifndef GRID_FUNCTIONS_H
#define GRID_FUNCTIONS_H

#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int spharm_grid_size_ord( int p, int& nu, int& nv );
int spharm_grid_size_tot( int ntot, int& nu, int& nv );
void g_grid( int n, gsl_vector* x, gsl_vector* w );
void gl_grid( int p, gsl_matrix * u, gsl_matrix * v );

#endif