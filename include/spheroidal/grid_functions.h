#ifndef GRID_FUNCTIONS_H_
#define GRID_FUNCTIONS_H_

#include <spheroidal/gsl_wrapper.h>

//! \brief Computes grid size of spheroidal harmonics grid given the order.
int spharm_grid_size_ord( int p, int& nu, int& nv );

//! \brief Computes grid size of spheroidal harmonics grid given the total number of points.
int spharm_grid_size_tot( int ntot, int& nu, int& nv );

//! \brief Populates the matrices \p U and \p V with a spheroidal harmonics grid of order \p
void gl_grid( size_t p, gsl::matrix& U, gsl::matrix& V);

#endif // GRID_FUNCTIONS_H_