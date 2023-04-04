#include "grid_functions.h"

// The version of this function in the repo was confusing, so I expanded
// it to make it into two different functions for clarity.
int spharm_grid_size_ord( int p, int& nu, int& nv )
{
    nu = p+1;
    nv = 2*p;
    return p;
}

int spharm_grid_size_tot( int ntot, int& nu, int& nv )
{
    int p = ( round( sqrt( 2*ntot + 1 ) ) - 1 ) / 2;
    nu = p+1;
    nv = 2*p;

    //assert( ntot == nu*nv );  // Would have to include this
    return p;
}

void g_grid( int n, gsl_vector* x, gsl_vector* w )
{
    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(n);
    for( int i = 0; i < n; ++i )
        gsl_integration_glfixed_point( -1.0, 1.0, i, gsl_vector_ptr(x, i), gsl_vector_ptr(w, i), t );
    
    gsl_integration_glfixed_table_free(t);
}

void gl_grid( int p, gsl_matrix * u, gsl_matrix * v )
{
    int nu, nv;
    spharm_grid_size_ord( p, nu, nv );

    gsl_vector* lambda = gsl_vector_alloc(nv);
    gsl_vector* theta = gsl_vector_alloc(nu);

    g_grid( nu, theta, nullptr );

    for( int i = 0; i < nv; ++i )
        gsl_vector_set( lambda, i, M_2_PI * i / nv );


}