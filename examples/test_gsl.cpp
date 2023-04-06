#include <spheroidal/grid_functions.h>
#include <spheroidal/gsl_utils.h>
#include <iostream>
#include <gsl/gsl_matrix.h>

int main()
{
    const int N = 15;
    gsl_vector* x = gsl_vector_alloc(N);
    gsl_vector* w = gsl_vector_alloc(N);

    g_grid(N, x, w);

    gsl_vector_print( x );
    gsl_vector_print( w );

    // Get the size of the grid
    int order = 5, nu, nv;
    spharm_grid_size_ord( order, nu, nv );

    gsl_matrix *u = gsl_matrix_alloc( nu, nv );
    gsl_matrix *v = gsl_matrix_alloc( nu, nv );
    gl_grid( u, v );

    fmt::print("u =\n");
    gsl_matrix_print( u );
    fmt::print("\nv =\n");
    gsl_matrix_print( v );

    gsl_vector_free(x);
    gsl_vector_free(w);
    
    return 0;
    
}