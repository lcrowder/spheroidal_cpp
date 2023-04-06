#include <spheroidal/grid_functions.h>
#include <spheroidal/gsl_utils.h>
#include <spheroidal/gsl_wrapper.h>
#include <fmt/core.h>

int main()
{
    gsl::vector vec(5);
    
    vec.print( );

    for( int i = 0; i < 5; ++i )
        vec(i) = i;

    vec.print( );

    gsl::vector vec2( vec );
    vec2.print( );

    const int p = 5;

    gsl::matrix U, V; 
    gl_grid( p, U, V );

    fmt::print("U =\n");
    U.print( );

    fmt::print("V =\n");
    V.print( );

    /*
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
    */
    return 0;
    
}