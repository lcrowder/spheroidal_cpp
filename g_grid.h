#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>

void g_grid( int n, gsl_vector* x, gsl_vector* w )
{
    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(n);
    for( int i = 0; i < n; ++i )
        gsl_integration_glfixed_point( -1.0, 1.0, i, gsl_vector_ptr(x, i), gsl_vector_ptr(w, i), t );
    
    gsl_integration_glfixed_table_free(t);
}