#include <iostream>
#include <gsl/gsl_matrix.h>
#include "grid_functions.h"

int test_gsl()
{
    const int N = 15;
    gsl_vector* x = gsl_vector_alloc(N);
    gsl_vector* w = gsl_vector_alloc(N);

    gsl_matrix *u, *v;

    gl_grid( 5, u, v );
    
    g_grid(N, x, w);

    for( int i = 0; i < N; ++i )
    {
        std::cout << "x[" << i << "] = " << gsl_vector_get(x, i) << ", w[" << i << "] = " << gsl_vector_get(w, i) << std::endl;
    }

    return 0;
    
}