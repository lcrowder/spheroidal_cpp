#include <spheroidal/grid_functions.h>
#include <iostream>
#include <gsl/gsl_matrix.h>

int main()
{
    const int N = 15;
    gsl_vector* x = gsl_vector_alloc(N);
    gsl_vector* w = gsl_vector_alloc(N);

    g_grid(N, x, w);

    for( int i = 0; i < N; ++i )
    {
        std::cout << "x[" << i << "] = " << gsl_vector_get(x, i) << ", w[" << i << "] = " << gsl_vector_get(w, i) << std::endl;
    }

    return 0;
    
}