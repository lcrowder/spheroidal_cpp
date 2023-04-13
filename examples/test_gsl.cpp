#include <gsl_wrapper/core.h>
#include <gsl_wrapper/utils.hpp>
#include <gsl_wrapper/fft.h>

int main()
{
    //gsl_vector * v = gsl_vector_alloc(10);
    gsl::cmatrix x(3, 4);
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 4; j++)
            x.set(i, j, gsl::complex(0.5 * i + 0.25 * j, 0.25 * i + 0.5 * j));

    gsl::cmatrix y = gsl::fft( x, 2 );
    y.print();
    return 0;
    
}