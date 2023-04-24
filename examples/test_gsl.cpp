#include <yawg/core.h>
#include <yawg/utils.hpp>
#include <yawg/fft.h>
#include <gsl/gsl_blas.h>

int main()
{
    gsl::matrix M(9, 9);
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            M(i, j) = i + j;
    
    gsl::matrix N(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            N(i, j) = i + j;

    gsl::matrix M2(3, 3);
    
    return 0; 
}