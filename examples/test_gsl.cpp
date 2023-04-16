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
    
    M.print();
    // Loop over each column of M, reshape it to a 3x3 array, and multiply it by N
    for (int i = 0; i < 9; i++)
    {
        // Get **Reference** to row i
        gsl::row_view row = M.row(i);

        // Reshape row to 3x3 matrix.
        //  Note: This only works for "Row Views" (i.e. not "Column Views")
        //   or vector_views that have stride 1. This is due to irreconcilable memory
        //   storage differences. Apologies
        gsl::matrix_view row_reshaped( row.reshape( 3, 3 ) );

        // Multiply row_reshaped by N.
        //  Another note: Keep in mind that matrix multiplication can't be done in place,
        //  so there is a temp created here. All in all, these views may not be super useful
        //  in every case. Consider carefully if it'd be simpler just to copy the column
        row_reshaped = row_reshaped * N;
    }
    M.print();

    return 0; 
}