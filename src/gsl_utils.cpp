#include <fmt/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*!
\brief Pretty(ish) print the elements of a gsl_vector to stdout
*/
void gsl_vector_print( gsl_vector* v )
{
    fmt::print( "[" );
    for( int i = 0; i < v->size; ++i )
        fmt::print( (i == 0) ? "{:g}" : ", {:g}", gsl_vector_get(v, i) );
    fmt::print( "]\n" );
}

/*!
\brief Pretty(ish) print the elements of a gsl_matrix to stdout
*/
void gsl_matrix_print( gsl_matrix* m )
{
    for( int i = 0; i < m->size1; ++i )
    {
        fmt::print( (i == 0) ? "[" : " " );
        for( int j = 0; j < m->size2; ++j )
            fmt::print( (j == 0) ? "{: 9g}" : " {: 9g}", gsl_matrix_get(m, i, j) );
        fmt::print( (i == (m->size1 - 1) ? "]\n" : "\n") );
    }
}
