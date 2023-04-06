#ifndef GSL_UTILS_H_
#define GSL_UTILS_H_

#include <fmt/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*!
\brief Print the elements of a gsl_vector to stdout
*/
void gsl_vector_print( gsl_vector* v );

/*!
\brief Print the elements of a gsl_matrix to stdout
*/
void gsl_matrix_print( gsl_matrix* m );

#endif // GSL_UTILS_H_