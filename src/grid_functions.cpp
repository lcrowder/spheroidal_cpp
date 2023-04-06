#include <spheroidal/gsl_utils.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <fmt/core.h>

/*!
 * \fn \brief Computes the number of points in each direction of a spheroidal harmonics grid, 
 *            given the order of the spheroidal harmonics.
 * \param p Order of the spheroidal harmonics
 * \param nu Number of points in the theta direction
 * \param nv Number of points in the lambda direction
 *
 * \return The order of the spheroidal harmonics (legacy usage)
 */
int spharm_grid_size_ord( int p, int& nu, int& nv )
{
    nu = p+1;
    nv = 2*p;
    return p;
}

/*!
 * \fn \brief Computes the number of points in each direction of a spheroidal harmonics grid,
 *            given the total number of points in the grid.
 * \param p The total number of points in the grid
 * \param nu Number of points in the theta direction
 * \param nv Number of points in the lambda direction
 *
 * \note Will cause an error if the number of points is not valid for a spheroidal harmonics grid
 * \return The order of the spheroidal harmonics (legacy usage)
 */
int spharm_grid_size_tot( int ntot, int& nu, int& nv )
{
    int p = ( round( sqrt( 2*ntot + 1 ) ) - 1 ) / 2;
    nu = p+1;
    nv = 2*p;

    //assert( ntot == nu*nv );  // TODO: Include this or equivalent
    return p;
}

/*!
 * \fn \brief Computes \p n Gauss-Legendre quadrature nodes on the interval [-1, 1]
 * \param n Number of nodes
 * \param x Pointer to a gsl_vector of nodes
 */
void g_grid( int n, gsl_vector* x )
{
    double w_temp;
    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(n);
    for( int i = 0; i < n; ++i )
        gsl_integration_glfixed_point( -1.0, 1.0, i, gsl_vector_ptr(x, i), &w_temp, t );
    
    gsl_integration_glfixed_table_free(t);
}

/*!
 * \fn \brief Computes \p n Gauss-Legendre quadrature nodes and weights for the interval [-1, 1]
 * \param n Number of nodes
 * \param x Pointer to a gsl_vector of nodes
 * \param w Pointer to a gsl_vector of weights
 */
void g_grid( int n, gsl_vector* x, gsl_vector* w )
{
    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(n);
    for( int i = 0; i < n; ++i )
        gsl_integration_glfixed_point( -1.0, 1.0, i, gsl_vector_ptr(x, i), gsl_vector_ptr(w, i), t );
    
    gsl_integration_glfixed_table_free(t);
}

/*!
 * \fn \brief Populates the matrices \p u and \p v with a meshgrid of Gauss-Uniform points on the spheroid.
 * \param u Pointer to a gsl_matrix of Gaussian spaced rows
 * \param v Pointer to a gsl_matrix of Uniform spaced columns
 *
 * \p u is uniformly spaced on [0, pi] and \p v is spaced according to Gaussian quadrature.
 * `[u, v]` is structured according to the MATLAB `meshgrid` convention.
 */
void gl_grid( gsl_matrix * u, gsl_matrix * v )
{
    // Store references to sizes for clarity
    auto &nu = u->size1, &nv = u->size2;

    // Allocate temporary vectors
    gsl_vector* theta = gsl_vector_alloc( nu );
    gsl_vector* lambda = gsl_vector_alloc( nv );

    // Put Gaussian quadrature points in theta
    g_grid( nu, theta );

    for( int i = 0; i < nv; ++i )
        gsl_vector_set( lambda, i, 2 * M_PI * i / nv );

    fmt::print("Lambda: ");
    gsl_vector_print( lambda );
    
    fmt::print("Theta: ");
    gsl_vector_print( theta );
    
    // TODO: Make sure this is working
    for( int i = 0; i < nu; ++i )
        for( int j = 0; j < nv; ++j )
        {
            gsl_matrix_set( u, i, j, gsl_vector_get(theta, i) );
            gsl_matrix_set( v, i, j, gsl_vector_get(lambda, j) );
        }

    gsl_vector_free(lambda);
    gsl_vector_free(theta);
}