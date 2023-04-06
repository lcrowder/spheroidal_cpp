#include <spheroidal/gsl_wrapper.h>
#include <gsl/gsl_integration.h>

/*!
 * \fn \brief Computes \p n Gauss-Legendre quadrature nodes and weights for the interval [a, b] (default [-1, 1])
 * \param n Number of nodes
 * \param x Reference to a gsl::vector of nodes
 * \param w Reference to a gsl::vector of weights
 * \param a Lower bound of the interval (Default -1)
 * \param b Upper bound of the interval (Default +1)
 *  
 * \note This could be made more efficient by caching the gsl_integration_glfixed_table,
 * rather than allocating and freeing it each time.
 */
void leggauss( size_t n, gsl::vector& x, gsl::vector& w, double a, double b )
{
    x.resize(n);
    w.resize(n);

    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(n);
    for( int i = 0; i < n; ++i )
        gsl_integration_glfixed_point( a, b, i, &x(i), &w(i), t );
    
    gsl_integration_glfixed_table_free(t);
}

/*!
 * \fn \brief Computes \p n Gauss-Legendre quadrature nodes and weights for the interval [a, b] (default [-1, 1])
 * \param n Number of nodes
 * \param a Lower bound of the interval (Default -1)
 * \param b Upper bound of the interval (Default +1)
 
 * \return gsl::vector of nodes
 */
gsl::vector leggauss( size_t n, double a, double b )
{
    gsl::vector x(n);
    double trash_w;

    gsl_integration_glfixed_table* t = gsl_integration_glfixed_table_alloc(n);
    for( int i = 0; i < n; ++i )
        gsl_integration_glfixed_point( a, b, i, &x(i), &trash_w, t );
    
    gsl_integration_glfixed_table_free(t);

    return x;
}

/*!
 * \fn \brief Computes \p n evenly spaced points on the interval [\p a, \p b] (inclusive)
 * \param a Lower bound of the interval
 * \param b Upper bound of the interval
 * \param n Number of points
 *  
 * \return gsl::vector of points
 */
gsl::vector linspace( double a, double b, size_t n )
{
    gsl::vector x(n);
    double dx = (b-a)/(n-1);
    for( int i = 0; i < n; ++i )
        x(i) = a + i*dx;
    return x;
}

/*! Do matlab thing */
void meshgrid( const gsl::vector& x, const gsl::vector& y, gsl::matrix& X, gsl::matrix& Y )
{
    X.resize( x.size(), y.size() );
    Y.resize( x.size(), y.size() );

    for( int i = 0; i < x.size(); ++i )
        for( int j = 0; j < y.size(); ++j )
        {
            X(i,j) = x(i);
            Y(i,j) = y(j);
        }
}

/*! Compute arccos of each element in array */
gsl::vector acos( const gsl::vector& x )
{
    gsl::vector y(x.size());
    for( int i = 0; i < x.size(); ++i )
        y(i) = std::acos(x(i));
    return y;
}