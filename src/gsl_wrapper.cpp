#include <spheroidal/gsl_wrapper.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <fmt/core.h>

/*------ Public Methods for gsl::vector ------*/
gsl::vector::vector() : gvec( nullptr ) {}
gsl::vector::vector( size_t n ) { this->calloc(n); }
    
gsl::vector::vector( gsl_vector * gvec_other )
{
    this->calloc( gvec_other->size );
    gsl_vector_memcpy( gvec, gvec_other );
}

gsl::vector::vector( const gsl::vector& gvec_other )
{
    this->calloc( gvec_other.size() );
    gsl_vector_memcpy( gvec, gvec_other.gvec );
}

gsl::vector::~vector() 
{  
    if( gvec == nullptr )
        return;
    gsl_vector_free( gvec ); 
}
    
double& gsl::vector::operator()( size_t i ) { return *gsl_vector_ptr( gvec, i ); }

const double& gsl::vector::operator()( size_t i ) const { return *gsl_vector_ptr( gvec, i ); }

size_t gsl::vector::size() const { return gvec->size; };

void gsl::vector::resize( size_t n )
{
    // Don't free an empty vector
    if( gvec != nullptr ) 
    {
        // Only allocate new memory if size is different
        if( gvec->size == n ) return;
        this->free();
    }
    this->calloc( n );
}

void gsl::vector::print() const 
{
    fmt::print( "[" );
    for( int i = 0; i < gvec->size; ++i )
        fmt::print( (i == 0) ? "{:g}" : ", {:g}", gsl_vector_get(gvec, i) );
    fmt::print( "]\n" );
}

/*------ Protected Methods for gsl::vector ------*/
void gsl::vector::free()
{
    if( gvec != nullptr ) return;
    gsl_vector_free( gvec );
    gvec = nullptr;
}

void gsl::vector::calloc( size_t n )
{
    if( gvec != nullptr ) this->free();
    gvec = gsl_vector_calloc( n );
}


/*------ Public Methods for gsl::matrix ------*/
gsl::matrix::matrix() : gmat( nullptr ) {}
gsl::matrix::matrix( size_t n, size_t m ) { this->calloc(n, m); }
    
gsl::matrix::matrix( gsl_matrix * gmat_other )
{
    this->calloc( gmat_other->size1, gmat_other->size2 );
    gsl_matrix_memcpy( gmat, gmat_other );
}

gsl::matrix::matrix( const gsl::matrix& gmat_other )
{
    this->calloc( gmat_other.nrows(), gmat_other.ncols() );
    gsl_matrix_memcpy( gmat, gmat_other.gmat );
}

gsl::matrix::~matrix() 
{  
    if( gmat == nullptr )
        return;
    gsl_matrix_free( gmat ); 
}
    

double& gsl::matrix::operator()( size_t i, size_t j ) { return *gsl_matrix_ptr( gmat, i, j ); }

const double& gsl::matrix::operator()( size_t i, size_t j ) const { return *gsl_matrix_ptr( gmat, i, j ); }

size_t gsl::matrix::size() const { return gmat->size1 * gmat->size2; };
size_t gsl::matrix::nrows() const { return gmat->size1; };
size_t gsl::matrix::ncols() const { return gmat->size2; };

void gsl::matrix::resize( size_t n, size_t m )
{
    // Don't free an empty vector
    if( gmat != nullptr ) 
    {
        // Only allocate new memory if size is different
        if( gmat->size1 == n && gmat->size2 == m ) return;
        this->free();
    }
    this->calloc( n, m );
}

void gsl::matrix::print() const 
{
    for( int i = 0; i < gmat->size1; ++i )
    {
        fmt::print( (i == 0) ? "[" : " " );
        for( int j = 0; j < gmat->size2; ++j )
            fmt::print( (j == 0) ? "{: 9g}" : " {: 9g}", gsl_matrix_get(gmat, i, j) );
        fmt::print( (i == (gmat->size1 - 1) ? "]\n" : "\n") );
    }
}

/*------ Protected Methods for gsl::matrix ------*/
void gsl::matrix::free()
{
    if( gmat != nullptr ) return;
    gsl_matrix_free( gmat );
    gmat = nullptr;
}

void gsl::matrix::calloc( size_t n, size_t m )
{
    if( gmat != nullptr ) this->free();
    gmat = gsl_matrix_calloc( n, m );
}
