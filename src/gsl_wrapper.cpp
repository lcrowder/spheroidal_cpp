#include <spheroidal/gsl_wrapper.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <fmt/core.h>

/*------ Public Methods for gsl::vector ------*/
//! \brief Default constructor
gsl::vector::vector() : gvec(nullptr) {}

//! \brief Constructor with size
gsl::vector::vector(size_t n) : gvec(nullptr) { if ( n != 0 ) this->calloc(n); }

//! \brief Build a gsl::vector from a gsl_vector
gsl::vector::vector(gsl_vector *gvec_other)
{
    if( gvec_other == nullptr ) return;
    this->calloc(gvec_other->size);
    gsl_vector_memcpy(gvec, gvec_other);
}

//! \brief Copy constructor
gsl::vector::vector(const gsl::vector &gvec_other)
{
    this->calloc(gvec_other.size());
    gsl_vector_memcpy(gvec, gvec_other.gvec);
}

//! \brief Move constructor
gsl::vector::vector(gsl::vector &&gvec_other) : gvec(gvec_other.gvec)
{
    gvec_other.gvec = nullptr;
}

//! \brief Assignment operator
gsl::vector &gsl::vector::operator=(const gsl::vector &gvec_other)
{
    if (this == &gvec_other)
        return *this;
    this->resize(gvec_other.size());
    gsl_vector_memcpy(gvec, gvec_other.gvec);
    return *this;
}

//! \brief Move assignment operator
gsl::vector &gsl::vector::operator=(gsl::vector &&gvec_other)
{
    if (this == &gvec_other)
        return *this;
    this->free();
    gvec = gvec_other.gvec;
    gvec_other.gvec = nullptr;
    return *this;
}

//! \brief Destructor
gsl::vector::~vector()
{
    if (gvec == nullptr)
        return;
    gsl_vector_free(gvec);
}

//! \brief Element getter
const double &gsl::vector::operator()(size_t i) const { return *gsl_vector_ptr(gvec, i); }

//! \brief Element setter
double &gsl::vector::operator()(size_t i) { return *gsl_vector_ptr(gvec, i); }

//! \brief Size accessor
size_t gsl::vector::size() const
{
    if (gvec == nullptr)
        return 0;
    return gvec->size;
};

//! \brief Resize the gsl::vector
//! \note If n == 0, the vector is cleared
void gsl::vector::resize(size_t n)
{
    if (n == 0)
    {
        this->clear();
        return;
    }
    if (gvec != nullptr)
    {
        if (gvec->size == n)
            return;
        this->free();
    }
    this->calloc(n);
}

//! \brief Clear the vector and free the underlying gsl_vector
void gsl::vector::clear()
{
    if (gvec == nullptr)
        return;
    this->free();
}

//! \brief Print the vector to stdout
void gsl::vector::print() const
{
    fmt::print("[");
    for (int i = 0; i < gvec->size; ++i)
        fmt::print((i == 0) ? "{:g}" : ", {:g}", gsl_vector_get(gvec, i));
    fmt::print("]\n");
}

/*------ Protected Methods for gsl::vector ------*/

//! \brief Free memory for underlying gsl_vector
void gsl::vector::free()
{
    gsl_vector_free(gvec);
    gvec = nullptr;
}

/*!
 * \brief Allocate memory for underlying gsl_vector
 * \note This method allocates contiguous, zero-initialized memory.
 *         This is slightly slower than using gsl_vector_alloc
 */
void gsl::vector::calloc(size_t n) { gvec = gsl_vector_calloc(n); }

/*------ Public Methods for gsl::matrix ------*/
//! \brief Default constructor
gsl::matrix::matrix() : gmat(nullptr) {}

//! \brief Constructor with size
gsl::matrix::matrix(size_t n, size_t m) : gmat(nullptr)
{ 
    if ( n == 0 && m == 0 ) return;
    this->calloc(n, m); 
}

//! \brief Build a gsl::matrix from a gsl_matrix
gsl::matrix::matrix(gsl_matrix *gmat_other)
{
    this->calloc(gmat_other->size1, gmat_other->size2);
    gsl_matrix_memcpy(gmat, gmat_other);
}

//! \brief Copy constructor
gsl::matrix::matrix(const gsl::matrix &gmat_other)
{
    this->calloc(gmat_other.nrows(), gmat_other.ncols());
    gsl_matrix_memcpy(gmat, gmat_other.gmat);
}

//! \brief Move constructor
gsl::matrix::matrix(gsl::matrix &&gmat_other) : gmat(gmat_other.gmat)
{
    gmat_other.gmat = nullptr;
}

//! \brief Assignment operator
gsl::matrix &gsl::matrix::operator=(const gsl::matrix &gmat_other)
{
    if (this == &gmat_other)
        return *this;
    this->resize(gmat_other.nrows(), gmat_other.ncols());
    gsl_matrix_memcpy(gmat, gmat_other.gmat);
    return *this;
}

//! \brief Move assignment operator
gsl::matrix &gsl::matrix::operator=(gsl::matrix &&gmat_other)
{
    if (this == &gmat_other)
        return *this;
    this->free();
    gmat = gmat_other.gmat;
    gmat_other.gmat = nullptr;
    return *this;
}

//! \brief Destructor
gsl::matrix::~matrix()
{
    if (gmat == nullptr)
        return;
    gsl_matrix_free(gmat);
}

//! \brief Element setter
double &gsl::matrix::operator()(size_t i, size_t j) { return *gsl_matrix_ptr(gmat, i, j); }

//! \brief Element getter
const double &gsl::matrix::operator()(size_t i, size_t j) const { return *gsl_matrix_ptr(gmat, i, j); }

//! \brief Size accessor
size_t gsl::matrix::size() const
{
    if (gmat == nullptr)
        return 0;
    return gmat->size1 * gmat->size2;
};

//! \brief Number of rows accessor
size_t gsl::matrix::nrows() const
{
    if (gmat == nullptr)
        return 0;
    return gmat->size1;
}

//! \brief Number of columns accessor
size_t gsl::matrix::ncols() const
{
    if (gmat == nullptr)
        return 0;
    return gmat->size2;
}

//! \brief Resize the gsl::matrix
void gsl::matrix::resize(size_t n, size_t m)
{
    if ((n == 0) || (m == 0))
    {
        this->clear();
        return;
    }
    // Don't free an empty vector
    if (gmat != nullptr)
    {
        // Only allocate new memory if size is different
        if (gmat->size1 == n && gmat->size2 == m)
            return;
        this->free();
    }
    this->calloc(n, m);
}

//! \brief Clear the gsl::matric and free the underlying gsl_matrix
void gsl::matrix::clear()
{
    if (gmat == nullptr)
        return;
    this->free();
}

//! \brief Print the matrix to stdout
void gsl::matrix::print() const
{
    for (int i = 0; i < gmat->size1; ++i)
    {
        fmt::print((i == 0) ? "[" : " ");
        for (int j = 0; j < gmat->size2; ++j)
            fmt::print((j == 0) ? "{: 9g}" : " {: 9g}", gsl_matrix_get(gmat, i, j));
        fmt::print((i == (gmat->size1 - 1) ? "]\n" : "\n"));
    }
}

/*------ Protected Methods for gsl::matrix ------*/
//! \brief Free memory for underlying gsl_matrix
void gsl::matrix::free()
{
    gsl_matrix_free(gmat);
    gmat = nullptr;
}

/*!
 * \brief Allocate memory for underlying gsl_matrix
 * \note This method allocates contiguous, zero-initialized memory.
 *       This is slightly slower than using gsl_matrix_alloc
 */
void gsl::matrix::calloc(size_t n, size_t m) { gmat = gsl_matrix_calloc(n, m); }