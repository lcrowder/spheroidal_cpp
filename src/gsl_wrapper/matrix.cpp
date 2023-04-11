#include <gsl_wrapper/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <fmt/core.h>

//! \brief Default constructor
gsl::matrix::matrix() : gmat(nullptr) {}

//! \brief Constructor with size
gsl::matrix::matrix(size_t n, size_t m) : gmat(nullptr)
{
    if (n == 0 && m == 0)
        return;
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