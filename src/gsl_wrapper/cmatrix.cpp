#include <gsl_wrapper/core.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <fmt/core.h>


//! \brief Default constructor
gsl::cmatrix::cmatrix() : gmat(nullptr) {}

//! \brief Constructor with size
gsl::cmatrix::cmatrix(size_t n, size_t m) : gmat(nullptr)
{
    if (n == 0 && m == 0)
        return;
    this->calloc(n, m);
}

//! \brief Build a gsl::cmatrix from a gsl_matrix_complex
gsl::cmatrix::cmatrix(gsl_matrix_complex *gmat_other)
{
    this->calloc(gmat_other->size1, gmat_other->size2);
    gsl_matrix_complex_memcpy(gmat, gmat_other);
}

//! \brief Copy constructor
gsl::cmatrix::cmatrix(const gsl::cmatrix &gmat_other)
{
    this->calloc(gmat_other.nrows(), gmat_other.ncols());
    gsl_matrix_complex_memcpy(gmat, gmat_other.gmat);
}

//! \brief Move constructor
gsl::cmatrix::cmatrix(gsl::cmatrix &&gmat_other) : gmat(gmat_other.gmat)
{
    gmat_other.gmat = nullptr;
}

//! \brief Assignment operator
gsl::cmatrix &gsl::cmatrix::operator=(const gsl::cmatrix &gmat_other)
{
    if (this == &gmat_other)
        return *this;
    this->resize(gmat_other.nrows(), gmat_other.ncols());
    gsl_matrix_complex_memcpy(gmat, gmat_other.gmat);
    return *this;
}

//! \brief Move assignment operator
gsl::cmatrix &gsl::cmatrix::operator=(gsl::cmatrix &&gmat_other)
{
    if (this == &gmat_other)
        return *this;
    this->free();
    gmat = gmat_other.gmat;
    gmat_other.gmat = nullptr;
    return *this;
}

//! \brief Destructor
gsl::cmatrix::~cmatrix()
{
    if (gmat == nullptr)
        return;
    gsl_matrix_complex_free(gmat);
}

//! \brief Element setter
void gsl::cmatrix::set(size_t i, size_t j, const gsl::complex &val)
{
    GSL_SET_COMPLEX( gsl_matrix_complex_ptr(gmat, i, j), val.real(), val.imag() );
}
//gsl_complex &gsl::cmatrix::operator()(size_t i, size_t j) { return *gsl_matrix_complex_ptr(gmat, i, j); }

//! \brief Element getter
gsl::complex gsl::cmatrix::get(size_t i, size_t j) const
{
    return gsl::complex( gsl_matrix_complex_get(gmat, i, j) );
}
//const gsl::complex &gsl::cmatrix::operator()(size_t i, size_t j) const { return *gsl_matrix_complex_ptr(gmat, i, j); }

//! \brief Size accessor
size_t gsl::cmatrix::size() const
{
    if (gmat == nullptr)
        return 0;
    return gmat->size1 * gmat->size2;
};

//! \brief Number of rows accessor
size_t gsl::cmatrix::nrows() const
{
    if (gmat == nullptr)
        return 0;
    return gmat->size1;
}

//! \brief Number of columns accessor
size_t gsl::cmatrix::ncols() const
{
    if (gmat == nullptr)
        return 0;
    return gmat->size2;
}

//! \brief Resize the gsl::cmatrix
void gsl::cmatrix::resize(size_t n, size_t m)
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

//! \brief Clear the gsl::matric and free the underlying gsl_matrix_complex
void gsl::cmatrix::clear()
{
    if (gmat == nullptr)
        return;
    this->free();
}

//! \brief Print the cmatrix to stdout
void gsl::cmatrix::print() const
{
    for (int i = 0; i < gmat->size1; ++i)
    {
        fmt::print((i == 0) ? "[" : " ");
        for (int j = 0; j < gmat->size2; ++j)
        {
            auto x = gsl_matrix_complex_get(gmat, i, j);
            fmt::print((j == 0) ? "{: 9g}{:+9g}j" : " {: 9g}{:+9g}j", GSL_REAL(x), GSL_IMAG(x));
        }
        fmt::print((i == (gmat->size1 - 1) ? "]\n" : "\n"));
    }
}

/*------ Protected Methods for gsl::cmatrix ------*/
//! \brief Free memory for underlying gsl_matrix_complex
void gsl::cmatrix::free()
{
    gsl_matrix_complex_free(gmat);
    gmat = nullptr;
}

/*!
 * \brief Allocate memory for underlying gsl_matrix_complex
 * \note This method allocates contiguous, zero-initialized memory.
 *       This is slightly slower than using gsl_matrix_alloc
 */
void gsl::cmatrix::calloc(size_t n, size_t m) { gmat = gsl_matrix_complex_calloc(n, m); }