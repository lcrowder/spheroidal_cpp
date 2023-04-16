#include <yawg/core.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

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
gsl::cmatrix::cmatrix(const gsl_matrix_complex *gmat_other)
{
    this->calloc(gmat_other->size1, gmat_other->size2);
    gsl_matrix_complex_memcpy(gmat, gmat_other);
}

gsl::cmatrix::cmatrix(const gsl::matrix &gmat_other)
{
    this->calloc(gmat_other.nrows(), gmat_other.ncols());
    for (size_t i = 0; i < gmat_other.nrows(); i++)
        for (size_t j = 0; j < gmat_other.ncols(); j++)
            GSL_SET_COMPLEX(gsl_matrix_complex_ptr(gmat, i, j), gsl_matrix_get(gmat_other.gmat, i, j), 0.0);
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
void gsl::cmatrix::set(size_t i, size_t j, gsl::complex val)
{
    GSL_SET_COMPLEX(gsl_matrix_complex_ptr(gmat, i, j), val.real(), val.imag());
}

//! \brief Element getter
gsl::complex gsl::cmatrix::get(size_t i, size_t j) const
{
    return gsl::complex(gsl_matrix_complex_get(gmat, i, j));
}

gsl::complex_ref gsl::cmatrix::operator()(size_t i, size_t j)
{
    return gsl::complex_ref(gsl_matrix_complex_ptr(gmat, i, j));
}

const gsl::complex_ref gsl::cmatrix::operator()(size_t i, size_t j) const
{
    return gsl::complex_ref(gsl_matrix_complex_ptr(gmat, i, j));
}

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

//! \brief Reshape the array
void gsl::cmatrix::reshape(size_t n, size_t m)
{
    // if (n * m != this->size())
    //     throw std::runtime_error("Cannot reshape cmatrix to new size");
    gmat->size1 = n;
    gmat->size2 = m;
}

//! \brief Clear the gsl::matric and free the underlying gsl_matrix_complex
void gsl::cmatrix::clear()
{
    if (gmat == nullptr)
        return;
    this->free();
}

//! \brief Print the cmatrix to stdout
void gsl::cmatrix::print(FILE *out) const
{
    for (int i = 0; i < gmat->size1; ++i)
    {
        fprintf(out, (i == 0) ? "[" : " ");
        for (int j = 0; j < gmat->size2; ++j)
        {
            auto x = gsl_matrix_complex_get(gmat, i, j);
            fprintf(out, "%s% 9g%+9gj", ((j == 0) ? "" : ", "), GSL_REAL(x), GSL_IMAG(x));
        }
        fprintf(out, (i == (gmat->size1 - 1) ? "]\n" : ",\n"));
    }
}

namespace gsl
{
    cmatrix operator*(const cmatrix &A, const cmatrix &B)
    {
        cmatrix C(A.nrows(), B.ncols());
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, complex(1, 0), A.gmat, B.gmat, complex(0, 0), C.gmat);
        return C;
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

gsl::cmatrix_view &gsl::cmatrix_view::operator=(const gsl::cmatrix &m)
{
    gsl_matrix_complex_memcpy(&gmat_view.matrix, m.gmat);
    return *this;
}

gsl::cmatrix_view &gsl::cmatrix_view::operator=(gsl::cmatrix_view m)
{
    gsl_matrix_complex_memcpy(&gmat_view.matrix, &m.gmat_view.matrix);
    return *this;
}

gsl::cmatrix_view gsl::cmatrix::submatrix( size_t i, size_t j, size_t n, size_t m )
{
    return gsl::cmatrix_view( gsl_matrix_complex_submatrix( gmat, i, j, n, m ) );
}

gsl::crow_view gsl::cmatrix::row( size_t i )
{
    return gsl::crow_view( gsl_matrix_complex_row( gmat, i ) );
}

gsl::ccolumn_view gsl::cmatrix::column( size_t j )
{
    return gsl::ccolumn_view( gsl_matrix_complex_column( gmat, j ) );
}
