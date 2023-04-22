#include <yawg/core.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

//! \brief Construct empty matrix
gsl::cmatrix::cmatrix() : gmat(nullptr) {}

//! \brief Construct zero matrix of size n x m
gsl::cmatrix::cmatrix(size_t n, size_t m) : gmat(nullptr)
{
    if (n == 0 && m == 0)
        return;
    this->calloc(n, m);
}

//! \brief Construct new gsl::matrix from gsl_matrix
gsl::cmatrix::cmatrix(const gsl_matrix_complex *gmat_other)
{
    this->calloc(gmat_other->size1, gmat_other->size2);
    gsl_matrix_complex_memcpy(gmat, gmat_other);
}

//! \brief Construct new n x 1 gsl::cmatrix from a gsl::cvector
gsl::cmatrix::cmatrix(const gsl::cvector &v)
{
    this->calloc(v.size(), 1);
    for (size_t i = 0; i < v.size(); i++)
        gsl_matrix_complex_set(gmat, i, 0, gsl_vector_complex_get(v.gvec, i));
}

//! \brief Construct new gsl::cmatrix from MATLAB's .csv file format
gsl::cmatrix::cmatrix(FILE *in)
{
    // Count the number of rows and columns
    size_t n = 0;
    size_t m = 0;
    char c;
    while ((c = fgetc(in)) != EOF)
    {
        if (c == '\n')
            ++n;
        else if (c == ',')
            ++m;
    }
    m = (m / n) + 1;
    this->calloc(n, m);

    // Rewind the file
    rewind(in);

    // Read the data
    double real, imag;
    char sep;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            fscanf(in, "%lf%lfi,", &real, &imag);
            GSL_SET_COMPLEX(gsl_matrix_complex_ptr(gmat, i, j), real, imag);
        }
}

/*! \brief Copy constructor creating n x m matrix
 * \param M gsl::matrix to copy
 * \param n Number of rows
 * \param m Number of columns
 */
gsl::cmatrix::cmatrix(const gsl::cmatrix &M, size_t n, size_t m)
{
    this->calloc(n, m);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < m; j++)
            gsl_matrix_complex_set(gmat, i, j, gsl_matrix_complex_get(M.gmat, i, j));
}

/*! \brief Construct new gsl::cmatrix from gsl::matrix
 *
 *  Copy the values from a real matrix into a complex matrix, setting
 *  the imaginary part to zero.
 */
gsl::cmatrix::cmatrix(const gsl::matrix &M)
{
    this->calloc(M.nrows(), M.ncols());
    for (size_t i = 0; i < M.nrows(); i++)
        for (size_t j = 0; j < M.ncols(); j++)
            GSL_SET_COMPLEX(gsl_matrix_complex_ptr(gmat, i, j), gsl_matrix_get(M.gmat, i, j), 0.0);
}

gsl::cmatrix::cmatrix(const gsl::cmatrix &gmat_other)
{
    this->calloc(gmat_other.nrows(), gmat_other.ncols());
    gsl_matrix_complex_memcpy(gmat, gmat_other.gmat);
}

gsl::cmatrix::cmatrix(gsl::cmatrix &&gmat_other) : gmat(gmat_other.gmat)
{
    gmat_other.gmat = nullptr;
}

gsl::cmatrix &gsl::cmatrix::operator=(const gsl::cmatrix &gmat_other)
{
    if (this == &gmat_other)
        return *this;
    this->resize(gmat_other.nrows(), gmat_other.ncols());
    gsl_matrix_complex_memcpy(gmat, gmat_other.gmat);
    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator=(const gsl::matrix &M)
{
    this->resize(M.nrows(), M.ncols());
    for (size_t i = 0; i < M.nrows(); i++)
        for (size_t j = 0; j < M.ncols(); j++)
            GSL_SET_COMPLEX(gsl_matrix_complex_ptr(gmat, i, j), gsl_matrix_get(M.gmat, i, j), 0.0);
    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator=(gsl::cmatrix &&gmat_other)
{
    if (this == &gmat_other)
        return *this;
    this->free();
    gmat = gmat_other.gmat;
    gmat_other.gmat = nullptr;
    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator+=(const gsl::cmatrix &M)
{
    gsl_matrix_complex_add(gmat, M.gmat);
    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator+=(const gsl::matrix &M)
{
    // Add the real element of M to the complex element of this
    for (size_t i = 0; i < M.nrows(); i++)
        for (size_t j = 0; j < M.ncols(); j++)
            gsl_matrix_complex_set(gmat, i, j, gsl_complex_add_real(gsl_matrix_complex_get(gmat, i, j), gsl_matrix_get(M.gmat, i, j)));

    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator-=(const gsl::cmatrix &M)
{
    gsl_matrix_complex_sub(gmat, M.gmat);
    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator-=(const gsl::matrix &M)
{
    // Subtract the real element of M from the complex element of this
    for (size_t i = 0; i < M.nrows(); i++)
        for (size_t j = 0; j < M.ncols(); j++)
            gsl_matrix_complex_set(gmat, i, j, gsl_complex_sub_real(gsl_matrix_complex_get(gmat, i, j), gsl_matrix_get(M.gmat, i, j)));

    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator*=(complex z)
{
    gsl_matrix_complex_scale(gmat, z);
    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator*=(double x)
{
    // Multiply each element by x
    for (size_t i = 0; i < gmat->size1; i++)
        for (size_t j = 0; j < gmat->size2; j++)
            gsl_matrix_complex_set(gmat, i, j, gsl_complex_mul_real(gsl_matrix_complex_get(gmat, i, j), x));

    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator/=(complex z)
{
    gsl_matrix_complex_scale(gmat, 1.0 / z);
    return *this;
}

gsl::cmatrix &gsl::cmatrix::operator/=(double x)
{
    // Divide each element by x
    for (size_t i = 0; i < gmat->size1; i++)
        for (size_t j = 0; j < gmat->size2; j++)
            gsl_matrix_complex_set(gmat, i, j, gsl_complex_div_real(gsl_matrix_complex_get(gmat, i, j), x));

    return *this;
}

// gsl::cmatrix gsl::cmatrix::operator-() const
// {
//     gsl::cmatrix result(*this);
//     return -1.0 * result;
// }

gsl::cmatrix::~cmatrix()
{
    if (gmat == nullptr)
        return;
    gsl_matrix_complex_free(gmat);
}

void gsl::cmatrix::set(size_t i, size_t j, gsl::complex val)
{
    GSL_SET_COMPLEX(gsl_matrix_complex_ptr(gmat, i, j), val.real(), val.imag());
}

gsl::complex gsl::cmatrix::get(size_t i, size_t j) const
{
    return gsl::complex(gsl_matrix_complex_get(gmat, i, j));
}

/*! \brief Return a reference to the element at position (i,j)
 *
 *  This function returns a complex_ref to the element
 *  at position (i,j) in the matrix. Allows setting
 */
gsl::complex_ref gsl::cmatrix::operator()(size_t i, size_t j)
{
    return gsl::complex_ref(gsl_matrix_complex_ptr(gmat, i, j));
}

/*! \brief Return a const reference to the element at position (i,j)
 *
 *  This function returns a constant complex_ref to the element
 *  at position (i,j) in the matrix. Allows getting.
 */
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

/*! \brief Resize the gsl::cmatrix, setting elements to zero
 * \param n Number of rows
 * \param m Number of columns
 *
 * \note This function will always free and reallocate memory,
 *  setting the elements to zero.
 */
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

//! \brief CLear the gsl::cmatrix, free underlying memory
void gsl::cmatrix::clear()
{
    if (gmat == nullptr)
        return;
    this->free();
}

/*! \brief Return a new n x m gsl::cmatrix with same elements
 * \param n Number of rows
 * \param m Number of columns
 *
 * \return New gsl::matrix with same elements
 */
gsl::cmatrix gsl::cmatrix::reshape(size_t n, size_t m) const
{
    gsl::cmatrix M_new(n, m);
    for (size_t t = 0; t < n * m; t++)
        M_new(t / m, t % m) = this->get(t / this->ncols(), t % this->ncols());

    return M_new;
}

//! \brief Replace the complex matrix with its transpose
gsl::cmatrix &gsl::cmatrix::T()
{
    gsl_matrix_complex_transpose(gmat);
    return *this;
}

//! \brief Replace the complex matrix with its conjugate transpose
gsl::cmatrix &gsl::cmatrix::H()
{
    return this->T().conj();
}

//! \brief Replace the matrix with its conjugate
gsl::cmatrix &gsl::cmatrix::conj()
{
    // Replace each element with its conjugate
    for (size_t i = 0; i < this->nrows(); i++)
        for (size_t j = 0; j < this->ncols(); j++)
            gsl_matrix_complex_set(gmat, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(gmat, i, j)));
    return *this;
}

//! \brief Pretty-print the complex matrix to file stream
//! \param out File stream to print to
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
/*! \brief Print the complex matrix to file stream in MATLAB's .csv format
 * \param out File stream to print to
 *
 * \note This function uses a complex valued .csv format
 * compatible with MATLAB's load/save functions, which has the following format
 * 1.0000+2.0000i,2.0000+3.0000i,3.0000+4.0000i
 * 2.0000+3.0000i,4.0000+5.0000i,6.0000+7.0000i
 */
void gsl::cmatrix::print_csv(FILE *out) const
{
    for (int i = 0; i < gmat->size1; ++i)
    {
        for (int j = 0; j < gmat->size2; ++j)
        {
            auto x = gsl_matrix_complex_get(gmat, i, j);
            fprintf(out, "%lf%+lfi%c", GSL_REAL(x), GSL_IMAG(x), ((j == gmat->size2 - 1) ? '\n' : ','));
        }
    }
}

//! \brief Load the complex matrix from file stream in MATLAB's .csv format
//! \param in File stream to load from
void gsl::cmatrix::load_csv(FILE *in)
{
    *this = gsl::cmatrix(in);
}

namespace gsl
{
    cmatrix operator*(const cmatrix &A, const cmatrix &B)
    {
        cmatrix C(A.nrows(), B.ncols());
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, A.gmat, B.gmat, GSL_COMPLEX_ZERO, C.gmat);
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
 * \brief Private function to (continuously) allocate memory
 * \note This method allocates contiguous, zero-initialized memory.
 *       This is slightly slower than using gsl_matrix_compllex_alloc, but allows
 *       for intuitive usage of row views.
 */
void gsl::cmatrix::calloc(size_t n, size_t m) { gmat = gsl_matrix_complex_calloc(n, m); }