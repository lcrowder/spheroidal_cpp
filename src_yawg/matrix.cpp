#include <yawg/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>

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
gsl::matrix::matrix(const gsl_matrix *gmat_other)
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

//! \brief Convert vector into nx1 matrix
gsl::matrix::matrix(const gsl::vector &gvec_other)
{
    this->calloc(gvec_other.size(), 1);
    for (size_t i = 0; i < gvec_other.size(); i++)
        gsl_matrix_set(gmat, i, 0, gsl_vector_get(gvec_other.gvec, i));
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
void gsl::matrix::set(size_t i, size_t j, double val) { *gsl_matrix_ptr(gmat, i, j) = val; }

//! \brief Element getter
double gsl::matrix::operator()(size_t i, size_t j) const { return *gsl_matrix_ptr(gmat, i, j); }
double gsl::matrix::get(size_t i, size_t j) const { return *gsl_matrix_ptr(gmat, i, j); }

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

//! \brief Make a copy of the gsl::matrix by size n1 x m1 instead of n x m (note the number of elements is unchanges)
gsl::matrix gsl::matrix::reshape(size_t n, size_t m)
{
    gsl::matrix gmat_new(n, m);

    for (size_t t = 0; t < n * m; t++)
        gmat_new(t / m, t % m) = this->get(t / this->ncols(), t % this->ncols());

    return gmat_new;
}

gsl::matrix::matrix(const matrix &M, size_t n, size_t m)
{
    this->calloc(M.gmat->size1, M.gmat->size2);

    for (size_t t = 0; t < n * m; t++)
        gsl_matrix_set(gmat, t / m, t % m,
                       gsl_matrix_get(M.gmat, t / M.gmat->size1, t % M.gmat->size2));
}

//! \brief Clear the gsl::matric and free the underlying gsl_matrix
void gsl::matrix::clear()
{
    if (gmat == nullptr)
        return;
    this->free();
}

//! \brief Print the matrix to stdout
void gsl::matrix::print(FILE *out) const
{
    for (int i = 0; i < gmat->size1; ++i)
    {
        fprintf(out, (i == 0) ? "[" : " ");
        for (int j = 0; j < gmat->size2; ++j)
            fprintf(out, "%s% 9g", ((j == 0) ? "" : ", "), gsl_matrix_get(gmat, i, j));
        fprintf(out, (i == (gmat->size1 - 1) ? "]\n" : ",\n"));
    }
}

void gsl::matrix::print_csv(FILE *out) const
{
    for (int i = 0; i < gmat->size1; ++i)
    {
        for (int j = 0; j < gmat->size2; ++j)
            fprintf(out, "%.17g,", gsl_matrix_get(gmat, i, j));
        fprintf(out, "\n");
    }
}

void gsl::matrix::load_csv(FILE *in)
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
    this->resize(n, m);

    // Rewind the file
    rewind(in);

    // Read the data
    double val;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            fscanf(in, "%lf,", &val);
            gsl_matrix_set(gmat, i, j, val);
        }
}

namespace gsl
{
    matrix operator*(const matrix &A, const matrix &B)
    {
        matrix C(A.nrows(), B.ncols());
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A.gmat, B.gmat, 0.0, C.gmat);
        return C;
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

gsl::matrix_view &gsl::matrix_view::operator=(const gsl::matrix &m)
{
    gsl_matrix_memcpy(&gmat_view.matrix, m.gmat);
    return *this;
}

gsl::matrix_view &gsl::matrix_view::operator=(matrix_view m)
{
    gsl_matrix_memcpy(&gmat_view.matrix, &m.gmat_view.matrix);
    return *this;
}

void gsl::matrix_view::print(FILE *out) const
{
    for (int i = 0; i < gmat_view.matrix.size1; ++i)
    {
        fprintf(out, (i == 0) ? "[" : " ");
        for (int j = 0; j < gmat_view.matrix.size2; ++j)
            fprintf(out, "%s% 9g", ((j == 0) ? "" : ", "), gsl_matrix_get(&gmat_view.matrix, i, j));
        fprintf(out, (i == (gmat_view.matrix.size1 - 1) ? "]\n" : ",\n"));
    }
}

gsl::matrix_view gsl::matrix::submatrix( size_t i, size_t j, size_t n, size_t m )
{
    return gsl::matrix_view( gsl_matrix_submatrix( gmat, i, j, n, m ) );
}

gsl::row_view gsl::matrix::row( size_t i )
{
    return gsl::row_view( gsl_matrix_row( gmat, i ) );
}

gsl::column_view gsl::matrix::column( size_t j )
{
    return gsl::column_view( gsl_matrix_column( gmat, j ) );
}
