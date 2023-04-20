#include <yawg/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>

//! \brief Construct empty matrix
gsl::matrix::matrix() : gmat(nullptr) {}

/*! \brief Construct zero matrix of size n x m
 * \param n Number of rows
 * \param m Number of columns
 *
 * \note By convention, all "empty" matrices have nullprt data
 */
gsl::matrix::matrix(size_t n, size_t m) : gmat(nullptr)
{
    if (n == 0 && m == 0)
        return;
    this->calloc(n, m);
}

/*! \brief Construct new gsl::matrix from gsl_matrix
 *  \param gmat_other Pointer to existing gsl_matrix
 *
 *  \note This constructor does not copy the pointer to avoid
 *      double-freeing the memory. Instead, it copies the data.
 */
gsl::matrix::matrix(const gsl_matrix *gsl_mat)
{
    this->calloc(gsl_mat->size1, gsl_mat->size2);
    gsl_matrix_memcpy(gmat, gsl_mat);
}

//! \brief Construct new gsl::matrix from .csv file
//! \param in stdio.h file handle
gsl::matrix::matrix(FILE *in)
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
    double val;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            fscanf(in, "%lf,", &val);
            gsl_matrix_set(gmat, i, j, val);
        }
}

//! \brief Copy constructor
//! \param M gsl::matrix to copy
gsl::matrix::matrix(const gsl::matrix &M)
{
    this->calloc(M.nrows(), M.ncols());
    gsl_matrix_memcpy(gmat, M.gmat);
}

//! \brief Move constructor
//! \param M gsl::matrix to move
gsl::matrix::matrix(gsl::matrix &&M) : gmat(M.gmat)
{
    M.gmat = nullptr;
}

//! \brief Construct new n x 1 gsl::matrix from gsl::vector
//! \param v Vector to copy
gsl::matrix::matrix(const gsl::vector &v)
{
    this->calloc(v.size(), 1);
    for (size_t i = 0; i < v.size(); i++)
        gsl_matrix_set(gmat, i, 0, gsl_vector_get(v.gvec, i));
}

gsl::matrix &gsl::matrix::operator=(const gsl::matrix &gmat_other)
{
    if (this == &gmat_other)
        return *this;
    this->resize(gmat_other.nrows(), gmat_other.ncols());
    gsl_matrix_memcpy(gmat, gmat_other.gmat);
    return *this;
}

gsl::matrix &gsl::matrix::operator=(gsl::matrix &&gmat_other)
{
    if (this == &gmat_other)
        return *this;
    this->free();
    gmat = gmat_other.gmat;
    gmat_other.gmat = nullptr;
    return *this;
}

gsl::matrix::~matrix()
{
    if (gmat == nullptr)
        return;
    gsl_matrix_free(gmat);
}

double &gsl::matrix::operator()(size_t i, size_t j) { return *gsl_matrix_ptr(gmat, i, j); }
void gsl::matrix::set(size_t i, size_t j, double val) { *gsl_matrix_ptr(gmat, i, j) = val; }

double gsl::matrix::operator()(size_t i, size_t j) const { return *gsl_matrix_ptr(gmat, i, j); }
double gsl::matrix::get(size_t i, size_t j) const { return *gsl_matrix_ptr(gmat, i, j); }

size_t gsl::matrix::size() const
{
    if (gmat == nullptr)
        return 0;
    return gmat->size1 * gmat->size2;
};

size_t gsl::matrix::nrows() const
{
    if (gmat == nullptr)
        return 0;
    return gmat->size1;
}

size_t gsl::matrix::ncols() const
{
    if (gmat == nullptr)
        return 0;
    return gmat->size2;
}

/*! \brief Resize the gsl::matrix, setting elements to zero
 * \param n Number of rows
 * \param m Number of columns
 *
 * \note This function will always free and reallocate memory,
 *  setting the elements to zero.
 */
void gsl::matrix::resize(size_t n, size_t m)
{
    if ((n == 0) || (m == 0))
    {
        this->clear();
        return;
    }
    // Don't free an empty vector
    if (gmat != nullptr)
        this->free();
    this->calloc(n, m);
}

/*! \brief Return a new n x m gsl::matrix with same elements
 * \param n Number of rows
 * \param m Number of columns
 *
 * \return New gsl::matrix with same elements
 */
gsl::matrix gsl::matrix::reshape(size_t n, size_t m) const
{
    gsl::matrix gmat_new(n, m);

    for (size_t t = 0; t < n * m; t++)
        gmat_new(t / m, t % m) = this->get(t / this->ncols(), t % this->ncols());

    return gmat_new;
}

gsl::matrix & gsl::matrix::T()
{
    gsl_matrix_transpose(gmat);
    return *this;
}

/*! \brief Copy constructor creating n x m matrix
 * \param M gsl::matrix to copy
 * \param n Number of rows
 * \param m Number of columns
 */
gsl::matrix::matrix(const matrix &M, size_t n, size_t m)
{
    this->calloc(M.gmat->size1, M.gmat->size2);

    for (size_t t = 0; t < n * m; t++)
        gsl_matrix_set(gmat, t / m, t % m,
                       gsl_matrix_get(M.gmat, t / M.gmat->size1, t % M.gmat->size2));
}

//! \brief CLear the gsl::matrix, free underlying memory
void gsl::matrix::clear()
{
    if (gmat == nullptr)
        return;
    this->free();
}

//! \brief Pretty-print the matrix to file stream
//! \param out File stream to print to
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

//! \brief Print the matrix to file stream in CSV format
//! \param out File stream to print to
void gsl::matrix::print_csv(FILE *out) const
{
    for (int i = 0; i < gmat->size1; ++i)
        for (int j = 0; j < gmat->size2; ++j)
            fprintf(out, "%.17g%c", gsl_matrix_get(gmat, i, j), ((j == gmat->size2 - 1) ? '\n' : ','));
}

//! \brief Load the matrix from a file stream in CSV format
//! \param in File stream to load from
void gsl::matrix::load_csv(FILE *in)
{
    *this = gsl::matrix( in );
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

//! \brief Private function to free allocated memory
void gsl::matrix::free()
{
    gsl_matrix_free(gmat);
    gmat = nullptr;
}

/*!
 * \brief Private function to (continuously) allocate memory
 * \note This method allocates contiguous, zero-initialized memory.
 *       This is slightly slower than using gsl_matrix_alloc, but allows
 *       for intuitive usage of row views.
 */
void gsl::matrix::calloc(size_t n, size_t m) { gmat = gsl_matrix_calloc(n, m); }

//! \brief Assignment to a matrix view from a matrix
gsl::matrix_view &gsl::matrix_view::operator=(const gsl::matrix &m)
{
    gsl_matrix_memcpy(&gmat_view.matrix, m.gmat);
    return *this;
}

//! \brief Assignment to a matrix view from another matrix view
gsl::matrix_view &gsl::matrix_view::operator=(matrix_view Mv)
{
    gsl_matrix_memcpy(&gmat_view.matrix, &Mv.gmat_view.matrix);
    return *this;
}

//! \brief Pretty-print the matrix to file stream
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

/*! \brief Return a view to a submatrix of the matrix
 * \param i Starting row index
 * \param j Starting column index
 * \param n Number of rows
 * \param m Number of columns
 *
 * \return Matrix view to submatrix
 */
gsl::matrix_view gsl::matrix::submatrix(size_t i, size_t j, size_t n, size_t m)
{
    return gsl::matrix_view(gsl_matrix_submatrix(gmat, i, j, n, m));
}

/*! \brief Return a view to a row of the matrix
 *  \param i Row index
 *  \return Row view
 */
gsl::row_view gsl::matrix::row(size_t i)
{
    return gsl::row_view(gsl_matrix_row(gmat, i));
}

/*! \brief Return a view to a column of the matrix
 *  \param j Column index
 *  \return Column view
 */
gsl::column_view gsl::matrix::column(size_t j)
{
    return gsl::column_view(gsl_matrix_column(gmat, j));
}
