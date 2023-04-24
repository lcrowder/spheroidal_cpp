#include <yawg/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>

//! \brief Construct empty matrix
gsl::matrix::matrix()
{
    gmat = new gsl_matrix;
    gmat->size1 = 0;
    gmat->size2 = 0;
    gmat->data = 0;
}

/*! \brief Construct zero matrix of size n x m
 * \param n Number of rows
 * \param m Number of columns
 *
 * \note By convention, all "empty" matrices have nullprt data
 */
gsl::matrix::matrix(size_t n, size_t m)
{
    this->galloc(n, m);
}

//! \brief Construct new gsl::vector from gsl_vector's data
gsl::matrix::matrix(gsl_matrix *gmat_other) : gmat(gmat_other) {}

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
    galloc(n, m);

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
    this->galloc(M.nrows(), M.ncols());
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
    galloc(v.size(), 1);
    for (size_t i = 0; i < v.size(); i++)
        set(i, 0, v(i));
}

gsl::matrix &gsl::matrix::operator=(const gsl::matrix &M)
{
    if (this == &M)
        return *this;
    if (nrows() != M.nrows() || ncols() != M.ncols())
        this->resize(M.nrows(), M.ncols());
    gsl_matrix_memcpy(gmat, M.gmat);
    return *this;
}

gsl::matrix &gsl::matrix::operator=(gsl::matrix &&M)
{
    if (this == &M)
        return *this;
    this->gfree();
    gmat = M.gmat;
    M.gmat = nullptr;
    return *this;
}

gsl::matrix &gsl::matrix::operator+=(const matrix &M)
{
    gsl_matrix_add(gmat, M.gmat);
    return *this;
}

gsl::matrix &gsl::matrix::operator-=(const matrix &M)
{
    gsl_matrix_sub(gmat, M.gmat);
    return *this;
}

gsl::matrix &gsl::matrix::operator*=(double x)
{
    gsl_matrix_scale(gmat, x);
    return *this;
}

gsl::matrix &gsl::matrix::operator/=(double x)
{
    gsl_matrix_scale(gmat, 1.0 / x);
    return *this;
}

gsl::matrix gsl::matrix::operator-() const
{
    gsl::matrix M(*this);
    return -1.0 * M;
}

gsl::matrix::~matrix()
{
    gfree();
}

double &gsl::matrix::operator()(size_t i, size_t j) { return *gsl_matrix_ptr(gmat, i, j); }
void gsl::matrix::set(size_t i, size_t j, double val) { *gsl_matrix_ptr(gmat, i, j) = val; }

void gsl::matrix::set_row(size_t i, const gsl::vector &v)
{
    gsl_matrix_set_row(gmat, i, v.get());
}

void gsl::matrix::set_col(size_t j, const gsl::vector &v)
{
    gsl_matrix_set_col(gmat, j, v.get());
}

double gsl::matrix::operator()(size_t i, size_t j) const { return *gsl_matrix_ptr(gmat, i, j); }
double gsl::matrix::get(size_t i, size_t j) const { return *gsl_matrix_ptr(gmat, i, j); }

gsl::vector gsl::matrix::get_row(size_t i) const
{
    gsl::vector v(gmat->size2);
    gsl_matrix_get_row(v.get(), gmat, i);
    return v;
}

gsl::vector gsl::matrix::get_col(size_t j) const
{
    gsl::vector v(gmat->size1);
    gsl_matrix_get_col(v.get(), gmat, j);
    return v;
}

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

bool gsl::matrix::is_square() const
{
    if (gmat == nullptr)
        return false;
    return (this->nrows() == this->ncols());
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
        clear();
        return;
    }
    // Don't free an empty vector
    if (gmat != nullptr)
    {
        if ((n == nrows()) && (m == ncols()))
            return;
        this->gfree();
    }
    galloc(n, m);
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
        gmat_new(t / n, t % m) = this->get(t / this->nrows(), t % this->ncols());

    return gmat_new;
}

//! \brief Compute the matrix transpose, in-place if square
gsl::matrix &gsl::matrix::T()
{
    if (this->is_square())
        gsl_matrix_transpose_memcpy(gmat, gmat);
    else
    {
        gsl::matrix M_new(this->nrows(), this->ncols());
        for (size_t i = 0; i < this->nrows(); i++)
            for (size_t j = 0; j < this->ncols(); j++)
                M_new(j, i) = this->get(i, j);
        *this = M_new;
    }
    return *this;
}

/*! \brief Copy constructor creating n x m matrix
 * \param M gsl::matrix to copy
 * \param n Number of rows
 * \param m Number of columns
 */
gsl::matrix::matrix(const matrix &M, size_t n, size_t m)
{
    galloc(M.gmat->size1, M.gmat->size2);

    for (size_t t = 0; t < n * m; t++)
        gsl_matrix_set(gmat, t / n, t % m,
                       gsl_matrix_get(M.gmat, t / M.gmat->size1, t % M.gmat->size2));
}

//! \brief CLear the gsl::matrix, free underlying memory
void gsl::matrix::clear()
{
    if (gmat == nullptr)
        return;
    this->gfree();
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
    *this = gsl::matrix(in);
}

namespace gsl
{
    matrix operator*(double a, const matrix &M)
    {
        matrix result(M);
        result *= a;
        return result;
    }

    matrix operator*(double a, matrix &&M)
    {
        M *= a;
        return M;
    }

    matrix operator*(const matrix &M, double a)
    {
        matrix result(M);
        result *= a;
        return result;
    }

    matrix operator*(matrix &&M, double a)
    {
        M *= a;
        return M;
    }

    matrix operator/(double a, const matrix &M)
    {
        matrix result(M);
        result /= a;
        return result;
    }

    matrix operator/(double a, matrix &&M)
    {
        M /= a;
        return M;
    }

    matrix operator/(const matrix &M, double a)
    {
        matrix result(M);
        result /= a;
        return result;
    }

    matrix operator/(matrix &&M, double a)
    {
        M /= a;
        return M;
    }

    matrix operator+(const matrix &M1, const matrix &M2)
    {
        matrix result(M1);
        result += M2;
        return result;
    }

    matrix operator+(const matrix &M1, matrix &&M2)
    {
        M2 += M1;
        return M2;
    }

    matrix operator+(matrix &&M1, const matrix &M2)
    {
        M1 += M2;
        return M1;
    }

    matrix operator+(matrix &&M1, matrix &&M2)
    {
        M1 += M2;
        return M1;
    }

    matrix operator-(const matrix &M1, const matrix &M2)
    {
        matrix result(M1);
        result -= M2;
        return result;
    }

    matrix operator-(const matrix &M1, matrix &&M2)
    {
        M2 -= M1;
        M2 *= -1.0;
        return M2;
    }

    matrix operator-(matrix &&M1, const matrix &M2)
    {
        M1 -= M2;
        return M1;
    }

    matrix operator-(matrix &&M1, matrix &&M2)
    {
        M1 -= M2;
        return M1;
    }

    matrix operator*(const matrix &A, const matrix &B)
    {
        matrix result(A.nrows(), B.ncols());
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A.gmat, B.gmat, 0.0, result.gmat);
        return result;
    }

    bool operator==(const matrix &M1, const matrix &M2)
    {
        return gsl_matrix_equal(M1.gmat, M2.gmat);
    }

    bool operator!=(const matrix &M1, const matrix &M2)
    {
        // Compare each element, returning true if any are not equal
        for (size_t i = 0; i < M1.nrows(); i++)
            for (size_t j = 0; j < M1.ncols(); j++)
                if (M1.get(i, j) != M2.get(i, j))
                    return true;
        return false;
    }
}

//! \brief Private function to free allocated memory
void gsl::matrix::gfree()
{
    if (gmat == nullptr)
        return;
    else if (size() == 0)
        delete gmat;
    else
        gsl_matrix_free(gmat);
    gmat = nullptr;
}

/*!
 * \brief Private function to (continuously) allocate memory
 * \note This method allocates contiguous, zero-initialized memory.
 *       This is slightly slower than using gsl_matrix_alloc, but allows
 *       for intuitive usage of row views.
 */
void gsl::matrix::galloc(size_t n, size_t m)
{
    if (n != 0 && m != 0)
        gmat = gsl_matrix_alloc(n, m);
    else
    {
        gmat = new gsl_matrix;
        gmat->size1 = 0;
        gmat->size2 = 0;
        gmat->data = 0;
    }
}

gsl::matrix_view gsl::matrix::view() const
{
    gsl_matrix *m = static_cast<gsl_matrix *>(malloc(sizeof(gsl_matrix)));
    *m = gsl_matrix_submatrix(get(), 0, 0, nrows(), ncols()).matrix;
    return gsl::matrix_view(m);
}

gsl::matrix_view gsl::matrix::submatrix(size_t i, size_t j, size_t n, size_t m) const
{
    gsl_matrix *v = static_cast<gsl_matrix *>(malloc(sizeof(gsl_matrix)));
    *v = gsl_matrix_submatrix(get(), i, j, n, m).matrix;
    return gsl::matrix_view(v);
}

gsl::row_view gsl::matrix::row(size_t i) const
{
    gsl_vector *w = static_cast<gsl_vector *>(malloc(sizeof(gsl_vector)));
    *w = gsl_matrix_row(gmat, i).vector;
    return gsl::row_view(w);
}

gsl::column_view gsl::matrix::column(size_t i) const
{
    gsl_vector *w = static_cast<gsl_vector *>(malloc(sizeof(gsl_vector)));
    *w = gsl_matrix_column(gmat, i).vector;
    return gsl::column_view(w);
}

gsl::matrix_view::matrix_view(gsl_matrix *gmat_other) : matrix(gmat_other) {}

gsl::matrix_view::~matrix_view()
{
    if (gmat != nullptr)
        free(gmat);
    gmat = nullptr;
}

void gsl::matrix_view::clear()
{
    printf("Warning: Attempting to clear a matrix view\n");
    gsl_matrix_set_zero(gmat);
}

void gsl::matrix_view::resize(size_t n, size_t m)
{
    printf("Warning: Attempting to resize a matrix view\n");
    gsl_matrix_set_zero(gmat);
}

gsl::matrix_view &gsl::matrix_view::operator=(const gsl::matrix &M)
{
    gsl_matrix_memcpy(gmat, M.get());
    return *this;
}