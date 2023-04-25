#include <yawg/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <utility>
#include <stdlib.h>

//! \brief Construct empty vector
gsl::vector::vector()
{
    gvec = new gsl_vector;
    gvec->size = 0;
    gvec->data = 0;
}

//! \brief Construct zero vector of size n
gsl::vector::vector(size_t n)
{
    this->galloc(n);
}

//! \brief Construct new gsl::vector from gsl_vector's data
gsl::vector::vector(gsl_vector *gvec_other) : gvec(gvec_other) {}

gsl::vector::vector(const gsl::vector &gvec_other)
{
    this->galloc(gvec_other.size());
    gsl_vector_memcpy(gvec, gvec_other.gvec);
}

gsl::vector::vector(gsl::vector &&v) : gvec(v.gvec)
{
    v.gvec = nullptr;
}

gsl::vector &gsl::vector::operator=(const gsl::vector &v)
{
    if (this == &v)
        return *this;
    if (size() != v.size())
        this->resize(v.size());
    gsl_vector_memcpy(gvec, v.gvec);
    return *this;
}

gsl::vector &gsl::vector::operator=(gsl::vector &&v)
{
    if (this == &v)
        return *this;
    this->gfree();
    gvec = v.gvec;
    v.gvec = nullptr;
    return *this;
}

gsl::vector &gsl::vector::operator+=(const gsl::vector &v)
{
    gsl_vector_add(gvec, v.gvec);
    return *this;
}

gsl::vector &gsl::vector::operator-=(const gsl::vector &gvec_other)
{
    gsl_vector_sub(gvec, gvec_other.gvec);
    return *this;
}

gsl::vector &gsl::vector::operator*=(double a)
{
    gsl_vector_scale(gvec, a);
    return *this;
}

gsl::vector &gsl::vector::operator/=(double a)
{
    gsl_vector_scale(gvec, 1.0 / a);
    return *this;
}

gsl::vector gsl::vector::operator-() const
{
    gsl::vector result(*this);
    return -1.0 * result;
}

gsl::vector::~vector()
{
    this->gfree();
}

void gsl::vector::set(size_t i, double val)
{
    *gsl_vector_ptr(gvec, i) = val;
}

double gsl::vector::get(size_t i) const
{
    return *gsl_vector_ptr(gvec, i);
}

double gsl::vector::operator()(size_t i) const
{
    return *gsl_vector_ptr(gvec, i);
}

double &gsl::vector::operator()(size_t i)
{
    return *gsl_vector_ptr(gvec, i);
}

size_t gsl::vector::size() const
{
    if (gvec == nullptr)
        return 0;
    return gvec->size;
};

/*! \brief Resize the gsl::vector, setting elements to zero
 * \param n Number of elements
 *
 * \note This function will always free and reallocate memory,
 *  setting the elements to zero.
 */
void gsl::vector::resize(size_t n)
{
    if (n == 0)
    {
        clear();
        return;
    }
    if (gvec != nullptr)
    {
        if (gvec->size == n)
            return;
        gfree();
    }
    galloc(n);
}

//! \brief Clear the gsl::vector, free underlying memory
void gsl::vector::clear()
{
    if (gvec == nullptr)
        return;
    this->gfree();
}

//! \brief Pretty-print the vector to file stream
//! \param out File stream to print to
void gsl::vector::print(FILE *out) const
{
    fprintf(out, "[");
    for (int i = 0; i < gvec->size; ++i)
        fprintf(out, "%s%g", (i == 0) ? "" : ", ", gsl_vector_get(gvec, i));
    fprintf(out, "]\n");
}

//! \brief Private function to free allocated memory
void gsl::vector::gfree()
{
    if (gvec == nullptr)
        return;
    else if (size() == 0)
        delete gvec;
    else
        gsl_vector_free(gvec);
    gvec = nullptr;
}

/*!
 * \brief Private function to (continuously) allocate memory
 * \note This method allocates contiguous, zero-initialized memory.
 *       This is slightly slower than using gsl_vector_alloc, but allows
 *       for intuitive usage of row views.
 */
void gsl::vector::galloc(size_t n)
{
    if (n >= 0)
        gvec = gsl_vector_alloc(n);
    else
    {
        gvec = new gsl_vector;
        gvec->size = 0;
        gvec->data = 0;
    }
}
namespace gsl
{
    vector operator*(double a, const vector &v)
    {
        vector result(v);
        result *= a;
        return result;
    }

    vector operator*(double a, vector &&v)
    {
        v *= a;
        return std::move(v);
    }

    vector operator*(const vector &v, double a)
    {
        vector result(v);
        result *= a;
        return result;
    }

    vector operator*(vector &&v, double a)
    {
        v *= a;
        return std::move(v);
    }

    vector operator/(double a, const vector &v)
    {
        vector result(v);
        result /= a;
        return result;
    }

    vector operator/(double a, vector &&v)
    {
        v /= a;
        return std::move(v);
    }

    vector operator/(const vector &v, double a)
    {
        vector result(v);
        result /= a;
        return result;
    }

    vector operator/(vector &&v, double a)
    {
        v /= a;
        return std::move(v);
    }

    vector operator+(const vector &v1, const vector &v2)
    {
        vector result(v1);
        result += v2;
        return result;
    }

    vector operator+(vector &&v1, const vector &v2)
    {
        v1 += v2;
        return std::move(v1);
    }

    vector operator+(const vector &v1, vector &&v2)
    {
        v2 += v1;
        return std::move(v2);
    }

    vector operator+(vector &&v1, vector &&v2)
    {
        v1 += v2;
        return std::move(v1);
    }

    cvector operator+(const vector &v1, const cvector &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator+(const vector &v1, cvector &&v2)
    {
        v2 += v1;
        return std::move(v2);
    }

    cvector operator+(const cvector &v1, const vector &v2)
    {
        cvector result(v2);
        result += v1;
        return result;
    }

    cvector operator+(cvector &&v1, const vector &v2)
    {
        v1 += v2;
        return std::move(v1);
    }

    vector operator-(const vector &v1, const vector &v2)
    {
        vector result(v1);
        result -= v2;
        return result;
    }

    vector operator-(vector &&v1, const vector &v2)
    {
        v1 -= v2;
        return std::move(v1);
    }

    vector operator-(const vector &v1, vector &&v2)
    {
        v2 -= v1;
        return std::move(v2);
    }

    vector operator-(vector &&v1, vector &&v2)
    {
        v1 -= v2;
        return std::move(v1);
    }

    // Compare vectors to vectors
    bool operator==(const vector &v1, const vector &v2)
    {
        return gsl_vector_equal(v1.gvec, v2.gvec);
    }

    bool operator!=(const vector &v1, const vector &v2)
    {
        // Compare each element, returning true if any are not equal
        for (size_t i = 0; i < v1.gvec->size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    vector operator*(const matrix& M, const vector& v)
    {
        vector result(M.nrows());
        gsl_blas_dgemv(CblasNoTrans, 1.0, M.get(), v.get(), 0.0, result.get());
        return result;
    }

}

gsl::vector_view gsl::vector::view() const
{
    gsl_vector *v = static_cast<gsl_vector *>(malloc(sizeof(gsl_vector)));
    *v = gsl_vector_subvector(get(), 0, size()).vector;
    return gsl::vector_view(v);
}

gsl::matrix_view gsl::row_view::reshape(size_t n, size_t m) const
{
    gsl_matrix *t = static_cast<gsl_matrix *>(malloc(sizeof(gsl_matrix)));
    *t = gsl_matrix_view_vector(get(), n, m).matrix;
    return gsl::matrix_view(t);
}

gsl::vector_view gsl::vector::subvector(size_t offset, size_t n) const
{
    gsl_vector *v = static_cast<gsl_vector *>(malloc(sizeof(gsl_vector)));
    *v = gsl_vector_subvector(get(), offset, n).vector;
    return gsl::vector_view(v);
}

//! \brief Construct new gsl::vector from gsl_vector
gsl::vector_view::vector_view(gsl_vector *gvec_other) : vector(gvec_other) {}

// Vector views should never own their memory, so we don't need to free it
// We only need to delete the pointer, since it was created on the heap
gsl::vector_view::~vector_view()
{
    if (gvec != nullptr)
        free(gvec);
    gvec = nullptr;
}

gsl::vector_view &gsl::vector_view::operator=(const gsl::vector &v)
{
    gsl_vector_memcpy(gvec, v.get());
    return *this;
}

//! Set all values in the view to zero
void gsl::vector_view::clear()
{
    printf("Warning: Attempting to clear a vector view\n");
    gsl_vector_set_zero(gvec);
}

//! Set all values in the view to zero
void gsl::vector_view::resize(size_t n)
{
    printf("Warning: Attempting to resize a vector view\n");
    gsl_vector_set_zero(gvec);
}

gsl::column_view &gsl::column_view::operator=(const gsl::vector &v)
{
    gsl_vector_memcpy(gvec, v.get());
    return *this;
}

gsl::row_view &gsl::row_view::operator=(const gsl::vector &v)
{
    gsl_vector_memcpy(gvec, v.get());
    return *this;
}