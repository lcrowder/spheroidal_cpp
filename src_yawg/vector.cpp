#include <yawg/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <utility>

//! \brief Construct empty vector
gsl::vector::vector() : gvec(nullptr) {}

//! \brief Construct zero vector of size n
gsl::vector::vector(size_t n) : gvec(nullptr)
{
    if (n != 0)
        this->calloc(n);
}

//! \brief Construct new gsl::vector from gsl_vector
gsl::vector::vector(const gsl_vector *gvec_other)
{
    if (gvec_other == nullptr)
        return;
    this->calloc(gvec_other->size);
    gsl_vector_memcpy(gvec, gvec_other);
}

gsl::vector::vector(const gsl::vector &gvec_other)
{
    this->calloc(gvec_other.size());
    gsl_vector_memcpy(gvec, gvec_other.gvec);
}

gsl::vector::vector(gsl::vector &&gvec_other) : gvec(gvec_other.gvec)
{
    gvec_other.gvec = nullptr;
}

gsl::vector &gsl::vector::operator=(const gsl::vector &gvec_other)
{
    if (gvec == gvec_other.gvec)
        return *this;
    this->resize(gvec_other.size());
    gsl_vector_memcpy(gvec, gvec_other.gvec);
    return *this;
}

gsl::vector &gsl::vector::operator=(const gsl::vector_view &vv)
{
    if (gvec == &vv.gvec_view.vector)
        return *this;
    this->resize(vv.size());
    gsl_vector_memcpy(gvec, &vv.gvec_view.vector);
    return *this;
}

gsl::vector &gsl::vector::operator=(gsl::vector &&v)
{
    if (gvec == v.gvec)
        return *this;
    this->free();
    gvec = v.gvec;
    v.gvec = nullptr;
    return *this;
}

gsl::vector &gsl::vector::operator+=(const gsl::vector &v)
{
    gsl_vector_add(gvec, v.gvec);
    return *this;
}

gsl::vector &gsl::vector::operator+=(const gsl::vector_view &vv)
{
    gsl_vector_add(gvec, &vv.gvec_view.vector);
    return *this;
}

gsl::vector &gsl::vector::operator-=(const gsl::vector &gvec_other)
{
    gsl_vector_sub(gvec, gvec_other.gvec);
    return *this;
}

gsl::vector &gsl::vector::operator-=(const gsl::vector_view &vv)
{
    gsl_vector_sub(gvec, &vv.gvec_view.vector);
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
    if (gvec == nullptr)
        return;
    gsl_vector_free(gvec);
}

void gsl::vector::set(size_t i, double val) { *gsl_vector_ptr(gvec, i) = val; }

double gsl::vector::get(size_t i) const { return *gsl_vector_ptr(gvec, i); }

double gsl::vector::operator()(size_t i) const { return *gsl_vector_ptr(gvec, i); }

double &gsl::vector::operator()(size_t i) { return *gsl_vector_ptr(gvec, i); }

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

//! \brief Clear the gsl::vector, free underlying memory
void gsl::vector::clear()
{
    if (gvec == nullptr)
        return;
    this->free();
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
void gsl::vector::free()
{
    gsl_vector_free(gvec);
    gvec = nullptr;
}

/*!
 * \brief Private function to (continuously) allocate memory
 * \note This method allocates contiguous, zero-initialized memory.
 *       This is slightly slower than using gsl_vector_alloc, but allows
 *       for intuitive usage of row views.
 */
void gsl::vector::calloc(size_t n) { gvec = gsl_vector_calloc(n); }

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

    vector operator+(const vector &v1, const vector_view &v2)
    {
        vector result(v1);
        result += v2;
        return result;
    }

    vector operator+(const vector_view &v1, const vector &v2)
    {
        vector result(v2);
        result += v1;
        return result;
    }

    vector operator+(vector &&v1, const vector_view &v2)
    {
        v1 += v2;
        return std::move(v1);
    }

    vector operator+(const vector_view &v1, vector &&v2)
    {
        v2 += v1;
        return std::move(v2);
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

    cvector operator+(const vector &v1, const cvector_view &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator+(const cvector_view &v1, const vector &v2)
    {
        cvector result(v2);
        result += v1;
        return result;
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

    vector operator-(const vector &v1, const vector_view &v2)
    {
        vector result(v1);
        result += v2;
        return result;
    }

    vector operator-(const vector_view &v1, const vector &v2)
    {
        vector result(v2);
        result -= v1;
        return result;
    }

    vector operator-(vector &&v1, const vector_view &v2)
    {
        v1 -= v2;
        return std::move(v1);
    }

    vector operator-(const vector_view &v1, vector &&v2)
    {
        v2 -= v1;
        return std::move(v2);
    }

    vector operator-(const vector_view &v1, const vector_view &v2)
    {
        vector result(v1);
        result -= v2;
        return result;
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

    // Compare vectors to vector views
    bool operator==(const vector_view &v1, const vector &v2)
    {
        return gsl_vector_equal(&v1.gvec_view.vector, v2.gvec);
    }

    bool operator==(const vector &v1, const vector_view &v2)
    {
        return gsl_vector_equal(v1.gvec, &v2.gvec_view.vector);
    }

    bool operator!=(const vector_view &v1, const vector &v2)
    {
        // Compare each element, returning true if any are not equal
        for (size_t i = 0; i < v1.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    bool operator!=(const vector &v1, const vector_view &v2)
    {
        // Compare each element, returning true if any are not equal
        for (size_t i = 0; i < v2.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    vector operator*(const vector_view &v, double a)
    {
        vector result(v);
        result *= a;
        return result;
    }

    vector operator*(double a, const vector_view &v)
    {
        vector result(v);
        result *= a;
        return result;
    }

    vector operator/(const vector_view &v, double a)
    {
        vector result(v);
        result /= a;
        return result;
    }

    vector operator/(double a, const vector_view &v)
    {
        vector result(v);
        result /= a;
        return result;
    }

    // Add vector views to vector views
    vector operator+(const vector_view &v1, const vector_view &v2)
    {
        vector result(v1);
        result += v2;
        return result;
    }

    bool operator==(const vector_view &v1, const vector_view &v2)
    {
        return gsl_vector_equal(&v1.gvec_view.vector, &v2.gvec_view.vector);
    }

    bool operator!=(const vector_view &v1, const vector_view &v2)
    {
        // Compare each element, returning true if any are not equal
        for (size_t i = 0; i < v1.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

}

/*! \brief Return a view to a subvector of the vector
 * \param offset Offset of the subvector
 * \param size Size of the subvector
 *
 * \note \return A gsl::vector_view to the subvector
 */
gsl::vector_view gsl::vector::subvector(size_t offset, size_t size)
{
    return gsl::vector_view(gsl_vector_subvector(gvec, offset, size));
}

//! \brief Return a view to the entire vector
gsl::vector_view gsl::vector::view()
{
    return gsl::vector_view(gsl_vector_subvector(gvec, 0, gvec->size));
}

//! \brief Assignment to a vector view from a vector
gsl::vector_view &gsl::vector_view::operator=(const gsl::vector &v)
{
    gsl_vector_memcpy(&gvec_view.vector, v.gvec);
    return *this;
}

//! \brief Assignment to a vector view from a vector view
gsl::vector_view &gsl::vector_view::operator=(vector_view v)
{
    gsl_vector_memcpy(&gvec_view.vector, &v.gvec_view.vector);
    return *this;
}

gsl::vector_view &gsl::vector_view::operator+=(const gsl::vector_view &v)
{
    gsl_vector_add(&gvec_view.vector, &v.gvec_view.vector);
    return *this;
}

gsl::vector_view &gsl::vector_view::operator-=(const gsl::vector_view &v)
{
    gsl_vector_sub(&gvec_view.vector, &v.gvec_view.vector);
    return *this;
}

gsl::vector_view &gsl::vector_view::operator+=(const gsl::vector &v)
{
    gsl_vector_add(&gvec_view.vector, v.gvec);
    return *this;
}

gsl::vector_view &gsl::vector_view::operator-=(const gsl::vector &v)
{
    gsl_vector_sub(&gvec_view.vector, v.gvec);
    return *this;
}

gsl::vector_view &gsl::vector_view::operator*=(double x)
{
    gsl_vector_scale(&gvec_view.vector, x);
    return *this;
}

gsl::vector_view &gsl::vector_view::operator/=(double x)
{
    gsl_vector_scale(&gvec_view.vector, 1.0 / x);
    return *this;
}

double gsl::vector_view::operator()(size_t i) const
{
// This is the only way to get around needing gsl_vector_ptr requiring
//   a non-const gsl_vector pointer. Would rather do anything else.
#if GSL_RANGE_CHECK
    if (GSL_RANGE_COND(i >= gvec_view.vector.size))
    {
        GSL_ERROR_NULL("index out of range", GSL_EINVAL);
    }
#endif
    return *(double *)(gvec_view.vector.data + i * gvec_view.vector.stride);
}

double &gsl::vector_view::operator()(size_t i)
{
    return *gsl_vector_ptr(&gvec_view.vector, i);
}

void gsl::vector_view::set(size_t i, double x)
{
    gsl_vector_set(&gvec_view.vector, i, x);
}

double gsl::vector_view::get(size_t i) const
{
    return gsl_vector_get(&gvec_view.vector, i);
}

//! \brief Pretty-print the viewed vector to file stream
//! \param out File stream to print to
void gsl::vector_view::print(FILE *out) const
{
    fprintf(out, "[");
    for (int i = 0; i < gvec_view.vector.size; ++i)
        fprintf(out, "%s%g", (i == 0) ? "" : ", ", gsl_vector_get(&gvec_view.vector, i));
    fprintf(out, "]\n");
}

//! \brief Return a matrix view out of the elements of the row
gsl::matrix_view gsl::row_view::reshape(size_t n, size_t m)
{
    return gsl::matrix_view(gsl_matrix_view_vector(&gvec_view.vector, n, m));
}