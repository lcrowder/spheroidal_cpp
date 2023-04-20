#include <yawg/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <stdio.h>
#include <utility>

//! \brief Construct empty cvector
gsl::cvector::cvector() : gvec(nullptr) {}

//! \brief Construct zero cvector of size n
gsl::cvector::cvector(size_t n) : gvec(nullptr)
{
    if (n != 0)
        this->calloc(n);
}

//! \brief Construct new gsl::cvector from gsl_vector_complex
gsl::cvector::cvector(const gsl_vector_complex *gvec_other)
{
    if (gvec_other == nullptr)
        return;
    this->calloc(gvec_other->size);
    gsl_vector_complex_memcpy(gvec, gvec_other);
}

gsl::cvector::cvector(const gsl::cvector &gvec_other)
{
    this->calloc(gvec_other.size());
    gsl_vector_complex_memcpy(gvec, gvec_other.gvec);
}

/*! \brief Construct new gsl::cmatrix from gsl::matrix
 *
 *  Copy the values from a real matrix into a complex matrix, setting
 *  the imaginary part to zero.
 */
gsl::cvector::cvector(const gsl::vector &gvec_other)
{
    this->calloc(gvec_other.size());
    for (size_t i = 0; i < gvec_other.size(); i++)
        GSL_SET_COMPLEX(gsl_vector_complex_ptr(gvec, i), gsl_vector_get(gvec_other.gvec, i), 0.0);
}

gsl::cvector::cvector(gsl::cvector &&gvec_other) : gvec(gvec_other.gvec)
{
    gvec_other.gvec = nullptr;
}

gsl::cvector &gsl::cvector::operator=(const gsl::cvector &gvec_other)
{
    if (this == &gvec_other)
        return *this;
    this->resize(gvec_other.size());
    gsl_vector_complex_memcpy(gvec, gvec_other.gvec);
    return *this;
}

gsl::cvector &gsl::cvector::operator=(gsl::cvector &&v)
{
    if (this == &v)
        return *this;
    this->free();
    gvec = v.gvec;
    v.gvec = nullptr;
    return *this;
}

gsl::cvector &gsl::cvector::operator+=(const gsl::cvector &v)
{
    gsl_vector_complex_add(gvec, v.gvec);
    return *this;
}

gsl::cvector &gsl::cvector::operator+=(const gsl::cvector_view &vv)
{
    gsl_vector_complex_add(gvec, &vv.gvec_view.vector);
    return *this;
}

gsl::cvector &gsl::cvector::operator+=(const gsl::vector &v)
{
    // Add the real element of v to each complex element of this
    for (size_t i = 0; i < gvec->size; i++)
        gsl_vector_complex_set(gvec, i, gsl_complex_add_real(gsl_vector_complex_get(gvec, i), gsl_vector_get(v.gvec, i)));
    return *this;
}

gsl::cvector &gsl::cvector::operator+=(const gsl::vector_view &vv)
{
    // Add the real element of v to each complex element of this
    for (size_t i = 0; i < gvec->size; i++)
        gsl_vector_complex_set(gvec, i, gsl_complex_add_real(gsl_vector_complex_get(gvec, i), gsl_vector_get(&vv.gvec_view.vector, i)));
    return *this;
}

gsl::cvector &gsl::cvector::operator-=(const gsl::cvector &gvec_other)
{
    gsl_vector_complex_sub(gvec, gvec_other.gvec);
    return *this;
}

gsl::cvector &gsl::cvector::operator-=(const gsl::cvector_view &vv)
{
    gsl_vector_complex_sub(gvec, &vv.gvec_view.vector);
    return *this;
}

gsl::cvector &gsl::cvector::operator-=(const gsl::vector &v)
{
    // Subtract the real element of v from each complex element of this
    for (size_t i = 0; i < gvec->size; i++)
        gsl_vector_complex_set(gvec, i, gsl_complex_sub_real(gsl_vector_complex_get(gvec, i), gsl_vector_get(v.gvec, i)));
    return *this;
}

gsl::cvector &gsl::cvector::operator-=(const gsl::vector_view &vv)
{
    // Subtract the real element of v from each complex element of this
    for (size_t i = 0; i < gvec->size; i++)
        gsl_vector_complex_set(gvec, i, gsl_complex_sub_real(gsl_vector_complex_get(gvec, i), gsl_vector_get(&vv.gvec_view.vector, i)));
    return *this;
}

gsl::cvector::~cvector()
{
    if (gvec == nullptr)
        return;
    gsl_vector_complex_free(gvec);
}

//! \brief Element getter (the nice C++ versions don't work)
void gsl::cvector::set(size_t i, gsl::complex val)
{
    GSL_SET_COMPLEX(gsl_vector_complex_ptr(gvec, i), val.real(), val.imag());
}

//! \brief Element setter (the nice C++ versions don't work)
gsl::complex gsl::cvector::get(size_t i) const
{
    return gsl::complex(*gsl_vector_complex_ptr(gvec, i));
}

gsl::cvector &gsl::cvector::operator*=(gsl::complex z)
{
    gsl_vector_complex_scale(gvec, z);
    return *this;
}

gsl::cvector &gsl::cvector::operator/=(gsl::complex z)
{
    gsl_vector_complex_scale(gvec, gsl_complex_inverse(z));
    return *this;
}

gsl::cvector gsl::cvector::operator-() const
{
    gsl::cvector v(*this);
    return -1.0 * v;
}

/*! \brief Return a reference to the element at position (i,j)
 *
 *  This function returns a complex_ref to the element
 *  at position (i,j) in the matrix. Allows setting.
 */
gsl::complex_ref gsl::cvector::operator()(size_t i)
{
    return gsl::complex_ref(gsl_vector_complex_ptr(gvec, i));
}

/*! \brief Return a const reference to the element at position (i,j)
 *
 *  This function returns a constant complex_ref to the element
 *  at position (i,j) in the matrix. Allows getting.
 */
const gsl::complex_ref gsl::cvector::operator()(size_t i) const
{
    return gsl::complex_ref(gsl_vector_complex_ptr(gvec, i));
}

size_t gsl::cvector::size() const
{
    if (gvec == nullptr)
        return 0;
    return gvec->size;
};

/*! \brief Resize the gsl::cvector, setting elements to zero
 * \param n Number of elements
 *
 * \note This function will always free and reallocate memory,
 *  setting the elements to zero.
 */
void gsl::cvector::resize(size_t n)
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
void gsl::cvector::clear()
{
    if (gvec == nullptr)
        return;
    this->free();
}

//! \brief Pretty-print the vector to file stream
//! \param out File stream to print to
void gsl::cvector::print(FILE *out) const
{
    fprintf(out, "[");
    for (int i = 0; i < gvec->size; ++i)
    {
        auto x = gsl_vector_complex_get(gvec, i);
        fprintf(out, "%s% 9g%+9gj", ((i == 0) ? "" : ", "), GSL_REAL(x), GSL_IMAG(x));
    }
    fprintf(out, "]\n");
}

//! \brief Private function to free allocated memory
void gsl::cvector::free()
{
    gsl_vector_complex_free(gvec);
    gvec = nullptr;
}

/*!
 * \brief Private function to (continuously) allocate memory
 * \note This method allocates contiguous, zero-initialized memory.
 *       This is slightly slower than using gsl_vector_complex_alloc, but allows
 *       for intuitive usage of row views.
 */
void gsl::cvector::calloc(size_t n) { gvec = gsl_vector_complex_calloc(n); }

namespace gsl
{
    cvector operator*(complex z, const cvector &v)
    {
        cvector result(v);
        result *= z;
        return result;
    }

    cvector operator*(complex z, cvector &&v)
    {
        v *= z;
        return std::move(v);
    }

    cvector operator*(const cvector &v, complex z)
    {
        cvector result(v);
        result *= z;
        return result;
    }

    cvector operator*(cvector &&v, complex z)
    {
        v *= z;
        return std::move(v);
    }

    cvector operator+(const cvector &v1, const cvector &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator+(cvector &&v1, const cvector &v2)
    {
        v1 += v2;
        return v1;
    }

    cvector operator+(const cvector &v1, cvector &&v2)
    {
        v2 += v1;
        return v2;
    }

    cvector operator+(cvector &&v1, cvector &&v2)
    {
        v1 += v2;
        return v1;
    }

    cvector operator+(const vector_view &v1, const cvector &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator+(const vector_view &v1, cvector &&v2)
    {
        v2 += v1;
        return v2;
    }

    cvector operator+(const cvector &v1, const vector_view &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator+(cvector &&v1, const vector_view &v2)
    {
        v1 += v2;
        return v1;
    }

    // Add vector views to complepx vector views
    cvector operator+(const vector_view &v1, const cvector_view &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator+(const cvector_view &v1, const vector_view &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator+(const cvector &v1, const cvector_view &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator+(cvector &&v1, const cvector_view &v2)
    {
        v1 += v2;
        return std::move(v1);
    }

    cvector operator+(const cvector_view &v1, const cvector &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator+(const cvector_view &v1, cvector &&v2)
    {
        v2 += v1;
        return std::move(v2);
    }

    cvector operator+(const cvector_view &v1, const cvector_view &v2)
    {
        cvector result(v1);
        result += v2;
        return result;
    }

    cvector operator-(const cvector &v1, const cvector &v2)
    {
        cvector result(v1);
        result -= v2;
        return result;
    }

    cvector operator-(cvector &&v1, const cvector &v2)
    {
        v1 -= v2;
        return v1;
    }

    cvector operator-(const cvector &v1, cvector &&v2)
    {
        v2 -= v1;
        return v2;
    }

    cvector operator-(cvector &&v1, cvector &&v2)
    {
        v1 -= v2;
        return v1;
    }

    cvector operator-(const vector &v1, const cvector &v2)
    {
        cvector result(v1);
        result -= v2;
        return result;
    }

    cvector operator-(const vector &v1, cvector &&v2)
    {
        v2 -= v1;
        return std::move(v2);
    }

    cvector operator-(const cvector &v1, const vector &v2)
    {
        cvector result(v2);
        result -= v1;
        return result;
    }

    cvector operator-(cvector &&v1, const vector &v2)
    {
        v1 -= v2;
        return std::move(v1);
    }

    cvector operator-(const vector &v1, const cvector_view &v2)
    {
        cvector result(v1);
        result -= v2;
        return result;
    }

    cvector operator-(const cvector_view &v1, const vector &v2)
    {
        cvector result(v2);
        result -= v1;
        return result;
    }

    cvector operator-(const vector_view &v1, const cvector_view &v2)
    {
        cvector result(v1);
        result -= v2;
        return result;
    }

    cvector operator-(const cvector_view &v1, const vector_view &v2)
    {
        cvector result(v2);
        result -= v1;
        return result;
    }

    cvector operator-(const vector_view &v1, const cvector &v2)
    {
        cvector result(v1);
        result -= v2;
        return result;
    }

    cvector operator-(const vector_view &v1, cvector &&v2)
    {
        v2 -= v1;
        return std::move(v2);
    }

    cvector operator-(const cvector &v1, const vector_view &v2)
    {
        cvector result(v1);
        result -= v2;
        return result;
    }

    cvector operator-(cvector &&v1, const vector_view &v2)
    {
        v1 -= v2;
        return std::move(v1);
    }

    cvector operator-(const cvector &v1, const cvector_view &v2)
    {
        cvector result(v1);
        result -= v2;
        return result;
    }

    cvector operator-(cvector &&v1, const cvector_view &v2)
    {
        v1 -= v2;
        return std::move(v1);
    }

    cvector operator-(const cvector_view &v1, const cvector &v2)
    {
        cvector result(v1);
        result -= v2;
        return result;
    }

    cvector operator-(const cvector_view &v1, cvector &&v2)
    {
        v2 -= v1;
        return std::move(v2);
    }

    bool operator==(const cvector &v1, const cvector &v2)
    {
        return gsl_vector_complex_equal(v1.gvec, v2.gvec);
    }

    bool operator==(const vector_view &v1, const cvector &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v1.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const vector_view &v1, const cvector &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v1.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    bool operator==(const cvector &v1, const vector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v2.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const cvector &v1, const vector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v2.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    bool operator==(const vector_view &v1, const cvector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v1.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const vector_view &v1, const cvector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v1.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    bool operator==(const cvector_view &v1, const vector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v2.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const cvector_view &v1, const vector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v2.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    bool operator==(const cvector_view &v1, const cvector &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v1.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const cvector_view &v1, const cvector &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v1.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    bool operator==(const cvector &v1, const cvector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v2.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const cvector &v1, const cvector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v2.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    bool operator==(const cvector_view &v1, const cvector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v2.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const cvector_view &v1, const cvector_view &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v2.gvec_view.vector.size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    // Compare vectors to complex vector views
    bool operator==(const vector &v1, const cvector_view &v2)
    {
        // Compare each element of v1 to v2, returning true if any are unequal
        for (size_t i = 0; i < v1.gvec->size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const vector &v1, const cvector_view &v2)
    {
        // Compare each element of v1 to v2, returning true if any are unequal
        for (size_t i = 0; i < v1.gvec->size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    bool operator==(const cvector_view &v1, const vector &v2)
    {
        // Compare each element of v1 to v2, returning true if any are unequal
        for (size_t i = 0; i < v2.gvec->size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const cvector_view &v1, const vector &v2)
    {
        // Compare each element of v1 to v2, returning true if any are unequal
        for (size_t i = 0; i < v2.gvec->size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    // Compare vectors to complex vectors
    bool operator==(const vector &v1, const cvector &v2)
    {
        // Compare each element of v1 to v2, returning false if any are unequal
        for (size_t i = 0; i < v1.gvec->size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const vector &v1, const cvector &v2)
    {
        // Compare each element of v1 to v2, returning true if any are unequal
        for (size_t i = 0; i < v1.gvec->size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

    bool operator==(const cvector &v1, const vector &v2)
    {
        // Compare each element of v1 to v2, returning false if any are unequal
        for (size_t i = 0; i < v1.gvec->size; ++i)
            if (v1(i) != v2(i))
                return false;
        return true;
    }

    bool operator!=(const cvector &v1, const vector &v2)
    {
        // Compare each element of v1 to v2, returning true if any are unequal
        for (size_t i = 0; i < v1.gvec->size; ++i)
            if (v1(i) != v2(i))
                return true;
        return false;
    }

}

/*! \brief Return a view to a subvector of the cvector
 * \param offset Offset of the subvector
 * \param size Size of the subvector
 *
 * \note \return A gsl::cvector_view to the subvector
 */
gsl::cvector_view gsl::cvector::subvector(size_t offset, size_t size)
{
    return gsl::cvector_view(gsl_vector_complex_subvector(gvec, offset, size));
}

//! \brief Return a view to the entire cvector
gsl::cvector_view gsl::cvector::view()
{
    return gsl::cvector_view(gsl_vector_complex_subvector(gvec, 0, gvec->size));
}

//! \brief Assignment to a cvector view from a cvector
gsl::cvector_view &gsl::cvector_view::operator=(const gsl::cvector &v)
{
    gsl_vector_complex_memcpy(&gvec_view.vector, v.gvec);
    return *this;
}

//! \brief Assignment to a cvector view from a cvector view
gsl::cvector_view &gsl::cvector_view::operator=(cvector_view v)
{
    gsl_vector_complex_memcpy(&gvec_view.vector, &v.gvec_view.vector);
    return *this;
}

const gsl::complex_ref gsl::cvector_view::operator()(size_t i) const
{
// This is the only way to get around needing gsl_vector_complex_ptr requiring
//   a non-const gsl_vector_complex pointer. Would rather do anything else.
#if GSL_RANGE_CHECK
    if (GSL_RANGE_COND(i >= gvec_view.vector.size))
    {
        GSL_ERROR_NULL("index out of range", GSL_EINVAL);
    }
#endif
    return gsl::complex_ref(GSL_COMPLEX_AT(&gvec_view.vector, i));
}

gsl::complex_ref gsl::cvector_view::operator()(size_t i)
{
    return gsl::complex_ref(gsl_vector_complex_ptr(&gvec_view.vector, i));
}

void gsl::cvector_view::set(size_t i, complex z)
{
    gsl_vector_complex_set(&gvec_view.vector, i, z);
}

gsl::complex gsl::cvector_view::get(size_t i) const
{
    return gsl_vector_complex_get(&gvec_view.vector, i);
}

//! \brief Pretty-print the viewed vvector to file stream
//! \param out File stream to print to
void gsl::cvector_view::print(FILE *out) const
{
    fprintf(out, "[");
    for (int i = 0; i < gvec_view.vector.size; ++i)
    {
        auto x = gsl_vector_complex_get(&gvec_view.vector, i);
        fprintf(out, "%s% 9g%+9gj", ((i == 0) ? "" : ", "), GSL_REAL(x), GSL_IMAG(x));
    }
    fprintf(out, "]\n");
}

//! \brief Return a cmatrix view out of the elements of the row
gsl::cmatrix_view gsl::crow_view::reshape(size_t n, size_t m)
{
    return gsl::cmatrix_view(gsl_matrix_complex_view_vector(&gvec_view.vector, n, m));
}