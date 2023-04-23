#include <yawg/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <stdio.h>
#include <utility>

//! \brief Construct empty cvector
gsl::cvector::cvector()
{
    gvec = new gsl_vector_complex;
    gvec->size = 0;
    gvec->data = 0;
}

//! \brief Construct zero cvector of size n
gsl::cvector::cvector(size_t n)
{
    this->galloc(n);
}

//! \brief Construct new gsl::cvector from gsl_vector_complex's data
gsl::cvector::cvector(gsl_vector_complex *gvec_other) : gvec(gvec_other) {}

gsl::cvector::cvector(const gsl::cvector &gvec_other)
{
    this->galloc(gvec_other.size());
    gsl_vector_complex_memcpy(gvec, gvec_other.gvec);
}

gsl::cvector::cvector(gsl::cvector &&gvec_other) : gvec(gvec_other.gvec)
{
    gvec_other.gvec = nullptr;
}

/*! \brief Construct new gsl::cmatrix from gsl::matrix
 *
 *  Copy the values from a real matrix into a complex matrix, setting
 *  the imaginary part to zero.
 */
gsl::cvector::cvector(const gsl::vector &v)
{
    galloc(v.size());
    for (size_t i = 0; i < v.size(); i++)
        set(i, gsl::complex(v(i), 0.0));
}


gsl::cvector &gsl::cvector::operator=(const gsl::cvector &v)
{
    if (this == &v)
        return *this;
    if (size() != v.size())
        this->resize(v.size());
    gsl_vector_complex_memcpy(gvec, v.gvec);
    return *this;
}

gsl::cvector &gsl::cvector::operator=(const gsl::vector &v)
{
    resize(v.size());
    for (size_t i = 0; i < v.size(); i++)
        set(i, gsl::complex(v(i), 0.0));
    return *this;
}

gsl::cvector &gsl::cvector::operator=(gsl::cvector &&v)
{
    if (this == &v)
        return *this;
    this->gfree();
    gvec = v.gvec;
    v.gvec = nullptr;
    return *this;
}

gsl::cvector &gsl::cvector::operator+=(const gsl::cvector &v)
{
    gsl_vector_complex_add(gvec, v.gvec);
    return *this;
}

gsl::cvector &gsl::cvector::operator+=(const gsl::vector &v)
{
    // Add the real element of v to each complex element of this
    for (size_t i = 0; i < gvec->size; i++)
        set(i, get(i) + v(i));
    return *this;
}

gsl::cvector &gsl::cvector::operator-=(const gsl::cvector &gvec_other)
{
    gsl_vector_complex_sub(gvec, gvec_other.gvec);
    return *this;
}

gsl::cvector &gsl::cvector::operator-=(const gsl::vector &v)
{
    // Subtract the real element of v from each complex element of this
    for (size_t i = 0; i < gvec->size; i++)
        set(i, get(i) - v(i));
    return *this;
}

gsl::cvector &gsl::cvector::operator*=(gsl::complex z)
{
    gsl_vector_complex_scale(gvec, z);
    return *this;
}

gsl::cvector &gsl::cvector::operator*=(double x)
{
    // For each element, multiply it by x
    for (size_t i = 0; i < gvec->size; i++)
        gsl_vector_complex_set(gvec, i, gsl_complex_mul_real(gsl_vector_complex_get(gvec, i), x));

    return *this;
}

gsl::cvector &gsl::cvector::operator/=(gsl::complex z)
{
    gsl_vector_complex_scale(gvec, gsl_complex_inverse(z));
    return *this;
}

gsl::cvector &gsl::cvector::operator/=(double x)
{
    // For each element, divide it by x
    for (size_t i = 0; i < gvec->size; i++)
        gsl_vector_complex_set(gvec, i, gsl_complex_div_real(gsl_vector_complex_get(gvec, i), x));

    return *this;
}

gsl::cvector gsl::cvector::operator-() const
{
    gsl::cvector v(*this);
    return -1.0 * v;
}

gsl::cvector::~cvector()
{
    this->gfree();
}

void gsl::cvector::set(size_t i, gsl::complex val)
{
    GSL_SET_COMPLEX(gsl_vector_complex_ptr(gvec, i), val.real(), val.imag());
}

gsl::complex gsl::cvector::get(size_t i) const
{
    return gsl::complex(*gsl_vector_complex_ptr(gvec, i));
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

/*! \brief Return a reference to the element at position (i,j)
 *
 *  This function returns a complex_ref to the element
 *  at position (i,j) in the matrix. Allows setting.
 */
gsl::complex_ref gsl::cvector::operator()(size_t i)
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
        this->gfree();
    }
    this->galloc(n);
}

//! \brief Clear the gsl::vector, free underlying memory
void gsl::cvector::clear()
{
    if (gvec == nullptr)
        return;
    this->gfree();
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
void gsl::cvector::gfree()
{
    if (gvec == nullptr)
        return;
    else if (size() == 0)
        delete gvec;
    else
        gsl_vector_complex_free(gvec);
    gvec = nullptr;
}

/*!
 * \brief Private function to (continuously) allocate memory
 * \note This method allocates contiguous, zero-initialized memory.
 *       This is slightly slower than using gsl_vector_complex_alloc, but allows
 *       for intuitive usage of row views.
 */
void gsl::cvector::galloc(size_t n)
{
    if (n >= 0)
        gvec = gsl_vector_complex_alloc(n);
    else
    {
        gvec = new gsl_vector_complex;
        gvec->size = 0;
        gvec->data = 0;
    }
}
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

    cvector operator*(complex a, const vector &v)
    {
        cvector result(v);
        result *= a;
        return result;
    }

    cvector operator*(const vector &v, complex a)
    {
        cvector result(v);
        result *= a;
        return result;
    }

    cvector operator*(const cvector &v, double x)
    {
        cvector result(v);
        result *= x;
        return result;
    }

    cvector operator*(double x, const cvector &v)
    {
        cvector result(v);
        result *= x;
        return result;
    }

    cvector operator*(cvector &&v, double x)
    {
        v *= x;
        return std::move(v);
    }

    cvector operator*(double x, cvector &&v)
    {
        v *= x;
        return std::move(v);
    }

    cvector operator/(const cvector &v, complex z)
    {
        cvector result(v);
        result /= z;
        return result;
    }

    cvector operator/(cvector &&v, complex z)
    {
        v /= z;
        return std::move(v);
    }

    cvector operator/(complex z, const cvector &v)
    {
        cvector result(v);
        result /= z;
        return result;
    }

    cvector operator/(complex z, cvector &&v)
    {
        v /= z;
        return std::move(v);
    }

    cvector operator/(const vector &v, complex z)
    {
        cvector result(v);
        result /= z;
        return result;
    }

    cvector operator/(complex z, const vector &v)
    {
        cvector result(v);
        result /= z;
        return result;
    }

    cvector operator/(const cvector &v, double x)
    {
        cvector result(v);
        result /= x;
        return result;
    }

    cvector operator/(cvector &&v, double x)
    {
        v /= x;
        return std::move(v);
    }

    cvector operator/(double x, const cvector &v)
    {
        cvector result(v);
        result /= x;
        return result;
    }

    cvector operator/(double x, cvector &&v)
    {
        v /= x;
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

    bool operator==(const cvector &v1, const cvector &v2)
    {
        return gsl_vector_complex_equal(v1.gvec, v2.gvec);
    }

    bool operator!=(const cvector &v1, const cvector &v2)
    {
        // Compare each element of v1 and v2 individually
        for (size_t i = 0; i < v1.gvec->size; ++i)
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

gsl::cvector::operator cvector_view() const
{   
    printf("Creating view from complex vector\n");
    return view();
}

gsl::cvector_view gsl::cvector::view() const
{
    gsl_vector_complex *v = static_cast<gsl_vector_complex *>(malloc(sizeof(gsl_vector_complex)));
    *v = gsl_vector_complex_subvector(get(), 0, size()).vector;
    return gsl::cvector_view(v);
}

gsl::cvector_view gsl::cvector::subvector(size_t offset, size_t n) const
{
    gsl_vector_complex *v = static_cast<gsl_vector_complex *>(malloc(sizeof(gsl_vector_complex)));
    *v = gsl_vector_complex_subvector(get(), offset, n).vector;
    return gsl::cvector_view(v);
}

//! \brief Construct new gsl::cvector from gsl_vector_complex
gsl::cvector_view::cvector_view(gsl_vector_complex *gvec_other) : cvector(gvec_other) {}

// Vector views should never own their memory, so we don't need to free it
// We only need to delete the pointer, since it was created on the heap
gsl::cvector_view::~cvector_view()
{
    if (gvec != nullptr)
        free(gvec);
    gvec = nullptr;
} 

//! Set all values in the view to zero
void gsl::cvector_view::clear()
{
    printf("Warning: Attempting to clear a vector view\n");
    gsl_vector_complex_set_zero(gvec);
}

//! Set all values in the view to zero
void gsl::cvector_view::resize(size_t n)
{
    printf("Warning: Attempting to resize a vector view\n");
    gsl_vector_complex_set_zero(gvec);
}

gsl::cmatrix_view gsl::crow_view::reshape(size_t n, size_t m) const
{
    gsl_matrix_complex *t = static_cast<gsl_matrix_complex *>(malloc(sizeof(gsl_matrix_complex)));
    *t = gsl_matrix_complex_view_vector(get(), n, m).matrix;
    return gsl::cmatrix_view(t);
}