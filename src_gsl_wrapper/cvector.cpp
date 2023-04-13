#include <gsl_wrapper/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <stdio.h>
#include <utility>

/*------ Public Methods for gsl::cvector ------*/
//! \brief Default constructor
gsl::cvector::cvector() : gvec(nullptr) {}

//! \brief Constructor with size
gsl::cvector::cvector(size_t n) : gvec(nullptr)
{
    if (n != 0)
        this->calloc(n);
}

//! \brief Build a gsl::cvector from a gsl_vector_complex
gsl::cvector::cvector(gsl_vector_complex *gvec_other)
{
    if (gvec_other == nullptr)
        return;
    this->calloc(gvec_other->size);
    gsl_vector_complex_memcpy(gvec, gvec_other);
}

//! \brief Copy constructor
gsl::cvector::cvector(const gsl::cvector &gvec_other)
{
    this->calloc(gvec_other.size());
    gsl_vector_complex_memcpy(gvec, gvec_other.gvec);
}

//! \brief Conversion constructor from gsl::vector
//! \note Move constructor isn't useful,
//   since contiguous memory can't be reclaimed
gsl::cvector::cvector(const gsl::vector &gvec_other)
{
    this->calloc(gvec_other.size());
    for (size_t i = 0; i < gvec_other.size(); i++)
        GSL_SET_COMPLEX(gsl_vector_complex_ptr(gvec, i), gsl_vector_get(gvec_other.gvec, i), 0.0);
}

//! \brief Move constructor
gsl::cvector::cvector(gsl::cvector &&gvec_other) : gvec(gvec_other.gvec)
{
    gvec_other.gvec = nullptr;
}

//! \brief Assignment operator
gsl::cvector &gsl::cvector::operator=(const gsl::cvector &gvec_other)
{
    if (this == &gvec_other)
        return *this;
    this->resize(gvec_other.size());
    gsl_vector_complex_memcpy(gvec, gvec_other.gvec);
    return *this;
}

//! \brief Move assignment operator
gsl::cvector &gsl::cvector::operator=(gsl::cvector &&gvec_other)
{
    if (this == &gvec_other)
        return *this;
    this->free();
    gvec = gvec_other.gvec;
    gvec_other.gvec = nullptr;
    return *this;
}

//! \brief Addition assignment operator
gsl::cvector &gsl::cvector::operator+=(const gsl::cvector &gvec_other)
{
    gsl_vector_complex_add(gvec, gvec_other.gvec);
    return *this;
}

//! \brief Subtraction assignment operator
gsl::cvector &gsl::cvector::operator-=(const gsl::cvector &gvec_other)
{
    gsl_vector_complex_sub(gvec, gvec_other.gvec);
    return *this;
}

//! \brief Destructor
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
// const gsl_complex & gsl::cvector::operator()(size_t i) const
// {
//     return *gsl_vector_complex_ptr(gvec, i);
// }

//! \brief Element setter (the nice C++ versions don't work)
gsl::complex gsl::cvector::get(size_t i) const
{
    return gsl::complex(*gsl_vector_complex_ptr(gvec, i));
}
// gsl_complex &gsl::cvector::operator()(size_t i)
// {
//     return *gsl_vector_complex_ptr(gvec, i);
// }

//! \brief Size accessor
size_t gsl::cvector::size() const
{
    if (gvec == nullptr)
        return 0;
    return gvec->size;
};

//! \brief Resize the gsl::cvector
//! \note If n == 0, the cvector is cleared
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

//! \brief Clear the cvector and free the underlying gsl_vector_complex
void gsl::cvector::clear()
{
    if (gvec == nullptr)
        return;
    this->free();
}

//! \brief Print the cvector to stdout
void gsl::cvector::print(FILE* out) const
{
    fprintf(out, "[");
    for (int i = 0; i < gvec->size; ++i)
    {
        auto x = gsl_vector_complex_get(gvec, i);
        fprintf(out, "%s% 9g%+9gj", ((i == 0) ? "" : ", "), GSL_REAL(x), GSL_IMAG(x));
    }
    fprintf(out, "]\n");
}

/*------ Protected Methods for gsl::cvector ------*/

//! \brief Free memory for underlying gsl_vector
void gsl::cvector::free()
{
    gsl_vector_complex_free(gvec);
    gvec = nullptr;
}

/*!
 * \brief Allocate memory for underlying gsl_vector
 * \nocomplexte This method allocates contiguous, zero-initialized memory.
 *         This is slightly slower than using gsl_vector_alloc
 */
void gsl::cvector::calloc(size_t n) { gvec = gsl_vector_complex_calloc(n); }

/*------ friend operators ------*/
namespace gsl
{
    cvector operator*(complex a, const cvector &v)
    {
        cvector result(v);
        gsl_vector_complex_scale(v.gvec, a.get_gsl_data());
        return result;
    }

    cvector operator*(complex a, cvector &&v)
    {
        gsl_vector_complex_scale(v.gvec, a.get_gsl_data());
        return v;
    }

    cvector operator*(const cvector &v, complex a)
    {
        return a * v;
    }

    cvector operator*(cvector &&v, complex a)
    {
        return a * std::move( v );
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
    
    bool operator==(const cvector &v1, const cvector &v2)
    {
        return gsl_vector_complex_equal( v1.gvec, v2.gvec );
    }
}