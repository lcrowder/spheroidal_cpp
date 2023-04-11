#include <gsl_wrapper/core.h>
#include <gsl/gsl_vector.h>
#include <fmt/core.h>
#include <complex.h>
#include <gsl/gsl_complex.h>

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
gsl::cvector::cvector(const gsl::vector& gvec_other)
{
    this->calloc(gvec_other.size());
    for (size_t i = 0; i < gvec_other.size(); i++)
        this->set(i, gsl::complex(gvec_other.get(i), 0.0));
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

//! \brief Destructor
gsl::cvector::~cvector()
{
    if (gvec == nullptr)
        return;
    gsl_vector_complex_free(gvec);
}

//! \brief Element getter (the nice C++ versions don't work)
void gsl::cvector::set(size_t i, const gsl::complex& val)
{
    GSL_SET_COMPLEX( gsl_vector_complex_ptr( gvec, i ), val.real(), val.imag() );
}
// const gsl_complex & gsl::cvector::operator()(size_t i) const
// {
//     return *gsl_vector_complex_ptr(gvec, i);
// }

//! \brief Element setter (the nice C++ versions don't work)
gsl::complex gsl::cvector::get(size_t i) const
{
    return gsl::complex( *gsl_vector_complex_ptr( gvec, i ) );
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
void gsl::cvector::print() const
{
    fmt::print("[");
    for (int i = 0; i < gvec->size; ++i)
    {
        auto x = gsl_vector_complex_get(gvec, i);
        fmt::print((i == 0) ? "{:g}{:+g}j" : ", {:g}{:+g}j", GSL_REAL(x), GSL_IMAG(x));
    }
    fmt::print("]\n");
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
 * \note This method allocates contiguous, zero-initialized memory.
 *         This is slightly slower than using gsl_vector_alloc
 */
void gsl::cvector::calloc(size_t n) { gvec = gsl_vector_complex_calloc(n); }