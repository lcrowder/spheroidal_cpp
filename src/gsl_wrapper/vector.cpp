#include <gsl_wrapper/core.h>
#include <gsl/gsl_vector.h>
#include <fmt/core.h>

/*------ Public Methods for gsl::vector ------*/
//! \brief Default constructor
gsl::vector::vector() : gvec(nullptr) {}

//! \brief Constructor with size
gsl::vector::vector(size_t n) : gvec(nullptr)
{
    if (n != 0)
        this->calloc(n);
}

//! \brief Build a gsl::vector from a gsl_vector
gsl::vector::vector(gsl_vector *gvec_other)
{
    if (gvec_other == nullptr)
        return;
    this->calloc(gvec_other->size);
    gsl_vector_memcpy(gvec, gvec_other);
}

//! \brief Copy constructor
gsl::vector::vector(const gsl::vector &gvec_other)
{
    this->calloc(gvec_other.size());
    gsl_vector_memcpy(gvec, gvec_other.gvec);
}

//! \brief Conversion constructor from gsl::cvector
//! \note Move constructor isn't useful,
//   since contiguous memory can't be reclaimed
gsl::vector::vector(const gsl::cvector &gvec_other)
{
    fmt::print( stderr, "Warning: Conversion from gsl::cvector to gsl::vector discards imaginary part\n" );
    this->calloc(gvec_other.size());
    for (size_t i = 0; i < gvec_other.size(); i++)
        gsl_vector_set( gvec, i, gsl_vector_complex_get( gvec_other.gvec, i ).dat[0] );
}

//! \brief Move constructor
gsl::vector::vector(gsl::vector &&gvec_other) : gvec(gvec_other.gvec)
{
    gvec_other.gvec = nullptr;
}

//! \brief Assignment operator
gsl::vector &gsl::vector::operator=(const gsl::vector &gvec_other)
{
    if (this == &gvec_other)
        return *this;
    this->resize(gvec_other.size());
    gsl_vector_memcpy(gvec, gvec_other.gvec);
    return *this;
}

//! \brief Move assignment operator
gsl::vector &gsl::vector::operator=(gsl::vector &&gvec_other)
{
    if (this == &gvec_other)
        return *this;
    this->free();
    gvec = gvec_other.gvec;
    gvec_other.gvec = nullptr;
    return *this;
}

//! \brief Destructor
gsl::vector::~vector()
{
    if (gvec == nullptr)
        return;
    gsl_vector_free(gvec);
}

//! \brief Element getter
double gsl::vector::operator()(size_t i) const { return *gsl_vector_ptr(gvec, i); }
double gsl::vector::get(size_t i) const { return *gsl_vector_ptr(gvec, i); }

//! \brief Element setter
double &gsl::vector::operator()(size_t i) { return *gsl_vector_ptr(gvec, i); } 
void gsl::vector::set(size_t i, double val) { *gsl_vector_ptr(gvec, i) = val; }

//! \brief Size accessor
size_t gsl::vector::size() const
{
    if (gvec == nullptr)
        return 0;
    return gvec->size;
};

//! \brief Resize the gsl::vector
//! \note If n == 0, the vector is cleared
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

//! \brief Clear the vector and free the underlying gsl_vector
void gsl::vector::clear()
{
    if (gvec == nullptr)
        return;
    this->free();
}

//! \brief Print the vector to stdout
void gsl::vector::print() const
{
    fmt::print("[");
    for (int i = 0; i < gvec->size; ++i)
        fmt::print((i == 0) ? "{:g}" : ", {:g}", gsl_vector_get(gvec, i));
    fmt::print("]\n");
}

/*------ Protected Methods for gsl::vector ------*/

//! \brief Free memory for underlying gsl_vector
void gsl::vector::free()
{
    gsl_vector_free(gvec);
    gvec = nullptr;
}

/*!
 * \brief Allocate memory for underlying gsl_vector
 * \note This method allocates contiguous, zero-initialized memory.
 *         This is slightly slower than using gsl_vector_alloc
 */
void gsl::vector::calloc(size_t n) { gvec = gsl_vector_calloc(n); }