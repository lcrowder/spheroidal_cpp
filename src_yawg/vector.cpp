#include <yawg/core.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <utility>

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
gsl::vector::vector(const gsl_vector *gvec_other)
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

//! \brief Addition assignment operator
gsl::vector &gsl::vector::operator+=(const gsl::vector &gvec_other)
{
    gsl_vector_add(gvec, gvec_other.gvec);
    return *this;
}

//! \brief Subtraction assignment operator
gsl::vector &gsl::vector::operator-=(const gsl::vector &gvec_other)
{
    gsl_vector_sub(gvec, gvec_other.gvec);
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

gsl::vector_view gsl::vector::subvector(size_t offset, size_t size)
{
    return gsl::vector_view(gsl_vector_subvector(gvec, offset, size));
}

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
void gsl::vector::print(FILE *out) const
{
    fprintf(out, "[");
    for (int i = 0; i < gvec->size; ++i)
        fprintf(out, "%s%g", (i == 0) ? "" : ", ", gsl_vector_get(gvec, i));
    fprintf(out, "]\n");
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

/*------ friend operators ------*/
namespace gsl
{
    vector operator*(double a, const vector &v)
    {
        vector result(v);
        gsl_vector_scale(result.gvec, a);
        return result;
    }

    vector operator*(double a, vector &&v)
    {
        gsl_vector_scale(v.gvec, a);
        return std::move(v);
    }

    vector operator*(const vector &v, double a)
    {
        return a * v;
    }

    vector operator*(vector &&v, double a)
    {
        return a * std::move(v);
    }

    cvector operator*(complex a, const vector &v)
    {
        return cvector(v) * a;
    }

    cvector operator*(const vector &v, complex a)
    {
        return cvector(v) * a;
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

    bool operator==(const vector &v1, const vector &v2)
    {
        return gsl_vector_equal(v1.gvec, v2.gvec);
    }
}

// gsl::vector_view::vector_view( const gsl::vector& v)
// {
//     gsl::vector_view(gsl_vector_subvector(v.gvec, 0, v.gvec->size));
// }

gsl::vector_view &gsl::vector_view::operator=(const gsl::vector &v)
{
    gsl_vector_memcpy(&gvec_view.vector, v.gvec);
    return *this;
}

gsl::vector_view &gsl::vector_view::operator=(vector_view v)
{
    gsl_vector_memcpy(&gvec_view.vector, &v.gvec_view.vector);
    return *this;
}

void gsl::vector_view::print(FILE *out) const
{
    fprintf(out, "[");
    for (int i = 0; i < gvec_view.vector.size; ++i)
        fprintf(out, "%s%g", (i == 0) ? "" : ", ", gsl_vector_get(&gvec_view.vector, i));
    fprintf(out, "]\n");
}

gsl::matrix_view gsl::row_view::reshape( size_t n, size_t m )
{
    return gsl::matrix_view(gsl_matrix_view_vector(&gvec_view.vector, n, m));
}