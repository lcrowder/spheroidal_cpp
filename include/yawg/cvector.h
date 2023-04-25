#ifndef YAWG_CVECTOR_H_
#define YAWG_CVECTOR_H_

#include <yawg/complex.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

namespace gsl
{
    class vector;
    class cvector_view;
    class cmatrix_view;

    /*! \class cvector
     *  \brief A wrapper class for gsl_vector_complex
     *
     * Stores and operates on a pointer to a gsl_vector_complex.
     */
    class cvector
    {
        // Scalar multiplication
        friend cvector operator*(complex a, const cvector &v);
        friend cvector operator*(complex a, cvector &&v);
        friend cvector operator*(const cvector &v, complex a);
        friend cvector operator*(cvector &&v, complex a);

        friend cvector operator*(const cvector &v, double x);
        friend cvector operator*(double x, const cvector &v);
        friend cvector operator*(cvector &&v, double x);
        friend cvector operator*(double x, cvector &&v);

        // Scalar division
        friend cvector operator/(complex a, const cvector &v);
        friend cvector operator/(complex a, cvector &&v);
        friend cvector operator/(const cvector &v, complex a);
        friend cvector operator/(cvector &&v, complex a);

        friend cvector operator/(const cvector &v, double x);
        friend cvector operator/(double x, const cvector &v);
        friend cvector operator/(cvector &&v, double x);
        friend cvector operator/(double x, cvector &&v);

        // Add complex vectors to complex vectors
        friend cvector operator+(const cvector &v1, const cvector &v2);
        friend cvector operator+(cvector &&v1, const cvector &v2);
        friend cvector operator+(const cvector &v1, cvector &&v2);
        friend cvector operator+(cvector &&v1, cvector &&v2);

        // Add complex vectors to vectors
        friend cvector operator+(const vector &v1, const cvector &v2);
        friend cvector operator+(const vector &v1, cvector &&v2);
        friend cvector operator+(const cvector &v1, const vector &v2);
        friend cvector operator+(cvector &&v1, const vector &v2);

        // Subtract complex vectors from complex vectors
        friend cvector operator-(const cvector &v1, const cvector &v2);
        friend cvector operator-(cvector &&v1, const cvector &v2);
        friend cvector operator-(const cvector &v1, cvector &&v2);
        friend cvector operator-(cvector &&v1, cvector &&v2);

        // Subtract complex vectors from vectors
        friend cvector operator-(const vector &v1, const cvector &v2);
        friend cvector operator-(const vector &v1, cvector &&v2);
        friend cvector operator-(const cvector &v1, const vector &v2);
        friend cvector operator-(cvector &&v1, const vector &v2);

        // Compare complex vectors to complex vectors
        friend bool operator==(const cvector &v1, const cvector &v2);
        friend bool operator!=(const cvector &v1, const cvector &v2);

        // Compare complex vectors to vectors
        friend bool operator==(const vector &v1, const cvector &v2);
        friend bool operator==(const cvector &v1, const vector &v2);

        friend bool operator!=(const vector &v1, const cvector &v2);
        friend bool operator!=(const cvector &v1, const vector &v2);

        friend cvector operator*(const cmatrix &M, const cvector &v);
    public:
        //! \brief Construct empty vector
        cvector();

        //! \brief Construct zero vector of size n
        explicit cvector(size_t n);

        cvector(gsl_vector_complex *gvec_other);

        cvector(const cvector &gvec_other);
        cvector(cvector &&gvec_other);

        //! \brief Construct new gsl::cvector from gsl::vector
        cvector(const vector &vec);

        cvector &operator=(const cvector &v);
        cvector &operator=(const vector &v);

        cvector &operator=(cvector &&gvec_other);

        cvector &operator+=(const cvector &v);
        cvector &operator+=(const vector &v);

        cvector &operator-=(const cvector &v);
        cvector &operator-=(const vector &v);

        cvector &operator*=(complex a);
        cvector &operator/=(complex a);
        cvector &operator*=(double a);
        cvector &operator/=(double a);
        cvector operator-() const;

        ~cvector();

        //! \brief Return a reference to the element at position (i,j)
        complex_ref operator()(size_t i);
        void set(size_t i, complex z);

        //! \brief Return a const reference to the element at position (i,j)
        const complex_ref operator()(size_t i) const;
        complex get(size_t i) const;

        size_t size() const;

        //! \brief Access the pointer to the underlying gsl_vector_complex
        gsl_vector_complex *get() const { return gvec; }

        //! \brief Resize the gsl::cvector, setting elements to zero
        void resize(size_t n);

        //! \brief Clear the gsl::cvector, free underlying memory
        void clear();

        //! \brief Pretty-print the complex vector to file stream
        void print(FILE *out = stdout) const;

        //! \brief Return the 2-norm of the vector
        double norm() const { return gsl_blas_dznrm2(gvec); }

        cvector_view view() const;
        cvector_view subvector(size_t offset, size_t n) const;

    protected:
        gsl_vector_complex *gvec;

        //! \brief Construct new gsl::cvector from gsl_vector_complex
        cvector(const gsl_vector_complex *gvec_other);

        //! \brief Private function to free allocated memory
        void gfree();

        //! \brief Private function to (continuously) allocate memory
        void galloc(size_t n);
    };

    class cvector_view : public cvector
    {
    public:
        //! \brief Constructor for vector_view pointing to data at gvec_other
        cvector_view(gsl_vector_complex *gvec_other);
        ~cvector_view();

        //! \brief Assign data from cvector to a view
        cvector_view &operator=(const cvector &v);

        // Override some nonconst member functions to be unusable
        void clear();
        void resize(size_t n);
    };

    /*! \class crow_view
     *  \brief A subclass of cvector_view for stride-1 cvectors
     */
    class crow_view : public cvector_view
    {
    public:
        //! \brief Construct row_view from existing vector view, checking that stride is 1
        crow_view(gsl_vector_complex *gvec_other) : cvector_view(gvec_other)
        {
            if (gvec_other->stride != 1)
                printf("Row view must have stride 1");
        };

        //! \brief Assign data from complex vector to view
        crow_view &operator=(const cvector &v);

        //! \brief Return a matrix view out of the elements of the row
        cmatrix_view reshape(size_t n, size_t m) const;
    };

    /*! \class column_view
     *  \brief A subclass of vector_view for non-stride-1 vectors
     */
    class ccolumn_view : public cvector_view
    {
    public:
        //! \brief Construct column_view from existing vector view
        ccolumn_view(gsl_vector_complex *gvec_other) : cvector_view(gvec_other) {}

        //! \brief Assign data from complex vector to view
        ccolumn_view &operator=(const cvector &v);
    };

} // namespace gsl

#endif // YAWG_CVECTOR_H_