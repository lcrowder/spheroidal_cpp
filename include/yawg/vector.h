#ifndef YAWG_VECTOR_H_
#define YAWG_VECTOR_H_

#include <yawg/complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

namespace gsl
{
    class vector_view;
    class matrix_view;
    
    class cvector;

    /*! \class vector
     *  \brief A wrapper class for gsl_vector
     *
     * Stores and operates on a pointer to a gsl_vector.
     */
    class vector
    {
        // Scalar multiplication
        friend vector operator*(double a, const vector &v);
        friend vector operator*(double a, vector &&v);
        friend vector operator*(const vector &v, double a);
        friend vector operator*(vector &&v, double a);

        friend cvector operator*(complex a, const vector &v);
        friend cvector operator*(const vector &v, complex a);

        // Scalar division
        friend vector operator/(double a, const vector &v);
        friend vector operator/(double a, vector &&v);
        friend vector operator/(const vector &v, double a);
        friend vector operator/(vector &&v, double a);

        friend cvector operator/(complex z, const vector &v);
        friend cvector operator/(const vector &v, complex z);

        // Add vectors to vectors
        friend vector operator+(const vector &v1, const vector &v2);
        friend vector operator+(vector &&v1, const vector &v2);
        friend vector operator+(const vector &v1, vector &&v2);
        friend vector operator+(vector &&v1, vector &&v2);

        // Add vectors to complex vectors
        friend cvector operator+(const vector &v1, const cvector &v2);
        friend cvector operator+(const vector &v1, cvector &&v2);
        friend cvector operator+(const cvector &v1, const vector &v2);
        friend cvector operator+(cvector &&v1, const vector &v2);

        // Subtract vectors from vectors
        friend vector operator-(const vector &v1, const vector &v2);
        friend vector operator-(vector &&v1, const vector &v2);
        friend vector operator-(const vector &v1, vector &&v2);
        friend vector operator-(vector &&v1, vector &&v2);

        // Subtract vectors from complex vectors
        friend cvector operator-(const vector &v1, const cvector &v2);
        friend cvector operator-(const vector &v1, cvector &&v2);
        friend cvector operator-(const cvector &v1, const vector &v2);
        friend cvector operator-(cvector &&v1, const vector &v2);

        // Compare vectors to vectors
        friend bool operator==(const vector &v1, const vector &v2);
        friend bool operator!=(const vector &v1, const vector &v2);

        // Compare vectors to complex vectors
        friend bool operator==(const vector &v1, const cvector &v2);
        friend bool operator==(const cvector &v1, const vector &v2);
        friend bool operator!=(const cvector &v1, const vector &v2);
        friend bool operator!=(const vector &v1, const cvector &v2);

    public:
        //! \brief Construct empty vector
        vector();

        //! \brief Construct zero vector of size n
        explicit vector(size_t n);

        vector(const vector &v);
        vector(vector &&v);

        vector &operator=(const vector &v);
        vector &operator=(vector &&v);

        vector &operator+=(const vector &v);

        vector &operator-=(const vector &v);

        vector &operator*=(double a);
        vector &operator/=(double a);
        vector operator-() const;
        ~vector();

        double &operator()(size_t i);
        void set(size_t i, double val);

        double operator()(size_t i) const;
        double get(size_t i) const;

        size_t size() const;

        //! \brief Access the pointer to the underlying gsl_vector
        gsl_vector *get() const { return gvec; }

        //! \brief Resize the gsl::vector, setting elements to zero
        void resize(size_t n);

        //! \brief Clear the gsl::vector, free underlying memory
        void clear();

        //! \brief Pretty-print the vector to file stream
        void print(FILE *out = stdout) const;

        //! \brief Return the 2-norm of the vector
        double norm() const { return gsl_blas_dnrm2(gvec); }

        operator vector_view() const;
        vector_view view() const;
        vector_view subvector(size_t offset, size_t n) const;

    protected:
        gsl_vector *gvec;

        //! \brief Construct new gsl::vector from gsl_vector
        vector(gsl_vector *gvec_other);

        //! \brief Private function to free allocated memory
        void gfree();

        //! \brief Private function to (continuously) allocate memory
        void galloc(size_t n);
    };

    class vector_view : public vector
    {
    public:
        //! \brief Constructor for vector_view pointing to data at gvec_other
        vector_view(gsl_vector *gvec_other);
        ~vector_view();

        // Override some nonconst member functions to be unusable
        void clear();
        void resize(size_t n);
    };

    /*! \class row_view
     *  \brief A subclass of cvector_view for stride-1 cvectors
     */
    class row_view : public vector_view
    {
    public:
        //! \brief Construct row_view from existing vector view, checking that stride is 1
        row_view(gsl_vector *gvec_other) : vector_view(gvec_other)
        {
            if (gvec_other->stride != 1)
                printf("Row view must have stride 1");
        };

        //! \brief Return a matrix view out of the elements of the row
        matrix_view reshape(size_t n, size_t m) const;
    };

    /*! \class column_view
     *  \brief A subclass of vector_view for non-stride-1 vectors
     */
    class column_view : public vector_view
    {
    public:
        //! \brief Construct column_view from existing vector view
        column_view(gsl_vector *gvec_other) : vector_view(gvec_other) {}
    };

}

#endif // YAWG_VECTOR_H_