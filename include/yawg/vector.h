#ifndef YAWG_VECTOR_H_
#define YAWG_VECTOR_H_

#include <yawg/complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

namespace gsl
{
    class cvector;

    /*! \class vector
     *  \brief A wrapper class for gsl_vector
     *
     * Stores and operates on a pointer to a gsl_vector.
     */
    class vector
    {
        friend class cvector;     
        friend class matrix;
        friend class cmatrix;

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
        friend bool operator!=(const vector &v1, const cvector &v2);
        friend bool operator==(const cvector &v1, const vector &v2);
        friend bool operator!=(const cvector &v1, const vector &v2);

    public:
        //! \brief Construct empty vector
        vector();

        //! \brief Construct zero vector of size n
        explicit vector(size_t n);

        //! \brief Construct new gsl::vector from gsl_vector
        vector(const gsl_vector *gvec_other);

        vector(const vector &gvec_other);
        vector(vector &&gvec_other);

        vector &operator=(const vector &gvec_other);
        vector &operator=(vector &&gvec_other);

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

    protected:
        gsl_vector *gvec;

        //! \brief Private function to free allocated memory
        void free();

        //! \brief Private function to (continuously) allocate memory
        void alloc(size_t n);
    };

}

#endif // YAWG_VECTOR_H_