#ifndef YAWG_VECTOR_H_
#define YAWG_VECTOR_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

namespace gsl
{
    // Forward Declarations
    class cvector;
    class complex;
    class complex_ref;

    class vector
    {
    public:
        // Ordinary constructors
        vector();
        explicit vector(size_t n);

        // Copy and move constructors
        vector(const vector &gvec_other);
        vector(vector &&gvec_other);

        // Conversion constructors
        explicit vector(gsl_vector *gvec_other);
        // vector( const std::vector<double>& svec );

        //! \brief Assignment operators
        vector &operator=(const vector &gvec_other);
        vector &operator=(vector &&gvec_other);

        //! \brief Addition/subtraction assignment operators
        vector &operator+=(const vector &gvec_other);
        vector &operator-=(const vector &gvec_other);

        // Destructor
        ~vector();

        // Element access
        double &operator()(size_t i); // Setter
        void set(size_t i, double val);

        double operator()(size_t i) const; // Getter
        double get(size_t i) const;

        size_t size() const;
        double norm() const { return gsl_blas_dnrm2(gvec); }

        gsl_vector *get_gsl_ptr() const { return gvec; }

        void resize(size_t n);
        void clear();

        void print(FILE *out = stdout) const;

        friend vector operator*(double a, const vector &v);
        friend vector operator*(double a, vector &&v);
        friend vector operator*(const vector &v, double a);
        friend vector operator*(vector &&v, double a);

        // Move versions unneeded; memory can't be reclaimed
        friend cvector operator*(complex a, const vector &v);
        friend cvector operator*(const vector &v, complex a);

        friend vector operator+(const vector &v1, const vector &v2);
        friend vector operator+(vector &&v1, const vector &v2);
        friend vector operator+(const vector &v1, vector &&v2);
        friend vector operator+(vector &&v1, vector &&v2);

        friend vector operator-(const vector &v1, const vector &v2);
        friend vector operator-(vector &&v1, const vector &v2);
        friend vector operator-(const vector &v1, vector &&v2);
        friend vector operator-(vector &&v1, vector &&v2);

        friend bool operator==(const vector &v1, const vector &v2);

    protected:
        gsl_vector *gvec;
        void free();
        void calloc(size_t n);

        friend class cvector;
        friend class matrix;
    };

    class cvector
    {
    public:
        // Ordinary constructors
        cvector();
        explicit cvector(size_t n);

        // Copy and move constructors
        cvector(const cvector &gvec_other);
        cvector(cvector &&gvec_other);

        // Conversion constructors
        explicit cvector(gsl_vector_complex *gvec_other);
        cvector(const vector &vec);
        // vector( const std::vector<double>& svec );

        //! \brief Assignment operators
        cvector &operator=(const cvector &gvec_other);
        cvector &operator=(cvector &&gvec_other);

        //! \brief Addition assignment operators
        cvector &operator+=(const cvector &gvec_other);
        cvector &operator-=(const cvector &gvec_other);

        // Destructor
        ~cvector();

        // Element access (the nice C++ versions don't work)
        void set(size_t i, complex z);
        complex get(size_t i) const;

        complex_ref operator()(size_t i);             // Setter
        const complex_ref operator()(size_t i) const; // Getter

        size_t size() const;
        double norm() const { return gsl_blas_dznrm2(gvec); }

        gsl_vector_complex *get_gsl_ptr() const { return gvec; }

        void resize(size_t n);
        void clear();

        void print(FILE *out = stdout) const;

        friend cvector operator*(complex a, const cvector &v);
        friend cvector operator*(complex a, cvector &&v);
        friend cvector operator*(const cvector &v, complex a);
        friend cvector operator*(cvector &&v, complex a);

        friend cvector operator+(const cvector &v1, const cvector &v2);
        friend cvector operator+(cvector &&v1, const cvector &v2);
        friend cvector operator+(const cvector &v1, cvector &&v2);
        friend cvector operator+(cvector &&v1, cvector &&v2);

        friend cvector operator-(const cvector &v1, const cvector &v2);
        friend cvector operator-(cvector &&v1, const cvector &v2);
        friend cvector operator-(const cvector &v1, cvector &&v2);
        friend cvector operator-(cvector &&v1, cvector &&v2);

        friend bool operator==(const cvector &v1, const cvector &v2);

    protected:
        gsl_vector_complex *gvec;
        void free();
        void calloc(size_t n);

        friend class vector;
    };
}

#endif // YAWG_VECTOR_H_