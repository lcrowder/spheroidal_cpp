#ifndef GSL_WRAPPER_CORE_H_
#define GSL_WRAPPER_CORE_H_

#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

namespace gsl
{
    // Forward declarations
    class vector;
    class matrix;
    class complex;
    class cvector;
    class cmatrix;

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

        void print() const;

        friend vector operator*( double a, const vector &v );
        friend vector operator*( double a, vector &&v );
        friend vector operator*( const vector &v, double a );
        friend vector operator*( vector &&v, double a );

        // Move versions unneeded; memory can't be reclaimed
        friend cvector operator*( complex a, const vector &v );
        friend cvector operator*( const vector &v, complex a );

        friend vector operator+( const vector &v1, const vector &v2 );
        friend vector operator+( vector &&v1, const vector &v2 );
        friend vector operator+( const vector &v1, vector &&v2 );
        friend vector operator+( vector &&v1, vector &&v2 );

        friend vector operator-( const vector &v1, const vector &v2 );
        friend vector operator-( vector &&v1, const vector &v2 );
        friend vector operator-( const vector &v1, vector &&v2 );
        friend vector operator-( vector &&v1, vector &&v2 );

        friend bool operator==( const vector &v1, const vector &v2 );
    protected:
        gsl_vector *gvec;
        void free();
        void calloc(size_t n);

        friend class cvector;
    };

    class matrix
    {
    public:
        // Ordinary constructors
        matrix();
        matrix(size_t n, size_t m);

        // Conversion constructors
        explicit matrix(gsl_matrix *gmat_other);
        // matrix( const std::vector<std::vector<double>>& smat );

        // Copy and move constructors
        matrix(const matrix &gmat_other);
        matrix(matrix &&gmat_other);

        //! \brief Assignment operators
        matrix &operator=(const matrix &gmat_other);
        matrix &operator=(matrix &&gmat_other);

        //! \brief Destructor
        ~matrix();

        //! \brief Element access
        double &operator()(size_t i, size_t j); // Setters
        void set(size_t i, size_t j, double val);

        double operator()(size_t i, size_t j) const; // Getters
        double get(size_t i, size_t j) const;

        size_t size() const;
        size_t nrows() const;
        size_t ncols() const;

        gsl_matrix *get_gsl_ptr() const { return gmat; }

        void resize(size_t n, size_t m);
        void clear();

        void print() const;
        void print2csv(FILE* out = stdout) const;
    protected:
        gsl_matrix *gmat;
        void free();
        void calloc(size_t n, size_t m);

        friend class cmatrix;
    };

    class complex
    {
    public:
        inline complex() { GSL_SET_COMPLEX(&z, 0.0, 0.0); }
        inline complex(double x) { GSL_SET_COMPLEX(&z, x, 0.0); }
        inline complex(double re, double im) { GSL_SET_COMPLEX(&z, re, im); }
        inline complex(gsl_complex gsl_complex_other) { z = gsl_complex_other; };

        inline double real() const { return GSL_REAL(z); };
        inline double imag() const { return GSL_IMAG(z); };
        inline double abs() const { return gsl_complex_abs(z); };
        inline double abs2() const { return gsl_complex_abs2(z); };
        inline double arg() const { return gsl_complex_arg(z); };

        complex &operator+=(complex gsl_complex_other);
        complex &operator-=(complex gsl_complex_other);
        complex &operator*=(complex gsl_complex_other);
        complex &operator/=(complex gsl_complex_other);

        complex &operator+=(double x);
        complex &operator-=(double x);
        complex &operator*=(double x);
        complex &operator/=(double x);

        complex &operator=(complex gsl_complex_other);
        inline void set(double re, double im) { GSL_SET_COMPLEX(&z, re, im); };
        inline void set(complex x) { z = x.z; };

        inline complex operator-() const { return gsl_complex_negative(z); };

        inline gsl_complex get_gsl_data() const { return z; };

        friend inline complex operator+(complex a, complex b) { return gsl_complex_add(a.z, b.z); };
        friend inline complex operator-(complex a, complex b) { return gsl_complex_sub(a.z, b.z); };
        friend inline complex operator*(complex a, complex b) { return gsl_complex_mul(a.z, b.z); };
        friend inline complex operator/(complex a, complex b) { return gsl_complex_div(a.z, b.z); };
        friend inline bool operator==(complex a, complex b) { return (GSL_REAL(a.z) == GSL_REAL(b.z)) && (GSL_IMAG(a.z) == GSL_IMAG(b.z)); }

    protected:
        gsl_complex z;
    };

    inline namespace complex_literals
    {
        inline complex operator""_i(long double y) { return complex(0.0, y); }
    }

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
        // gsl_complex &operator()(size_t i);             // Setter
        // const gsl_complex &operator()(size_t i) const; // Getter
        void set(size_t i, complex z);
        complex get(size_t i) const;

        size_t size() const;
        double norm() const { return gsl_blas_dznrm2(gvec); }

        gsl_vector_complex *get_gsl_ptr() const { return gvec; }

        void resize(size_t n);
        void clear();

        void print() const;

        friend cvector operator*( complex a, const cvector &v );
        friend cvector operator*( complex a, cvector &&v );
        friend cvector operator*( const cvector &v, complex a );
        friend cvector operator*( cvector &&v, complex a );

        friend cvector operator+( const cvector &v1, const cvector &v2 );
        friend cvector operator+( cvector &&v1, const cvector &v2 );
        friend cvector operator+( const cvector &v1, cvector &&v2 );
        friend cvector operator+( cvector &&v1, cvector &&v2 );

        friend cvector operator-( const cvector &v1, const cvector &v2 );
        friend cvector operator-( cvector &&v1, const cvector &v2 );
        friend cvector operator-( const cvector &v1, cvector &&v2 );
        friend cvector operator-( cvector &&v1, cvector &&v2 );

        friend bool operator==( const cvector &v1, const cvector &v2 );
    protected:
        gsl_vector_complex *gvec;
        void free();
        void calloc(size_t n);

        friend class vector;
    };

    class cmatrix
    {
    public:
        // Ordinary constructors
        cmatrix();
        explicit cmatrix(size_t n, size_t m);

        // Copy and move constructors
        cmatrix(const cmatrix &gvec_other);
        cmatrix(cmatrix &&gvec_other);

        // Conversion constructors
        explicit cmatrix(gsl_matrix_complex *gvec_other);
        cmatrix(const matrix &mat);
        // vector( const std::vector<double>& svec );

        //! \brief Assignment operators
        cmatrix &operator=(const cmatrix &gvec_other);
        cmatrix &operator=(cmatrix &&gvec_other);

        // Destructor
        ~cmatrix();

        // Element access
        void set(size_t i, size_t j, complex z);
        complex get(size_t i, size_t j) const;
        // complex &operator()(size_t i, size_t j);             // Getter
        // complex operator()(size_t i, size_t j) const; // Setter

        size_t size() const;
        size_t nrows() const;
        size_t ncols() const;

        gsl_matrix_complex *get_gsl_ptr() const { return gmat; }

        void resize(size_t n, size_t m);
        void clear();

        void print() const;

    protected:
        gsl_matrix_complex *gmat;
        void free();
        void calloc(size_t n, size_t m);

        friend class matrix;
    };

} // namespace gsl

#endif // GSL_WRAPPER_CORE_H_