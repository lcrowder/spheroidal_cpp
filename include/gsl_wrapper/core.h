#ifndef GSL_WRAPPER_H_
#define GSL_WRAPPER_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

namespace gsl
{
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
        vector( const cvector& gvec_other);
        // vector( const std::vector<double>& svec );

        //! \brief Assignment operators
        vector &operator=(const vector &gvec_other);
        vector &operator=(vector &&gvec_other);

        // Destructor
        ~vector();

        // Element access
        double &operator()(size_t i);      // Setter
        void set(size_t i, double val);

        double operator()(size_t i) const; // Getter
        double get(size_t i) const;


        size_t size() const;

        gsl_vector *get_gsl_ptr() const { return gvec; }

        void resize(size_t n);
        void clear();

        void print() const;

    protected:
        gsl_vector *gvec;
        void free();
        void calloc(size_t n);
    };

    class matrix
    {
    public:
        // Ordinary constructors
        matrix();
        matrix(size_t n, size_t m);

        // Conversion constructors
        explicit matrix(gsl_matrix *gmat_other);
        matrix( const cmatrix& gmat_other);
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
        double &operator()(size_t i, size_t j);      // Setters
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

        // Other functions to write later
        // vector( initialization_list )
        // Do printing better
    protected:
        gsl_matrix *gmat;
        void free();
        void calloc(size_t n, size_t m);
    };

    class complex
    {
    public:
        complex(double x);
        complex(double re, double im);
        complex(const gsl_complex &gsl_complex_other);

        inline double real() const { return GSL_REAL(z); };
        inline double imag() const { return GSL_IMAG(z); };
        inline double abs() const { return gsl_complex_abs(z); };
        inline double abs2() const { return gsl_complex_abs2(z); };
        inline double arg() const { return gsl_complex_arg(z); };

        complex &operator+=(const complex &gsl_complex_other);
        complex &operator-=(const complex &gsl_complex_other);
        complex &operator*=(const complex &gsl_complex_other);
        complex &operator/=(const complex &gsl_complex_other);

        complex &operator+=(double x);
        complex &operator-=(double x);
        complex &operator*=(double x);
        complex &operator/=(double x);

        complex &operator=(const complex &gsl_complex_other);
        complex &operator=(double x);

        complex operator-() const;

        complex operator+(const complex &gsl_complex_other) const;
        complex operator-(const complex &gsl_complex_other) const;
        complex operator*(const complex &gsl_complex_other) const;
        complex operator/(const complex &gsl_complex_other) const;

    protected:
        gsl_complex z;

        friend class cvector;
    };

    inline namespace literals
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

        // Destructor
        ~cvector();

        // Element access (the nice C++ versions don't work)
        // gsl_complex &operator()(size_t i);             // Setter
        // const gsl_complex &operator()(size_t i) const; // Getter
        void set(size_t i, const complex &z);
        complex get(size_t i) const;

        size_t size() const;

        gsl_vector_complex *get_gsl_ptr() const { return gvec; }

        void resize(size_t n);
        void clear();

        void print() const;

    protected:
        gsl_vector_complex *gvec;
        void free();
        void calloc(size_t n);
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
        void set(size_t i, size_t j, const complex &z);
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
    };

} // namespace gsl

#endif // GSL_WRAPPER_H_