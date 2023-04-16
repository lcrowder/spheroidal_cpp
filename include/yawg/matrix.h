#ifndef YAWG_MATRIX_H_
#define YAWG_MATRIX_H_

#include <gsl/gsl_matrix.h>
#include <stdio.h>

namespace gsl
{
    // Forward Declarations
    class vector;
    class complex;
    class complex_ref;
    
    class matrix
    {
    public:
        // Ordinary constructors
        matrix();
        matrix(size_t n, size_t m);

        // Conversion constructors
        explicit matrix(gsl_matrix *gmat_other);
        matrix(const vector &v);
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
        matrix reshape(size_t n, size_t m);
        void clear();

        void print(FILE *out = stdout) const;
        void print_csv(FILE *out = stdout) const;
        void load_csv(FILE *in = stdin);

    protected:
        gsl_matrix *gmat;
        void free();
        void calloc(size_t n, size_t m);

        friend class cmatrix;
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

        complex_ref operator()(size_t i, size_t j);             // Setter
        const complex_ref operator()(size_t i, size_t j) const; // Getter

        size_t size() const;
        size_t nrows() const;
        size_t ncols() const;

        gsl_matrix_complex *get_gsl_ptr() const { return gmat; }

        void resize(size_t n, size_t m);
        void reshape(size_t n, size_t m);

        void clear();

        void print(FILE *out = stdout) const;

    protected:
        gsl_matrix_complex *gmat;
        void free();
        void calloc(size_t n, size_t m);

        friend class matrix;
    };

} // namespace gsl

#endif // YAWG_MATRIX_H_