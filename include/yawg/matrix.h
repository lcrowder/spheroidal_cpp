#ifndef YAWG_MATRIX_H_
#define YAWG_MATRIX_H_

#include <gsl/gsl_matrix.h>
#include <stdio.h>

namespace gsl
{
    // Forward Declarations
    class matrix;
    class vector;
    class complex;
    class complex_ref;

    class matrix_view;
    class row_view;
    class column_view;

    class cmatrix_view;
    class crow_view;
    class ccolumn_view;


    class matrix
    {
        friend class matrix_view;

        friend class vector_view;
        friend class row_view;
        friend class column_view;
     
        friend class cmatrix;

    public:
        // Ordinary constructors
        matrix();
        matrix(size_t n, size_t m);

        // Conversion constructors
        matrix(const gsl_matrix *gmat_other);
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
        void clear();

        matrix reshape(size_t n, size_t m);
        matrix(const matrix &gmat_other, size_t n, size_t m);

        void print(FILE *out = stdout) const;
        void print_csv(FILE *out = stdout) const;
        void load_csv(FILE *in = stdin);

        // Matrix multiplication can't be done in place
        friend matrix operator*(const matrix &A, const matrix &B);

        matrix_view submatrix( size_t i, size_t j, size_t n, size_t m );
        row_view row( size_t i );
        column_view column( size_t j );

    protected:
        gsl_matrix *gmat;
        void free();
        void calloc(size_t n, size_t m);
    };

    class matrix_view
    {
    protected:
        gsl_matrix_view gmat_view;

    public:
        // Create a new object from a vector_view
        operator matrix() const { return matrix(&gmat_view.matrix); }

        // Constructors
        matrix_view(gsl_matrix_view gmat_view) : gmat_view(gmat_view){};
        matrix_view(const matrix &m) : gmat_view(gsl_matrix_submatrix( m.gmat, 0, 0, m.gmat->size1, m.gmat->size2)) {}

        matrix_view &operator=(const matrix &v);
        matrix_view &operator=(matrix_view v);

        void print(FILE *out = stdout) const;
    };

    class cmatrix
    {
        friend class cmatrix_view;
        friend class matrix;

    public:
        // Ordinary constructors
        cmatrix();
        explicit cmatrix(size_t n, size_t m);

        // Copy and move constructors
        cmatrix(const cmatrix &gvec_other);
        cmatrix(cmatrix &&gvec_other);

        // Conversion constructors
        cmatrix(const gsl_matrix_complex *gvec_other);
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

        friend cmatrix operator*(const cmatrix &A, const cmatrix &B);

        cmatrix_view submatrix( size_t i, size_t j, size_t n, size_t m );
        crow_view row( size_t i );
        ccolumn_view column( size_t j );

    protected:
        gsl_matrix_complex *gmat;
        void free();
        void calloc(size_t n, size_t m);
    };

    class cmatrix_view
    {
    protected:
        gsl_matrix_complex_view gmat_view;

    public:
        // Create a new object from a vector_view
        operator cmatrix() const { return cmatrix(&gmat_view.matrix); }

        // Constructors
        cmatrix_view(gsl_matrix_complex_view gmat_view) : gmat_view(gmat_view){};
        cmatrix_view(const cmatrix &m) : gmat_view(gsl_matrix_complex_submatrix( m.gmat, 0, 0, m.gmat->size1, m.gmat->size2)) {}

        cmatrix_view &operator=(const cmatrix &v);
        cmatrix_view &operator=(cmatrix_view v);

        void print(FILE *out = stdout) const;
    };

} // namespace gsl

#endif // YAWG_MATRIX_H_