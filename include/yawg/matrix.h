#ifndef YAWG_MATRIX_H_
#define YAWG_MATRIX_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

namespace gsl
{
    class vector;
    class row_view;
    class column_view;

    class cvector_view;
    class crow_view;
    class ccolumn_view;

    class matrix;
    class matrix_view;

    class cmatrix;
    class cmatrix_view;
    
    /*! \class matrix
     *  \brief A wrapper class for gsl_matrix
     *
     *  Stores and operates on a pointer to a gsl_matrix.
     */
    class matrix
    {
        friend class cmatrix;

        friend class matrix_view;
        friend class row_view;
        friend class column_view;

    public:
        //! \brief Construct empty matrix
        matrix();
        //! \brief Construct zero matrix of size n x m
        matrix(size_t n, size_t m);

        //! \brief Construct new gsl::matrix from gsl_matrix
        matrix(const gsl_matrix *gsl_mat);

        //! \brief Construct new gsl::matrix from .csv file
        matrix(FILE *in);

        //! \brief Construct new n x 1 gsl::matrix from gsl::vector
        matrix(const vector &v);

        //! \brief Copy constructor creating n x m matrix
        matrix(const matrix &M, size_t n, size_t m);

        matrix(const matrix &M);
        matrix(matrix &&M);

        matrix &operator=(const matrix &M);
        matrix &operator=(matrix &&M);

        ~matrix();

        double &operator()(size_t i, size_t j);
        void set(size_t i, size_t j, double val);

        double operator()(size_t i, size_t j) const;
        double get(size_t i, size_t j) const;

        size_t size() const;
        size_t nrows() const;
        size_t ncols() const;

        //! \brief Access the pointer to the underlying gsl_matrix
        gsl_matrix *get_gsl_ptr() const { return gmat; }

        //! \brief Resize the gsl::matrix, setting elements to zero
        void resize(size_t n, size_t m);

        //! \brief CLear the gsl::matrix, free underlying memory
        void clear();

        //! \brief Return a new n x m gsl::matrix with same elements
        matrix reshape(size_t n, size_t m) const;

        //! \brief Pretty-print the matrix to file stream
        void print(FILE *out = stdout) const;

        //! \brief Print the matrix to file stream in CSV format
        void print_csv(FILE *out = stdout) const;

        //! \brief Load the matrix from a file stream in CSV format
        void load_csv(FILE *in = stdin);

        friend matrix operator*(const matrix &A, const matrix &B);

        //! \brief Return a view to a submatrix of the matrix
        matrix_view submatrix(size_t i, size_t j, size_t n, size_t m);

        //! \brief Return a view to a row of the matrix
        row_view row(size_t i);

        //! \brief Return a view to a column of the matrix
        column_view column(size_t j);

    protected:
        gsl_matrix *gmat;

        //! \brief Private function to free allocated memory
        void free();

        //! \brief Private function to (continuously) allocate memory
        void calloc(size_t n, size_t m);
    };

    /*! \class matrix_view
     *  \brief A wrapper class for gsl_matrix_view
     *
     * Stores a gsl_matrix_view and uses it to access original member data.
     */
    class matrix_view
    {
    protected:
        gsl_matrix_view gmat_view;

    public:
        //! \brief "Dereferences" a matrix_view into independent gsl::matrix object
        operator matrix() const { return matrix(&gmat_view.matrix); }

        //! \brief Construct a view of a gsl::matrix through another matrix_view
        matrix_view(gsl_matrix_view gmat_view) : gmat_view(gmat_view){};

        //! \brief Construct a view of the given gsl::matrix
        matrix_view(const matrix &m) : gmat_view(gsl_matrix_submatrix(m.gmat, 0, 0, m.gmat->size1, m.gmat->size2)) {}

        //! \brief Assignment to a matrix view from a matrix
        matrix_view &operator=(const matrix &M);

        //! \brief Assignment to a matrix view from another matrix view
        matrix_view &operator=(matrix_view Mv);

        //! \brief Pretty-print the viewed matrix to file stream
        void print(FILE *out = stdout) const;

        //! \brief Return a constant pointer to the underlying gsl_matrix
        const gsl_matrix *get_gsl_ptr() const { return &gmat_view.matrix; }
    };

} // namespace gsl

#endif // YAWG_MATRIX_H_