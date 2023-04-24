#ifndef YAWG_MATRIX_H_
#define YAWG_MATRIX_H_

#include <yawg/complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

namespace gsl
{
    class vector;
    class row_view;
    class column_view;

    class cmatrix;
    class matrix_view;

    /*! \class matrix
     *  \brief A wrapper class for gsl_matrix
     *
     *  Stores and operates on a pointer to a gsl_matrix.
     */
    class matrix
    {
        // Scalar multiplication
        friend matrix operator*(double a, const matrix &M);
        friend matrix operator*(double a, matrix &&M);
        friend matrix operator*(const matrix &M, double a);
        friend matrix operator*(matrix &&M, double a);

        friend cmatrix operator*(complex a, const matrix &M);
        friend cmatrix operator*(const matrix &M, complex a);

        // Scalar division
        friend matrix operator/(double a, const matrix &M);
        friend matrix operator/(double a, matrix &&M);
        friend matrix operator/(const matrix &M, double a);
        friend matrix operator/(matrix &&M, double a);

        friend cmatrix operator/(complex z, const matrix &M);
        friend cmatrix operator/(const matrix &M, complex z);

        // Add matrices to matrices
        friend matrix operator+(const matrix &M1, const matrix &M2);
        friend matrix operator+(matrix &&M1, const matrix &M2);
        friend matrix operator+(const matrix &M1, matrix &&M2);
        friend matrix operator+(matrix &&M1, matrix &&M2);

        // Add matrices to complex matrices
        friend cmatrix operator+(const matrix &M1, const cmatrix &M2);
        friend cmatrix operator+(const matrix &M1, cmatrix &&M2);
        friend cmatrix operator+(const cmatrix &M1, const matrix &M2);
        friend cmatrix operator+(cmatrix &&M1, const matrix &M2);

        // Subtract matrices from matrices
        friend matrix operator-(const matrix &M1, const matrix &M2);
        friend matrix operator-(matrix &&M1, const matrix &M2);
        friend matrix operator-(const matrix &M1, matrix &&M2);
        friend matrix operator-(matrix &&M1, matrix &&M2);

        // Subtract matrices from complex matrices
        friend cmatrix operator-(const matrix &M1, const cmatrix &M2);
        friend cmatrix operator-(const matrix &M1, cmatrix &&M2);
        friend cmatrix operator-(const cmatrix &M1, const matrix &M2);
        friend cmatrix operator-(cmatrix &&M1, const matrix &M2);

        // Multiply matrix by matrix
        friend matrix operator*(const matrix &A, const matrix &B);

        // Compare matrix to matrix
        friend bool operator==(const matrix &M1, const matrix &M2);
        friend bool operator!=(const matrix &M1, const matrix &M2);

        // Compare matrix to complex matrix
        friend bool operator==(const matrix &M1, const cmatrix &M2);
        friend bool operator==(const cmatrix &M1, const matrix &M2);
        friend bool operator!=(const matrix &M1, const cmatrix &M2);
        friend bool operator!=(const cmatrix &M1, const matrix &M2);

    public:
        //! \brief Construct empty matrix
        matrix();

        //! \brief Construct zero matrix of size n x m
        matrix(size_t n, size_t m);

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

        matrix &operator+=(const matrix &M);
        matrix &operator-=(const matrix &M);

        matrix &operator*=(double x);
        matrix &operator/=(double x);

        matrix operator-() const;

        ~matrix();

        double &operator()(size_t i, size_t j);
        void set(size_t i, size_t j, double val);
        void set_col(size_t j, const vector &v);
        void set_row(size_t i, const vector &v);

        double operator()(size_t i, size_t j) const;
        double get(size_t i, size_t j) const;
        vector get_col(size_t j) const;
        vector get_row(size_t i) const;

        size_t size() const;
        size_t nrows() const;
        size_t ncols() const;

        bool is_square() const;

        //! \brief Access the pointer to the underlying gsl_matrix
        gsl_matrix *get() const { return gmat; }

        //! \brief CLear the gsl::matrix, free underlying memory
        void clear();

        //! \brief Resize the gsl::matrix, freeing and allocating new memory
        void resize(size_t n, size_t m);

        //! \brief Return a new n x m gsl::matrix with same elements
        matrix reshape(size_t n, size_t m) const;

        //! \brief Replace the matrix with its transpose
        matrix &T();

        //! \brief Pretty-print the matrix to file stream
        void print(FILE *out = stdout) const;

        //! \brief Print the matrix to file stream in CSV format
        void print_csv(FILE *out = stdout) const;

        //! \brief Load the matrix from a file stream in CSV format
        void load_csv(FILE *in = stdin);

        matrix_view view() const;
        matrix_view submatrix(size_t i, size_t j, size_t n, size_t m) const;

        row_view row(size_t i) const;
        column_view column(size_t j) const;

    protected:
        gsl_matrix *gmat;

        matrix(gsl_matrix *gmat);

        //! \brief Private function to free allocated memory
        void gfree();

        //! \brief Private function to (continuously) allocate memory
        void galloc(size_t n, size_t m);
    };

    class matrix_view : public matrix
    {
    public:
        //! \brief Constructor for vector_view pointing to data at gvec_other
        matrix_view(gsl_matrix *gvec_other);
        ~matrix_view();

        //! \brief Assign data from matrix to view
        matrix_view& operator=(const matrix &M);

        // Override some nonconst member functions to be unusable
        void clear();
        void resize(size_t n, size_t m);
    };

} // namespace gsl

#endif // YAWG_MATRIX_H_