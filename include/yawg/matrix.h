#ifndef YAWG_MATRIX_H_
#define YAWG_MATRIX_H_

#include <yawg/complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

namespace gsl
{
    class vector;
    class matrix;
    class cmatrix;

    /*! \class matrix
     *  \brief A wrapper class for gsl_matrix
     *
     *  Stores and operates on a pointer to a gsl_matrix.
     */
    class matrix
    {
        friend class cmatrix;

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

        //! \brief Resize the gsl::matrix, setting elements to zero
        void resize(size_t n, size_t m);

        //! \brief CLear the gsl::matrix, free underlying memory
        void clear();

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

    protected:
        gsl_matrix *gmat;

        //! \brief Private function to free allocated memory
        void free();

        //! \brief Private function to (continuously) allocate memory
        void calloc(size_t n, size_t m);
    };

} // namespace gsl

#endif // YAWG_MATRIX_H_