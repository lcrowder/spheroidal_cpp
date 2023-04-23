#ifndef YAWG_CMATRIX_H_
#define YAWG_CMATRIX_H_

#include <yawg/complex.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>

namespace gsl
{
    class cvector;
    class matrix;

    class cmatrix_view;
    class crow_view;
    class ccolumn_view;
    class cmatrix
    {
        // Scalar multiplication
        friend cmatrix operator*(double a, const cmatrix &M);
        friend cmatrix operator*(double a, cmatrix &&M);
        friend cmatrix operator*(const cmatrix &M, double a);
        friend cmatrix operator*(cmatrix &&M, double a);

        friend cmatrix operator*(complex z, const cmatrix &M);
        friend cmatrix operator*(const cmatrix &M, complex z);
        friend cmatrix operator*(cmatrix &&M, complex z);
        friend cmatrix operator*(complex z, cmatrix &&M);

        // Scalar multiplication
        friend cmatrix operator/(double a, const cmatrix &M);
        friend cmatrix operator/(double a, cmatrix &&M);
        friend cmatrix operator/(const cmatrix &M, double a);
        friend cmatrix operator/(cmatrix &&M, double a);

        friend cmatrix operator/(complex z, const cmatrix &M);
        friend cmatrix operator/(const cmatrix &M, complex z);
        friend cmatrix operator/(cmatrix &&M, complex z);
        friend cmatrix operator/(complex z, cmatrix &&M);

        // Add complex matrices to complex matrices
        friend cmatrix operator+(const cmatrix &M1, const cmatrix &M2);
        friend cmatrix operator+(cmatrix &&M1, const cmatrix &M2);
        friend cmatrix operator+(const cmatrix &M1, cmatrix &&M2);
        friend cmatrix operator+(cmatrix &&M1, cmatrix &&M2);

        // Add complex matrices to matrices
        friend cmatrix operator+(const cmatrix &M1, const matrix &M2);
        friend cmatrix operator+(cmatrix &&M1, const matrix &M2);
        friend cmatrix operator+(const matrix &M1, const cmatrix &M2);
        friend cmatrix operator+(const matrix &M1, cmatrix &&M2);

        // Subtract complex matrices and complex matrices
        friend cmatrix operator-(const cmatrix &M1, const cmatrix &M2);
        friend cmatrix operator-(cmatrix &&M1, const cmatrix &M2);
        friend cmatrix operator-(const cmatrix &M1, cmatrix &&M2);
        friend cmatrix operator-(cmatrix &&M1, cmatrix &&M2);

        // Subtract complex matrices from matrices
        friend cmatrix operator-(const cmatrix &M1, const matrix &M2);
        friend cmatrix operator-(cmatrix &&M1, const matrix &M2);
        friend cmatrix operator-(const matrix &M1, const cmatrix &M2);
        friend cmatrix operator-(const matrix &M1, cmatrix &&M2);

        // Multiply complex matrices
        friend cmatrix operator*(const cmatrix &A, const cmatrix &B);

        // Compare complex matrices to complex matrices
        friend bool operator==(const cmatrix &M1, const cmatrix &M2);
        friend bool operator!=(const cmatrix &M1, const cmatrix &M2);

        // Compare complex matrices to matrices
        friend bool operator==(const cmatrix &M1, const matrix &M2);
        friend bool operator==(const matrix &M1, const cmatrix &M2);

        friend bool operator!=(const cmatrix &M1, const matrix &M2);
        friend bool operator!=(const matrix &M1, const cmatrix &M2);

    public:
        //! \brief Construct empty matrix
        cmatrix();
        
        //! \brief Construct zero matrix of size n x m
        cmatrix(size_t n, size_t m);

        //! \brief Construct new gsl::cmatrix from MATLAB's .csv file format
        cmatrix(FILE *in);

        //! \brief Construct new n x 1 gsl::cmatrix from a gsl::cvector
        cmatrix(const cvector &v);

        //! \brief Copy constructor creating n x m complex matrix
        cmatrix(const cmatrix &M, size_t n, size_t m);

        cmatrix(const cmatrix &M);
        cmatrix(cmatrix &&M);

        //! \brief Construct new gsl::cmatrix from gsl::matrix
        cmatrix(const matrix &M);

        cmatrix &operator=(const cmatrix &M);
        cmatrix &operator=(const matrix &M);
        cmatrix &operator=(cmatrix &&v);

        cmatrix &operator+=(const cmatrix &M);
        cmatrix &operator+=(const matrix &M);

        cmatrix &operator-=(const cmatrix &M);
        cmatrix &operator-=(const matrix &M);

        cmatrix &operator*=(complex z);
        cmatrix &operator*=(double x);

        cmatrix &operator/=(complex z);
        cmatrix &operator/=(double x);

        cmatrix operator-() const;

        ~cmatrix();

        //! \brief Return a reference to the element at position (i,j)
        complex_ref operator()(size_t i, size_t j);
        void set(size_t i, size_t j, complex z);
        void set_col(size_t j, const cvector &v);
        void set_row(size_t i, const cvector &v);

        //! \brief Return a const reference to the element at position (i,j)
        const complex_ref operator()(size_t i, size_t j) const;
        complex get(size_t i, size_t j) const;
        cvector get_col(size_t j) const;
        cvector get_row(size_t i) const;

        size_t size() const;
        size_t nrows() const;
        size_t ncols() const;

        bool is_square() const;

        //! \brief Access the pointer to the underlying gsl_matrix_complex
        gsl_matrix_complex *get() const { return gmat; }

        //! \brief Clear the gsl::cmatrix, free underlying memory
        void clear();

        //! \brief Resize the gsl::cmatrix
        void resize(size_t n, size_t m);

        //! \brief Return a new n x m gsl::cmatrix with same elements
        cmatrix reshape(size_t n, size_t m) const;

        //! \brief Replace the complex matrix with its transpose
        cmatrix &T();

        //! \brief Replace the complex matrix with its conjugate transpose
        cmatrix &H();

        //! \brief Return the conjugate of the complex matrix
        cmatrix &conj();

        //! \brief Pretty-print the complex matrix to file stream
        void print(FILE *out = stdout) const;

        //! \brief Print the complex matrix to file stream in MATLAB's .csv format
        void print_csv(FILE *out = stdout) const;

        //! \brief Load the complex matrix from file stream in MATLAB's .csv format
        void load_csv(FILE *in);

        friend cmatrix operator*(const cmatrix &A, const cmatrix &B);

        operator cmatrix_view() const;
        cmatrix_view view() const;
        cmatrix_view submatrix(size_t i, size_t j, size_t n, size_t m) const;

        crow_view row(size_t i) const;
        ccolumn_view column(size_t j) const;
    protected:
        gsl_matrix_complex *gmat;

        //! \brief Construct new gsl::matrix from gsl_matrix
        cmatrix(gsl_matrix_complex *gmat_other);

        //! \brief Private function to free allocated memory
        void gfree();

        //! \brief Private function to (continuously) allocate memory
        void galloc(size_t n, size_t m);
    };

    class cmatrix_view : public cmatrix
    {
    public:
        //! \brief Constructor for vector_view pointing to data at gvec_other
        cmatrix_view(gsl_matrix_complex *gvec_other);
        ~cmatrix_view();

        // Override some nonconst member functions to be unusable
        void clear();
        void resize(size_t n, size_t m);
    };

} // namespace gsl
#endif // YAWG_CVECTOR_H