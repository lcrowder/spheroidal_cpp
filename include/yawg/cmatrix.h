#ifndef YAWG_CMATRIX_H_
#define YAWG_CMATRIX_H_

#include <yawg/complex.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>

namespace gsl
{
    class cmatrix_view; 
    class crow_view;
    class ccolumn_view;

    class cvector;
    class matrix;

    class cmatrix
    {
        friend class matrix;

        friend class cmatrix_view;
        friend class crow_view;
        friend class ccolumn_view;

    public:
        //! \brief Construct empty matrix
        cmatrix();
        //! \brief Construct zero matrix of size n x m
        cmatrix(size_t n, size_t m);

        //! \brief Construct new gsl::matrix from gsl_matrix
        cmatrix(const gsl_matrix_complex *gmat_other);

        //! \brief Construct new n x 1 gsl::cmatrix from a gsl::cvector
        cmatrix(const cvector &v);

        //! \brief Copy constructor creating n x m complex matrix
        cmatrix(const cmatrix &M, size_t n, size_t m);

        cmatrix(const cmatrix &M);
        cmatrix(cmatrix &&M);

        //! \brief Construct new gsl::cmatrix from gsl::matrix
        cmatrix(const matrix &M);

        cmatrix &operator=(const cmatrix &gvec_other);
        cmatrix &operator=(cmatrix &&gvec_other);

        ~cmatrix();

        //! \brief Return a reference to the element at position (i,j)
        complex_ref operator()(size_t i, size_t j);
        void set(size_t i, size_t j, complex z);

        //! \brief Return a const reference to the element at position (i,j)
        complex get(size_t i, size_t j) const;
        const complex_ref operator()(size_t i, size_t j) const;

        size_t size() const;
        size_t nrows() const;
        size_t ncols() const;

        //! \brief Access the pointer to the underlying gsl_matrix_complex
        gsl_matrix_complex *get_gsl_ptr() const { return gmat; }

        //! \brief Resize the gsl::cmatrix, setting elements to zero
        void resize(size_t n, size_t m);

        //! \brief CLear the gsl::cmatrix, free underlying memory
        void clear();

        //! \brief Return a new n x m gsl::cmatrix with same elements
        cmatrix reshape(size_t n, size_t m) const;

        //! \brief Pretty-print the complex matrix to file stream
        void print(FILE *out = stdout) const;

        friend cmatrix operator*(const cmatrix &A, const cmatrix &B);

        //! \brief Return a view to a submatrix of the complex matrix
        cmatrix_view submatrix(size_t i, size_t j, size_t n, size_t m);

        //! \brief Return a view to a row of the complex matrix
        crow_view row(size_t i);

        //! \brief Return a view to a column of the complex matrix
        ccolumn_view column(size_t j);

    protected:
        gsl_matrix_complex *gmat;

        //! \brief Private function to free allocated memory
        void free();

        //! \brief Private function to (continuously) allocate memory
        void calloc(size_t n, size_t m);
    };

    /*! \class cmatrix_view
     *  \brief A wrapper class for gsl_matrix_complex_view
     *
     * Stores a gsl_matrix_complex_view and uses it to access original member data.
     */
    class cmatrix_view
    {
    protected:
        gsl_matrix_complex_view gmat_view;

    public:
        //! \brief "Dereferences" a matrix_view into independent gsl::cmatrix object
        operator cmatrix() const { return cmatrix(&gmat_view.matrix); }

        //! \brief Construct a view of a gsl::cmatrix through another cmatrix_view
        cmatrix_view(gsl_matrix_complex_view gmat_view) : gmat_view(gmat_view){};

        //! \brief Construct a view of the given gsl::cmatrix
        cmatrix_view(const cmatrix &m) : gmat_view(gsl_matrix_complex_submatrix(m.gmat, 0, 0, m.gmat->size1, m.gmat->size2)) {}

        //! \brief Assignment to a complex matrix view from a complex matrix
        cmatrix_view &operator=(const cmatrix &v);

        //! \brief Assignment to a complex matrix view from another complex matrix view
        cmatrix_view &operator=(cmatrix_view v);

        //! \brief Pretty-print the viewed complex matrix to file stream
        void print(FILE *out = stdout) const;

        //! \brief Return a constant pointer to the underlying gsl_matrix_complex
        const gsl_matrix_complex *get_gsl_ptr() const { return &gmat_view.matrix; }
    };

} // namespace gsl
#endif // YAWG_CVECTOR_H