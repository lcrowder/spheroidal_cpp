#ifndef YAWG_CVECTOR_H_
#define YAWG_CVECTOR_H_

#include <yawg/complex.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

namespace gsl
{
    class vector;

    class cvector_view;
    class crow_view;
    class ccolumn_view;

    class cmatrix_view;

    class vector;

    /*! \class cvector
     *  \brief A wrapper class for gsl_vector_complex
     *
     * Stores and operates on a pointer to a gsl_vector_complex.
     */
    class cvector
    {
        friend class vector;

        friend class cvector_view;
        friend class crow_view;
        friend class ccolumn_view;

        friend class cmatrix;

    public:
        //! \brief Construct empty vector
        cvector();
        //! \brief Construct zero vector of size n
        explicit cvector(size_t n);

        //! \brief Construct new gsl::cvector from gsl_vector_complex
        cvector(const gsl_vector_complex *gvec_other);

        cvector(const cvector &gvec_other);
        cvector(cvector &&gvec_other);

        //! \brief Construct new gsl::cvector from gsl::vector
        cvector(const vector &vec);

        cvector &operator=(const cvector &gvec_other);
        cvector &operator=(cvector &&gvec_other);

        ~cvector();

        //! \brief Return a reference to the element at position (i,j)
        complex_ref operator()(size_t i);
        void set(size_t i, complex z);

        const complex_ref operator()(size_t i) const;
        complex get(size_t i) const;

        size_t size() const;

        //! \brief Access the pointer to the underlying gsl_vector_complex
        gsl_vector_complex *get_gsl_ptr() const { return gvec; }

        //! \brief Resize the gsl::cvector, setting elements to zero
        void resize(size_t n);

        //! \brief Clear the gsl::cvector, free underlying memory
        void clear();

        //! \brief Pretty-print the complex vector to file stream
        void print(FILE *out = stdout) const;

        //! \brief Return the 2-norm of the vector
        double norm() const { return gsl_blas_dznrm2(gvec); }

        cvector &operator+=(const cvector &gvec_other);
        cvector &operator-=(const cvector &gvec_other);

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

        //! \brief Return a view to a subvector of the vector
        cvector_view subvector(size_t offset, size_t size);

        //! \brief Return a view to the entire vector
        cvector_view view();

    protected:
        gsl_vector_complex *gvec;

        //! \brief Private function to free allocated memory
        void free();

        //! \brief Private function to (continuously) allocate memory
        void calloc(size_t n);
    };

    /*! \class cvector_view
     *  \brief A wrapper class for gsl_vector_complex_view
     *
     * Stores a gsl_vector_complex_view and uses it to access original member data.
     */
    class cvector_view
    {
    protected:
        gsl_vector_complex_view gvec_view;

    public:
        //! \brief "Dereferences" a cvector_view into independent gsl::cvector object
        operator cvector() const { return cvector(&gvec_view.vector); }

        //! \brief Construct a view of a gsl::cvector through another cvector_view
        cvector_view(gsl_vector_complex_view gvec_view) : gvec_view(gvec_view){};

        //! \brief Construct a view of the given gsl::cvector
        cvector_view(const cvector &v) : gvec_view(gsl_vector_complex_subvector(v.gvec, 0, v.gvec->size)){};

        //! \brief Assignment to a cvector view from a cvector
        cvector_view &operator=(const cvector &v);

        //! \brief Assignment to a cvector view from another cvector view
        cvector_view &operator=(cvector_view v);

        //! \brief Pretty-print the viewed vector to file stream
        void print(FILE *out = stdout) const;

        //! \brief Return a constant pointer to the underlying gsl_vector_complex
        const gsl_vector_complex *get_gsl_ptr() const { return &gvec_view.vector; };
    };

    /*! \class crow_view
     *  \brief A subclass of cvector_view for stride-1 vectors
     */
    class crow_view : public cvector_view
    {
    public:
        //! \brief Construct crow_view from existing cvector view, checking that stride is 1
        crow_view(gsl_vector_complex_view gvec_view) : cvector_view(gvec_view)
        {
            if (gvec_view.vector.stride != 1)
                printf("Row view must have stride 1");
        };

        //! \brief Construct crow_view from cvector, checking that stride is 1
        crow_view(const cvector &v) : cvector_view(v)
        {
            if (v.gvec->stride != 1)
                printf("Row view must have stride 1");
        };

        //! \brief Return a cmatrix view out of the elements of the row
        cmatrix_view reshape(size_t n, size_t m);
    };

    /*! \class ccolumn_view
     *  \brief A subclass of cvector_view for non-stride-1 cvectors
     */
    class ccolumn_view : public cvector_view
    {
    public:
        //! \brief Construct ccolumn_view from existing cvector view
        ccolumn_view(gsl_vector_complex_view gvec_view) : cvector_view(gvec_view){};

        //! \brief Construct ccolumn_view from cvector
        ccolumn_view(const cvector &v) : cvector_view(v){};
    };

} // namespace gsl

#endif // YAWG_CVECTOR_H_