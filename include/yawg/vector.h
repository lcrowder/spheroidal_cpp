#ifndef YAWG_VECTOR_H_
#define YAWG_VECTOR_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

namespace gsl
{
    // Forward Declarations
    class vector_view;
    class row_view;
    class column_view;

    class cvector_view;
    class crow_view;
    class ccolumn_view;

    class matrix_view;
    class cmatrix_view;

    class cvector;

    class complex;
    class complex_ref;

    class vector
    {
        friend class vector_view;
        friend class row_view;
        friend class cvector;
        friend class matrix;

    public:
        // Ordinary constructors
        vector();
        explicit vector(size_t n);

        // Copy and move constructors
        vector(const vector &gvec_other);
        vector(vector &&gvec_other);

        // Conversion constructors
        vector(const gsl_vector *gvec_other);
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

        void print(FILE *out = stdout) const;

        friend vector operator*(double a, const vector &v);
        friend vector operator*(double a, vector &&v);
        friend vector operator*(const vector &v, double a);
        friend vector operator*(vector &&v, double a);

        // Move versions unneeded; memory can't be reclaimed
        friend cvector operator*(complex a, const vector &v);
        friend cvector operator*(const vector &v, complex a);

        friend vector operator+(const vector &v1, const vector &v2);
        friend vector operator+(vector &&v1, const vector &v2);
        friend vector operator+(const vector &v1, vector &&v2);
        friend vector operator+(vector &&v1, vector &&v2);

        friend vector operator-(const vector &v1, const vector &v2);
        friend vector operator-(vector &&v1, const vector &v2);
        friend vector operator-(const vector &v1, vector &&v2);
        friend vector operator-(vector &&v1, vector &&v2);

        friend bool operator==(const vector &v1, const vector &v2);

        vector_view subvector(size_t offset, size_t size);

    protected:
        gsl_vector *gvec;
        void free();
        void calloc(size_t n);
    };

    class vector_view
    {
    protected:
        gsl_vector_view gvec_view;

    public:
        // Create a new object from a vector_view
        operator vector() const { return vector(&gvec_view.vector); }

        // Constructors
        vector_view(gsl_vector_view gvec_view) : gvec_view(gvec_view){};
        vector_view(const vector &v) : gvec_view(gsl_vector_subvector(v.gvec, 0, v.gvec->size)){};

        vector_view &operator=(const vector &v);
        vector_view &operator=(vector_view v);

        void print(FILE *out = stdout) const;
    };

    // Because C++ is row-major, row view have stride one, and can be reshaped into matrices
    class row_view : public vector_view
    {
    public:
        row_view(gsl_vector_view gvec_view) : vector_view(gvec_view)
        {
            if (gvec_view.vector.stride != 1)
                printf("Row view must have stride 1");
        };

        row_view(const vector &v) : vector_view(v)
        {
            if (v.gvec->stride != 1)
                printf("Row view must have stride 1");
        };

        matrix_view reshape(size_t n, size_t m);
    };

    class column_view : public vector_view
    {
    public:
        column_view(gsl_vector_view gvec_view) : vector_view(gvec_view){};
        column_view(const vector &v) : vector_view(v){};
    };

    class cvector
    {
        friend class vector;
        friend class cvector_view;
        friend class crow_view;
        friend class ccolumn_view;

    public:
        // Ordinary constructors
        cvector();
        explicit cvector(size_t n);

        // Copy and move constructors
        cvector(const cvector &gvec_other);
        cvector(cvector &&gvec_other);

        // Conversion constructors
        cvector(const gsl_vector_complex *gvec_other);
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
        void set(size_t i, complex z);
        complex get(size_t i) const;

        complex_ref operator()(size_t i);             // Setter
        const complex_ref operator()(size_t i) const; // Getter

        size_t size() const;
        double norm() const { return gsl_blas_dznrm2(gvec); }

        gsl_vector_complex *get_gsl_ptr() const { return gvec; }

        void resize(size_t n);
        void clear();

        void print(FILE *out = stdout) const;

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

    protected:
        gsl_vector_complex *gvec;
        void free();
        void calloc(size_t n);
    };

    class cvector_view
    {
    protected:
        gsl_vector_complex_view gvec_view;

    public:
        // Create a new object from a vector_view
        operator cvector() const { return cvector(&gvec_view.vector); }

        // Constructors
        cvector_view(gsl_vector_complex_view gvec_view) : gvec_view(gvec_view){};
        cvector_view(const cvector &v) : gvec_view(gsl_vector_complex_subvector(v.gvec, 0, v.gvec->size)){};

        cvector_view &operator=(const cvector &v);
        cvector_view &operator=(cvector_view v);

        // Print method. Inefficient, but printing shouldn't be production ready
        void print(FILE *out = stdout) const;
    };

    // Because C++ is row-major, row view have stride one, and can be reshaped into matrices
    class crow_view : public cvector_view
    {
    public:
        crow_view(gsl_vector_complex_view gvec_view) : cvector_view(gvec_view)
        {
            if (gvec_view.vector.stride != 1)
                printf("Row view must have stride 1");
        };

        crow_view(const cvector &v) : cvector_view(v)
        {
            if (v.gvec->stride != 1)
                printf("Row view must have stride 1");
        };

        cmatrix_view reshape(size_t n, size_t m);
    };

    class ccolumn_view : public cvector_view
    {
    public:
        ccolumn_view(gsl_vector_complex_view gvec_view) : cvector_view(gvec_view){};
        ccolumn_view(const cvector &v) : cvector_view(v){};
    };

}

#endif // YAWG_VECTOR_H_