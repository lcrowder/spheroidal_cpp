#ifndef YAWG_VECTOR_H_
#define YAWG_VECTOR_H_

#include <yawg/complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

namespace gsl
{
    class cvector;

    class vector_view;
    class row_view;
    class column_view;

    class cvector_view;

    class matrix_view;

    /*! \class vector
     *  \brief A wrapper class for gsl_vector
     *
     * Stores and operates on a pointer to a gsl_vector.
     */
    class vector
    {
        friend class cvector;

        friend class vector_view;
        friend class row_view;
        friend class column_view;

        friend class cvector_view;
        
        friend class matrix;

        // Scalar multiplication
        friend vector operator*(double a, const vector &v);
        friend vector operator*(double a, vector &&v);
        friend vector operator*(const vector &v, double a);
        friend vector operator*(vector &&v, double a);

        friend cvector operator*(complex a, const vector &v);
        friend cvector operator*(const vector &v, complex a);

        // Scalar division
        friend vector operator/(double a, const vector &v);
        friend vector operator/(double a, vector &&v);
        friend vector operator/(const vector &v, double a);
        friend vector operator/(vector &&v, double a);

        friend cvector operator/(complex z, const vector &v);
        friend cvector operator/(const vector &v, complex z);

        // Add vectors to vectors
        friend vector operator+(const vector &v1, const vector &v2);
        friend vector operator+(vector &&v1, const vector &v2);
        friend vector operator+(const vector &v1, vector &&v2);
        friend vector operator+(vector &&v1, vector &&v2);

        // Add vector views to vectors
        friend vector operator+(const vector &v1, const vector_view &v2);
        friend vector operator+(const vector_view &v1, const vector &v2);
        friend vector operator+(vector &&v1, const vector_view &v2);
        friend vector operator+(const vector_view &v1, vector &&v2);

        // Add vectors to complex vectors
        friend cvector operator+(const vector &v1, const cvector &v2);
        friend cvector operator+(const vector &v1, cvector &&v2);
        friend cvector operator+(const cvector &v1, const vector &v2);
        friend cvector operator+(cvector &&v1, const vector &v2);

        // Add vectors to complex vector views
        friend cvector operator+(const vector &v1, const cvector_view &v2);
        friend cvector operator+(const cvector_view &v1, const vector &v2);

        // Subtract vectors from vectors
        friend vector operator-(const vector &v1, const vector &v2);
        friend vector operator-(vector &&v1, const vector &v2);
        friend vector operator-(const vector &v1, vector &&v2);
        friend vector operator-(vector &&v1, vector &&v2);

        // Subtract vector views from vectors
        friend vector operator-(const vector &v1, const vector_view &v2);
        friend vector operator-(const vector_view &v1, const vector &v2);
        friend vector operator-(vector &&v1, const vector_view &v2);
        friend vector operator-(const vector_view &v1, vector &&v2);

        // Subtract vectors from complex vectors
        friend cvector operator-(const vector &v1, const cvector &v2);
        friend cvector operator-(const vector &v1, cvector &&v2);
        friend cvector operator-(const cvector &v1, const vector &v2);
        friend cvector operator-(cvector &&v1, const vector &v2);

        // Subtract vectors from complex vector views
        friend cvector operator-(const vector &v1, const cvector_view &v2);
        friend cvector operator-(const cvector_view &v1, const vector &v2);

        // Compare vectors to vectors
        friend bool operator==(const vector &v1, const vector &v2);
        friend bool operator!=(const vector &v1, const vector &v2);

        // Compare vectors to vector views
        friend bool operator==(const vector_view &v1, const vector &v2);
        friend bool operator==(const vector &v1, const vector_view &v2);
        friend bool operator!=(const vector_view &v1, const vector &v2);
        friend bool operator!=(const vector &v1, const vector_view &v2);

        // Compare vectors to complex vectors
        friend bool operator==(const vector &v1, const cvector &v2);
        friend bool operator!=(const vector &v1, const cvector &v2);
        friend bool operator==(const cvector &v1, const vector &v2);
        friend bool operator!=(const cvector &v1, const vector &v2);

        // Compare vectors to complex vector views
        friend bool operator==(const vector &v1, const cvector_view &v2);
        friend bool operator!=(const vector &v1, const cvector_view &v2);
        friend bool operator==(const cvector_view &v1, const vector &v2);
        friend bool operator!=(const cvector_view &v1, const vector &v2);

    public:
        //! \brief Construct empty vector
        vector();

        //! \brief Construct zero vector of size n
        explicit vector(size_t n);

        //! \brief Construct new gsl::vector from gsl_vector
        vector(const gsl_vector *gvec_other);

        vector(const vector &gvec_other);
        vector(vector &&gvec_other);

        vector &operator=(const vector &gvec_other);
        vector &operator=(const vector_view &gvec_other);
        vector &operator=(vector &&gvec_other);

        vector &operator+=(const vector &v);
        vector &operator+=(const vector_view &vv);

        vector &operator-=(const vector &v);
        vector &operator-=(const vector_view &vv);

        vector &operator*=(double a);
        vector &operator/=(double a);
        vector operator-() const;
        ~vector();

        double &operator()(size_t i);
        void set(size_t i, double val);

        double operator()(size_t i) const;
        double get(size_t i) const;

        size_t size() const;

        //! \brief Access the pointer to the underlying gsl_vector
        gsl_vector *get_gsl_ptr() const { return gvec; }

        //! \brief Resize the gsl::vector, setting elements to zero
        void resize(size_t n);

        //! \brief Clear the gsl::vector, free underlying memory
        void clear();

        //! \brief Pretty-print the vector to file stream
        void print(FILE *out = stdout) const;

        //! \brief Return the 2-norm of the vector
        double norm() const { return gsl_blas_dnrm2(gvec); }



        //! \brief Return a view to a subvector of the vector
        vector_view subvector(size_t offset, size_t size);

        //! \brief Return a view to the entire vector
        vector_view view();

    protected:
        gsl_vector *gvec;

        //! \brief Private function to free allocated memory
        void free();

        //! \brief Private function to (continuously) allocate memory
        void calloc(size_t n);
    };

    /*! \class vector_view
     *  \brief A wrapper class for gsl_vector_view
     *
     * Stores a gsl_vector_view and uses it to access original member data.
     */
    class vector_view
    {
        friend class vector;
        friend class cvector;
        friend class cvector_view;

        // Scalar multiplication
        friend vector operator*(const vector_view &v, double a);
        friend vector operator*(double a, const vector_view &v);

        // Scalar division
        friend vector operator/(const vector_view &v, double a);
        friend vector operator/(double a, const vector_view &v);
        
        // Add vector views to vector views
        friend vector operator+(const vector_view &v1, const vector_view &v2);

        // Add vector view to complex vector view
        friend cvector operator+(const vector_view &v1, const cvector_view &v2);
        friend cvector operator+(const cvector_view &v1, const vector_view &v2);

        // Add vector views to vectors
        friend vector operator+(const vector &v1, const vector_view &v2);
        friend vector operator+(const vector_view &v1, const vector &v2);
        friend vector operator+(vector &&v1, const vector_view &v2);
        friend vector operator+(const vector_view &v1, vector &&v2);

        // Add vector views to complex vectors
        friend cvector operator+(const vector_view &v1, const cvector &v2);
        friend cvector operator+(const vector_view &v1, cvector &&v2);
        friend cvector operator+(const cvector &v1, const vector_view &v2);
        friend cvector operator+(cvector &&v1, const vector_view &v2);

        // Subtract vector views from vector views
        friend vector operator-(const vector_view &v1, const vector_view &v2);

        // Subtract vector view from complex vector view
        friend cvector operator-(const vector_view &v1, const cvector_view &v2);
        friend cvector operator-(const cvector_view &v1, const vector_view &v2);

        // Subtract vector views from vectors
        friend vector operator-(const vector &v1, const vector_view &v2);
        friend vector operator-(const vector_view &v1, const vector &v2);
        friend vector operator-(vector &&v1, const vector_view &v2);
        friend vector operator-(const vector_view &v1, vector &&v2);

        // Subtract vector views from complex vectors
        friend cvector operator-(const vector_view &v1, const cvector &v2);
        friend cvector operator-(const vector_view &v1, cvector &&v2);
        friend cvector operator-(const cvector &v1, const vector_view &v2);
        friend cvector operator-(cvector &&v1, const vector_view &v2);

        // Compare vector views to vector views
        friend bool operator==(const vector_view &v1, const vector_view &v2);
        friend bool operator!=(const vector_view &v1, const vector_view &v2);

        // Compare vector views to complex vector views
        friend bool operator==(const vector_view &v1, const cvector_view &v2);
        friend bool operator!=(const vector_view &v1, const cvector_view &v2);
        friend bool operator==(const cvector_view &v1, const vector_view &v2);
        friend bool operator!=(const cvector_view &v1, const vector_view &v2);

        // Compare vector views to vectors
        friend bool operator==(const vector_view &v1, const vector &v2);
        friend bool operator!=(const vector_view &v1, const vector &v2);
        friend bool operator==(const vector &v1, const vector_view &v2);
        friend bool operator!=(const vector &v1, const vector_view &v2);

        // Compare vector views to complex vectors
        friend bool operator==(const vector_view &v1, const cvector &v2);
        friend bool operator!=(const vector_view &v1, const cvector &v2);
        friend bool operator==(const cvector &v1, const vector_view &v2);
        friend bool operator!=(const cvector &v1, const vector_view &v2);

    protected:
        gsl_vector_view gvec_view;

    public:
        //! \brief "Dereferences" a vector_view into independent gsl::vector object
        operator vector() const { return vector(&gvec_view.vector); }

        //! \brief Construct a view of a gsl::vector through another vector_view
        vector_view(gsl_vector_view gvec_view) : gvec_view(gvec_view){};

        //! \brief Construct a view of the given gsl::vector
        vector_view(const vector &v) : gvec_view(gsl_vector_subvector(v.gvec, 0, v.gvec->size)){};

        //! \brief Assignment to a vector view from a vector
        vector_view &operator=(const vector &v);

        //! \brief Assignment to a vector view from a vector view
        vector_view &operator=(vector_view v);

        vector_view &operator+=(const vector_view &vv);
        vector_view &operator+=(const vector &v);

        vector_view &operator-=(const vector_view &vv);
        vector_view &operator-=(const vector &v);
        
        vector_view &operator*=(double x);
        vector_view &operator/=(double x);
        
        double &operator()(size_t i);
        void set(size_t i, double val);

        double operator()(size_t i) const;
        double get(size_t i) const;

        size_t size() const { return gvec_view.vector.size; };

        //! \brief Pretty-print the viewed vector to file stream
        void print(FILE *out = stdout) const;

        //! \brief Return a constant pointer to the underlying gsl_vector
        const gsl_vector *get_gsl_ptr() const { return &gvec_view.vector; };
    };

    /*! \class row_view
     *  \brief A subclass of vector_view for stride-1 vectors
     */
    class row_view : public vector_view
    {
    public:
        //! \brief Construct row_view from existing vector view, checking that stride is 1
        row_view(gsl_vector_view gvec_view) : vector_view(gvec_view)
        {
            if (gvec_view.vector.stride != 1)
                printf("Row view must have stride 1");
        };

        //! \brief Construct row_view from vector, checking that stride is 1
        row_view(const vector &v) : vector_view(v)
        {
            if (v.gvec->stride != 1)
                printf("Row view must have stride 1");
        };

        //! \brief Return a matrix view out of the elements of the row
        matrix_view reshape(size_t n, size_t m);
    };

    /*! \class column_view
     *  \brief A subclass of vector_view for non-stride-1 vectors
     */
    class column_view : public vector_view
    {
    public:
        //! \brief Construct column_view from existing vector view
        column_view(gsl_vector_view gvec_view) : vector_view(gvec_view){};

        //! \brief Construct column_view from vector
        column_view(const vector &v) : vector_view(v){};
    };
}

#endif // YAWG_VECTOR_H_