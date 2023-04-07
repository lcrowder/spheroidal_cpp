#ifndef GSL_WRAPPER_H_
#define GSL_WRAPPER_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <iostream>

namespace gsl
{

    class vector
    {
    public:
        // Ordinary constructors
        vector();
        vector(size_t n);

        // Copy constructors
        vector(gsl_vector *gvec_other);
        vector(const vector &gvec_other);

        // Destructor
        ~vector();

        // Element access
        double &operator()(size_t i);             // Getter
        const double &operator()(size_t i) const; // Setter

        //! \brief Assignment operator
        vector &operator=(const vector &gvec_other);

        size_t size() const;

        gsl_vector *get_gsl_ptr() const { return gvec; }

        void resize(size_t n);
        void clear();

        void print() const;

        // Other functions to write later
        // vector( initialization_list )
        // Fix printing lol
    protected:
        gsl_vector *gvec;
        void free();
        void calloc(size_t n);
    };

    class matrix
    {
    public:
        // Ordinary constructors
        matrix();
        matrix(size_t n, size_t m);

        // Copy constructors
        matrix(gsl_matrix *gmat_other);
        matrix(const matrix &gmat_other);

        // Destructor
        ~matrix();

        // Element access
        double &operator()(size_t i, size_t j);             // Getter
        const double &operator()(size_t i, size_t j) const; // Setter

        //! \brief Assignment operator
        matrix &operator=(const matrix &gmat_other);

        size_t size() const;
        size_t nrows() const;
        size_t ncols() const;

        gsl_matrix *get_gsl_ptr() const { return gmat; }

        void resize(size_t n, size_t m);
        void clear();

        void print() const;

        // Other functions to write later
        // vector( initialization_list )
        // Do printing better
    protected:
        gsl_matrix *gmat;
        void free();
        void calloc(size_t n, size_t m);
    };

}

#endif // GSL_WRAPPER_H_