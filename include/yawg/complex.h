#ifndef YAWG_COMPLEX_H_
#define YAWG_COMPLEX_H_

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdio.h>

namespace gsl
{
    class complex : public gsl_complex
    {
        friend class complex_ref;

    public:
        inline complex() { GSL_SET_COMPLEX(this, 0.0, 0.0); }
        inline complex(double x) { GSL_SET_COMPLEX(this, x, 0.0); }
        inline complex(double re, double im) { GSL_SET_COMPLEX(this, re, im); }
        inline complex(gsl_complex z) { GSL_SET_COMPLEX(this, GSL_REAL(z), GSL_IMAG(z)); };

        inline double real() const { return GSL_REAL(*this); };
        inline double imag() const { return GSL_IMAG(*this); };
        inline double abs() const { return gsl_complex_abs(*this); };
        inline double abs2() const { return gsl_complex_abs2(*this); };
        inline double arg() const { return gsl_complex_arg(*this); };

        complex &operator+=(complex gsl_complex_other);
        complex &operator-=(complex gsl_complex_other);
        complex &operator*=(complex gsl_complex_other);
        complex &operator/=(complex gsl_complex_other);

        complex &operator+=(double x);
        complex &operator-=(double x);
        complex &operator*=(double x);
        complex &operator/=(double x);

        complex &operator=(complex gsl_complex_other);
        inline void set(double re, double im) { GSL_SET_COMPLEX(this, re, im); };
        inline void set(complex z) { GSL_SET_COMPLEX(this, GSL_REAL(z), GSL_IMAG(z)); };

        inline complex operator-() const { return complex(gsl_complex_negative(*this)); };

        friend inline complex operator+(complex a, complex b) { return gsl_complex_add(a, b); };
        friend inline complex operator-(complex a, complex b) { return gsl_complex_sub(a, b); };
        friend inline complex operator*(complex a, complex b) { return gsl_complex_mul(a, b); };
        friend inline complex operator/(complex a, complex b) { return gsl_complex_div(a, b); };
        friend inline bool operator==(complex a, complex b) { return (GSL_REAL(a) == GSL_REAL(b)) && (GSL_IMAG(a) == GSL_IMAG(b)); }
    
        void print() const { printf("%f + %fi", GSL_REAL(*this), GSL_IMAG(*this)); };
    };

    inline namespace complex_literals
    {
        inline complex operator""_i(long double y) { return complex(0.0, y); }
    }

    // Class to store a reference to a gsl::complex object.
    //  This is necessary for the () operator overloads of cvector and cmatrix to work.
    class complex_ref
    {
    protected:
        double *dat;
        complex_ref() : dat(nullptr){};

    public:
        // Conversion from reference to own object
        operator complex() const { return complex(dat[0], dat[1]); };

        // Constructors
        complex_ref(complex_ref& z) : dat(z.dat){};
        complex_ref(gsl_complex *z) : dat(z->dat){};
        complex_ref(complex &z) : dat(z.dat){};

        // Creates a reference to gsl::complex object z
        complex_ref &operator=(complex z)
        {
            dat[0] = z.dat[0];
            dat[1] = z.dat[1];
            return *this;
        };

        // For unclear reasons, the default assignment operator doesn't work.
        complex_ref &operator=(complex_ref z)
        {
            dat[0] = z.dat[0];
            dat[1] = z.dat[1];
            return *this;
        };

        // Allows the user to get these values without dereferencing
        inline double real() const { return dat[0]; };
        inline double imag() const { return dat[1]; };
    };

}

#endif // YAWG_COMPLEX_H_