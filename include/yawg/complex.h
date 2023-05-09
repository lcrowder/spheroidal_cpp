#ifndef YAWG_COMPLEX_H_
#define YAWG_COMPLEX_H_

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

namespace gsl
{
    class complex_ref;

    /*! \class complex
     * \brief Wrapper class for gsl_complex structs
     *
     * Inherits `double dat[2]` from gsl_complex and provides a number of convenience functions.
     */
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

        inline complex operator-() const { return gsl_complex_negative(*this); };

        // Do math between complex and complex
        friend inline complex operator+(const complex &a, const complex &b) { return gsl_complex_add(a, b); };
        friend inline complex operator-(const complex &a, const complex &b) { return gsl_complex_sub(a, b); };
        friend inline complex operator*(const complex &a, const complex &b) { return gsl_complex_mul(a, b); };
        friend inline complex operator/(const complex &a, const complex &b) { return gsl_complex_div(a, b); };
        friend inline bool operator==(const complex &a, const complex &b) { return (a.real() == b.real()) && (a.imag() == b.imag()); };
        friend inline bool operator!=(const complex &a, const complex &b) { return (a.real() != b.real()) || (a.imag() != b.imag()); };

        // Do math between complex and double
        friend inline complex operator+(const complex &a, double b) { return gsl_complex_add_real(a, b); };
        friend inline complex operator-(const complex &a, double b) { return gsl_complex_sub_real(a, b); };
        friend inline complex operator*(const complex &a, double b) { return gsl_complex_mul_real(a, b); };
        friend inline complex operator/(const complex &a, double b) { return gsl_complex_div_real(a, b); };
        friend inline bool operator==(const complex &a, double b) { return (a.real() == b) && (a.imag() == 0.0); };
        friend inline bool operator!=(const complex &a, double b) { return (a.imag() != 0.0) || (a.real() != b); };

        // Do math between double and complex
        friend inline complex operator+(double a, const complex &b) { return gsl_complex_add_real(b, a); };
        friend inline complex operator-(double a, const complex &b) { return gsl_complex_negative(gsl_complex_sub_real(b, a)); };
        friend inline complex operator*(double a, const complex &b) { return gsl_complex_mul_real(b, a); };
        friend inline complex operator/(double a, const complex &b) { return gsl_complex_mul_real(gsl_complex_inverse(b), a); };
        friend inline bool operator==(double a, const complex &b) { return (b.real() == a) && (b.imag() == 0.0); };
        friend inline bool operator!=(double a, const complex &b) { return (b.imag() != 0.0) || (b.real() != a); };

        void print() const;
    };

    inline namespace complex_literals
    {
        //! \brief User defined literal overload for complex numbers
        inline complex operator""_i(long double y) { return complex(0.0, y); }
    }

    /*! \class complex_ref
     * \brief Stores a refernce to a gsl::complex object
     *
     * This class is necessary to communicate between gsl_complex and gsl::complex
     * so that overloads of () work gsl::cvector and gsl::cmatrix.
     *
     * Implementation heavily inspired by ccgsl (https://ccgsl.sourceforge.net/)
     */
    class complex_ref
    {
        friend class complex;

        // Do math between complex_ref and complex_ref
        friend inline complex operator+(const complex_ref &a, const complex_ref &b) { return gsl_complex_add(a, b); };
        friend inline complex operator-(const complex_ref &a, const complex_ref &b) { return gsl_complex_sub(a, b); };
        friend inline complex operator*(const complex_ref &a, const complex_ref &b) { return gsl_complex_mul(a, b); };
        friend inline complex operator/(const complex_ref &a, const complex_ref &b) { return gsl_complex_div(a, b); };
        friend inline bool operator==(const complex_ref &a, const complex_ref &b) { return (a.real() == b.real()) && (a.imag() == b.imag()); };
        friend inline bool operator!=(const complex_ref &a, const complex_ref &b) { return (a.real() != b.real()) || (a.imag() != b.imag()); };

        // Do math between complex and double
        friend inline complex operator+(const complex_ref &a, double b) { return gsl_complex_add_real(a, b); };
        friend inline complex operator-(const complex_ref &a, double b) { return gsl_complex_sub_real(a, b); };
        friend inline complex operator*(const complex_ref &a, double b) { return gsl_complex_mul_real(a, b); };
        friend inline complex operator/(const complex_ref &a, double b) { return gsl_complex_div_real(a, b); };
        friend inline bool operator==(const complex_ref &a, double b) { return (a.real() == b) && (a.imag() == 0.0); };
        friend inline bool operator!=(const complex_ref &a, double b) { return (a.imag() != 0.0) || (a.real() != b); };

        // Do math between double and complex ref
        friend inline complex operator+(double a, const complex_ref &b) { return gsl_complex_add_real(b, a); };
        friend inline complex operator-(double a, const complex_ref &b) { return gsl_complex_sub_real(b, a); };
        friend inline complex operator*(double a, const complex_ref &b) { return gsl_complex_mul_real(b, a); };
        friend inline complex operator/(double a, const complex_ref &b) { return gsl_complex_div_real(b, a); };
        friend inline bool operator==(double a, const complex_ref &b) { return (b.real() == a) && (b.imag() == 0.0); };
        friend inline bool operator!=(double a, const complex_ref &b) { return (b.imag() != 0.0) || (b.real() != a); };

    protected:
        double *dat;
        complex_ref() : dat(nullptr){};

    public:
        //! \brief "Dereferences" a complex_ref into independent gsl::complex object
        operator complex() const { return complex(dat[0], dat[1]); };

        //! \brief Constructs a reference to another compelx_ref object
        complex_ref(complex_ref &z) : dat(z.dat){};

        //! \brief Constructs a reference to a gsl_complex struct
        complex_ref(gsl_complex *z) : dat(z->dat){};

        //! \brief Constructs a reference to a gsl::complex object
        complex_ref(complex &z) : dat(z.dat){};

        complex_ref &operator+=(const complex &gsl_complex_other);
        complex_ref &operator-=(const complex &gsl_complex_other);
        complex_ref &operator*=(const complex &gsl_complex_other);
        complex_ref &operator/=(const complex &gsl_complex_other);

        complex_ref &operator+=(double x);
        complex_ref &operator-=(double x);
        complex_ref &operator*=(double x);
        complex_ref &operator/=(double x);

        //! \brief Assigns values of gsl::complex object to the reference
        complex_ref &operator=(const complex& z)
        {
            dat[0] = z.dat[0];
            dat[1] = z.dat[1];
            return *this;
        };

        //! \brief Assigns values of one complex_ref object to another.
        //! \note This is an alternative to the default assignment operator,
        //! which does not work for unknown reasons.
        complex_ref &operator=(const complex_ref& z)
        {
            dat[0] = z.dat[0];
            dat[1] = z.dat[1];
            return *this;
        };

        inline double real() const { return dat[0]; };
        inline double imag() const { return dat[1]; };
    };

}

#endif // YAWG_COMPLEX_H_