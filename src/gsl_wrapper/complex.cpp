#include <gsl_wrapper/core.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

gsl::complex::complex(double x)
{
    GSL_SET_COMPLEX(&z, x, 0.0);
}

gsl::complex::complex(double re, double im)
{
    GSL_SET_COMPLEX(&z, re, im);
}

gsl::complex::complex(const gsl_complex &gsl_complex_other)
{
    z = gsl_complex_other;
}

//! \brief Addition assignment operator
gsl::complex &gsl::complex::operator+=(const complex &gsl_complex_other)
{
    z = gsl_complex_add(z, gsl_complex_other.z);
    return *this;
}

//! \brief Subtraction assignment operator
gsl::complex &gsl::complex::operator-=(const complex &gsl_complex_other)
{
    z = gsl_complex_sub(z, gsl_complex_other.z);
    return *this;
}

//! \brief Multiplication assignment operator
gsl::complex &gsl::complex::operator*=(const complex &gsl_complex_other)
{
    z = gsl_complex_mul(z, gsl_complex_other.z);
    return *this;
}

//! \brief Division assignment operator
gsl::complex &gsl::complex::operator/=(const complex &gsl_complex_other)
{
    z = gsl_complex_div(z, gsl_complex_other.z);
    return *this;
}

//! \brief Addition assignment operator
gsl::complex &gsl::complex::operator+=(double x)
{
    z = gsl_complex_add_real(z, x);
    return *this;
}

//! \brief Subtraction assignment operator
gsl::complex &gsl::complex::operator-=(double x)
{
    z = gsl_complex_sub_real(z, x);
    return *this;
}

//! \brief Multiplication assignment operator
gsl::complex &gsl::complex::operator*=(double x)
{
    z = gsl_complex_mul_real(z, x);
    return *this;
}

//! \brief Division assignment operator
gsl::complex &gsl::complex::operator/=(double x)
{
    z = gsl_complex_div_real(z, x);
    return *this;
}

//! \brief Assignment operator
gsl::complex &gsl::complex::operator=(double x)
{
    GSL_SET_COMPLEX(&z, x, 0.0);
    return *this;
}

//! \brief Assignment operator
gsl::complex &gsl::complex::operator=(const complex &gsl_complex_other)
{
    z = gsl_complex_other.z;
    return *this;
}

//! \brief Unary minus operator
gsl::complex gsl::complex::operator-() const
{
    return gsl_complex_negative(z);
}

//! \brief Addition operator
gsl::complex gsl::complex::operator+(const complex &gsl_complex_other) const
{
    return gsl_complex_add(z, gsl_complex_other.z);
}

//! \brief Subtraction operator
gsl::complex gsl::complex::operator-(const complex &gsl_complex_other) const
{
    return gsl_complex_sub(z, gsl_complex_other.z);
}

//! \brief Multiplication operator
gsl::complex gsl::complex::operator*(const complex &gsl_complex_other) const
{
    return gsl_complex_mul(z, gsl_complex_other.z);
}

//! \brief Division operator
gsl::complex gsl::complex::operator/(const complex &gsl_complex_other) const
{
    return gsl_complex_div(z, gsl_complex_other.z);
}