#include <yawg/core.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

gsl::complex &gsl::complex::operator+=(complex z)
{
    *this = gsl_complex_add(*this, z);
    return *this;
}

gsl::complex &gsl::complex::operator-=(complex z)
{
    *this = gsl_complex_sub(*this, z);
    return *this;
}

gsl::complex &gsl::complex::operator*=(complex z)
{
    *this = gsl_complex_mul(*this, z);
    return *this;
}

gsl::complex &gsl::complex::operator/=(complex z)
{
    *this = gsl_complex_div(*this, z);
    return *this;
}

gsl::complex &gsl::complex::operator+=(double x)
{
    *this = gsl_complex_add_real(*this, x);
    return *this;
}

gsl::complex &gsl::complex::operator-=(double x)
{
    *this = gsl_complex_sub_real(*this, x);
    return *this;
}

gsl::complex &gsl::complex::operator*=(double x)
{
    *this = gsl_complex_mul_real(*this, x);
    return *this;
}

gsl::complex &gsl::complex::operator/=(double x)
{
    *this = gsl_complex_div_real(*this, x);
    return *this;
}

gsl::complex &gsl::complex::operator=(gsl::complex z)
{
    dat[0] = z.dat[0];
    dat[1] = z.dat[1];
    return *this;
}

void gsl::complex::print( ) const
{
    printf("%f + %fi", GSL_REAL(*this), GSL_IMAG(*this));
}