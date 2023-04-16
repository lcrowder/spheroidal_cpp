#include <yawg/core.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

//! \brief Addition assignment operator
gsl::complex &gsl::complex::operator+=(complex z)
{
    *this = gsl_complex_add(*this, z);
    return *this;
}

//! \brief Subtraction assignment operator
gsl::complex &gsl::complex::operator-=(complex z)
{
    *this = gsl_complex_sub(*this, z);
    return *this;
}

//! \brief Multiplication assignment operator
gsl::complex &gsl::complex::operator*=(complex z)
{
    *this = gsl_complex_mul(*this, z);
    return *this;
}

//! \brief Division assignment operator
gsl::complex &gsl::complex::operator/=(complex z)
{
    *this = gsl_complex_div(*this, z);
    return *this;
}

//! \brief Addition assignment operator
gsl::complex &gsl::complex::operator+=(double x)
{
    *this = gsl_complex_add_real(*this, x);
    return *this;
}

//! \brief Subtraction assignment operator
gsl::complex &gsl::complex::operator-=(double x)
{
    *this = gsl_complex_sub_real(*this, x);
    return *this;
}

//! \brief Multiplication assignment operator
gsl::complex &gsl::complex::operator*=(double x)
{
    *this = gsl_complex_mul_real(*this, x);
    return *this;
}

//! \brief Division assignment operator
gsl::complex &gsl::complex::operator/=(double x)
{
    *this = gsl_complex_div_real(*this, x);
    return *this;
}

//! \brief Assignment operator
gsl::complex &gsl::complex::operator=(gsl::complex z)
{
    *this = z;
    return *this;
}