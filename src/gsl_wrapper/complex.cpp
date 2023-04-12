#include <gsl_wrapper/core.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

//! \brief Addition assignment operator
gsl::complex &gsl::complex::operator+=(complex gsl_complex_other)
{
    z = gsl_complex_add(z, gsl_complex_other.z);
    return *this;
}

//! \brief Subtraction assignment operator
gsl::complex &gsl::complex::operator-=(complex gsl_complex_other)
{
    z = gsl_complex_sub(z, gsl_complex_other.z);
    return *this;
}

//! \brief Multiplication assignment operator
gsl::complex &gsl::complex::operator*=(complex gsl_complex_other)
{
    z = gsl_complex_mul(z, gsl_complex_other.z);
    return *this;
}

//! \brief Division assignment operator
gsl::complex &gsl::complex::operator/=(complex gsl_complex_other)
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
gsl::complex &gsl::complex::operator=(gsl::complex gsl_complex_other)
{
    z = gsl_complex_other.z;
    return *this;
}