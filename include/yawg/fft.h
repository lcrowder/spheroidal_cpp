#ifndef GSL_WRAPPER_FFT_H_
#define GSL_WRAPPER_FFT_H_

#include <yawg/core.h>

namespace gsl
{
    /*! \brief Compute fft of a gsl::vector */
    gsl::cvector fft(gsl::cvector &&x);
    gsl::cvector fft(const gsl::cvector &x);

    gsl::cvector ifft(gsl::cvector &&x);
    gsl::cvector ifft(const gsl::cvector &x);

    gsl::cmatrix fft(gsl::cmatrix &&x, int dim = 1);
    gsl::cmatrix fft(const gsl::cmatrix &x, int dim = 1);

    gsl::cmatrix ifft(gsl::cmatrix &&x, int dim = 1);
    gsl::cmatrix ifft(const gsl::cmatrix &x, int dim = 1);
}

#endif // GSL_WRAPPER_FFT_H_