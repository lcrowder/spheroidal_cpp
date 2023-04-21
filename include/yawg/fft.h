#ifndef YAWG_FFT_H_
#define YAWG_FFT_H_

#include <yawg/core.h>

namespace gsl
{
    //! \brief Compute in-place fft of a gsl::(c)vector
    gsl::cvector fft(gsl::cvector &&x);

    //! \brief Compute fft of a gsl::(c)vector
    gsl::cvector fft(const gsl::cvector &x);

    //! \brief Compute in-place inverse fft of a gsl::(c)vector
    gsl::cvector ifft(gsl::cvector &&x);

    //! \brief Compute inverse fft of a gsl::(c)vector
    gsl::cvector ifft(const gsl::cvector &x);

    //! \brief Compute in-place 1D fft of each column/row of a gsl::(c)matrix
    gsl::cmatrix fft(gsl::cmatrix &&x, int dim = 1);

    //! \brief Compute 1D fft of each column/row of a gsl::(c)matrix
    gsl::cmatrix fft(const gsl::cmatrix &x, int dim = 1);

    //! \brief Compute in-place 1D inverse fft of each column/row of a gsl::(c)matrix
    gsl::cmatrix ifft(gsl::cmatrix &&x, int dim = 1);

    //! \brief Compute 1D inverse fft of each column/row of a gsl::(c)matrix
    gsl::cmatrix ifft(const gsl::cmatrix &x, int dim = 1);
}

#endif // YAWG_FFT_H_