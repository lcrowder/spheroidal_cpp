#include <yawg/core.h>
#include <yawg/lls.h>
#include <gsl/gsl_fit.h>

/*! \brief Compute in-place fft of a gsl::(c)vector
 * Computes the in-place complex Fast Fourier Transform of the data in x.
 * \param x The cvector to be transformed.
 *
 * \note More efficient implementations of this function (and the others in this file)
 * would use gsl_fft_real_* functions for real arguments, or cache the wavetable
 * and workspace for repeated calls in different parts of execution.
 *
 * \return The transformed vector, reclaiming the memory of the original vector.
 */
void gsl::fit_linear(gsl::vector &x, gsl::vector& y, double &c0, double &c1, double &sumsq)
{
    gsl_vector *xv = x.get();
    gsl_vector *yv = y.get();

    double cov00, cov01, cov11;

    gsl_fit_linear(xv->data, xv->stride, yv->data, yv->stride, xv->size, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
}