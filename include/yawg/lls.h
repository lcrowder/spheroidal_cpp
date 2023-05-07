#ifndef YAWG_LLS_H_
#define YAWG_LLS_H_

#include <yawg/core.h>

namespace gsl
{
    //! \brief Compute the intercept \p c0 and slope \p c1 of best fit  
    void fit_linear(gsl::vector &x, gsl::vector& y, double &c0, double &c1, double &sumsq);
}

#endif // YAWG_LLS_H_