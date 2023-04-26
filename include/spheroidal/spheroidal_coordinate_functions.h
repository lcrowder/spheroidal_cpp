#ifndef SPHEROIDAL_COORDINATE_FUNCTIONS_H_
#define SPHEROIDAL_COORDINATE_FUNCTIONS_H_

#include <yawg/core.h>

gsl::matrix spheroidal_to_cart(gsl::matrix S, double a);
gsl::matrix cart_to_spheroidal(gsl::matrix X, double a);

#endif // _SPHEROIDAL_COORDINATE_FUNCTIONS_H_