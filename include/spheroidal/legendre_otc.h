#ifndef LEGENDRE_OTC_H_
#define LEGENDRE_OTC_H_

#include <yawg/vector.h>
#include <yawg/matrix.h>
using namespace std;

gsl::vector  cont_frac(int n, int m, gsl::vector  u);
int geti(int n,int m);
void legendre_otc(int p, gsl::vector u, gsl::matrix &P);
void legendre_otc(int p, gsl::vector u, gsl::matrix &P, gsl::matrix &Q);
void Dlegendre_otc(int p, gsl::vector u, gsl::matrix &P, gsl::matrix &dP);
void Dlegendre_otc(int p, gsl::vector u, gsl::matrix &P, gsl::matrix &Q, gsl::matrix &dP, gsl::matrix &dQ);

#endif // _LEGENDGRE_OTC_H_