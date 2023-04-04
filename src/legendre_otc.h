#ifndef LEGENDRE_OTC_H
#define LEGENDRE_OTC_H

#include <vector>

std::vector<double> cont_frac(int n, int m, std::vector<double> u);
int geti( int n, int m );
int double_factorial( int n );
int factorial( int n );
std::vector<std::vector<std::vector<double> > > legendre_otc(int p,std::vector<double> u, int n, int m, int flag);

#endif