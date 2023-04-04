#ifndef LEGENDRE_OTC_H
#define LEGENDRE_OTC_H

vector<double> cont_frac(int n, int m, vector<double> u);
int geti( int n, int m, int p );
int double_factorial( int n );
int factorial( int n );
vector<vector<vector<double> > > legendre_otc(int p, vector<double> u, int n, int m, int flag);

#endif