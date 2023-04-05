#ifndef LEGENDRE_OTC_H_
#define LEGENDRE_OTC_H_

#include <vector>
using namespace std;

vector<double> cont_frac(int n, int m, vector<double> u);
int geti(int n,int m);
vector<vector<vector<double> > > legendre_otc(int p, vector<double> u, int Qoption = 1, int dPoption = 0, int dQoption = 0);

#endif // _LEGENDGRE_OTC_H_