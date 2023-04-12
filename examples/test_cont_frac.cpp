#include <spheroidal/legendre_otc.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>  
#include <math.h> 
#include <gsl_wrapper/core.h>
#include <gsl_wrapper/utils.hpp>
using namespace std;

int main()
{
    int n=1;
    int m=0;
    gsl::vector u;
    u(0)=(1.5);

    gsl::vector cf = cont_frac(n,m,u);
    cout << "cf=" << cf(0) << endl;
}