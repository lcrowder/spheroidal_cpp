#include <spheroidal/legendre_otc.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>  
#include <math.h> 
using namespace std;

int main()
{
    int n=1;
    int m=0;
    vector<double> u;
    u.push_back(1.5);

    vector<double> cf = cont_frac(n,m,u);
    cout << "cf=" << cf[0] << endl;
}