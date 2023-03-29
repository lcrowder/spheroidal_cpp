#include "legendre_otc.h"
#include <iostream>
#include <vector>
#include <stdio.h>  
#include <math.h> 
using namespace std;


int main()
{
    int p=1;
    vector<double> u;
    for (int i = 2; i <= 4; i++)
        {u.push_back(i);}

    // vector<double> h=cont_frac(1,1,u);

    // cout << endl;
    // for (int i=0; i< 3; i++) {cout << h[i] << endl;}

    
    vector<vector<vector<double> > > PQ;
    PQ=legendre_otc(p,u,1,0,0);
    // legendre_otc(p,u,1,0,0);


    return 0;
}

