#include <spheroidal/spheroidal_analysis.h>
#include <spheroidal/grid_functions.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>  
#include <math.h> 
#include <yawg/utils.hpp>
#include <yawg/core.h>
using namespace std;


int main(int argc, char *argv[])
{
    int p=2;
    int sp=(p+1)*(p+1);
    
    if (argc==2){p=atoi(argv[1]);}
    
    cout << "p=" << p <<endl;

    gsl::matrix theta, phi;
    gl_grid(p, theta, phi);
    theta.print();
    phi.print();

    gsl::matrix f;
    f.resize(2*p*(p+1),1);
    f.column(0) = gsl::linspace(1. , 2.*p*(p+1), 2*p*(p+1));

    gsl::matrix shc=spheroidal_analysis(f);
 

    // cout << "Last entry of P is " << P.get((p+1)*(p+1)-1,N-1) << endl;
    // cout << "Last entry of Q is " << Q.get((p+1)*(p+1)-1,N-1) << endl;
    // cout << "Last entry of dP is " << dP.get((p+1)*(p+1)-1,N-1) << endl;
    // cout << "Last entry of dQ is " << dQ.get((p+1)*(p+1)-1,N-1) << endl;

    // FILE* Pfile = fopen("../../data/P.csv","w");
    // FILE* Qfile = fopen("../../data/Q.csv","w");
    // FILE* dPfile = fopen("../../data/dP.csv","w");
    // FILE* dQfile = fopen("../../data/dQ.csv","w");

    // P.print_csv(Pfile); 
    // Q.print_csv(Qfile); 
    // dP.print_csv(dPfile); 
    // dQ.print_csv(dQfile);

    // fclose(Pfile);
    // fclose(Qfile);
    // fclose(dPfile);
    // fclose(dQfile);

    return 0;
}

