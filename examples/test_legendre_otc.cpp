#include <spheroidal/legendre_otc.h>
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
    int p=2, N=3;
    double ua=1.1, ub=2.1;
    int sp=(p+1)*(p+1);
    
    if (argc==5)
    {
        p=atoi(argv[1]);
        ua=atof(argv[2]);
        ub=atof(argv[3]);
        N=atoi(argv[4]);
    }
    

    cout << "p=" << p <<endl;
    cout << "ua=" << ua <<endl;
    cout << "ub=" << ub <<endl;
    cout << "N=" << N <<endl;

    gsl::vector u = gsl::linspace(ua,ub,N);

    cout << "length of u is " << u.size() << endl;
    cout << "Last u-value is " << u(u.size()-1) << endl;

    gsl::matrix P, dP;
    gsl::matrix Q, dQ;

    legendre_otc(p, u, P);
    // legendre_otc(p, u, P, Q);
    // Dlegendre_otc(p, u, P, dP);
    // Dlegendre_otc(p, u, P, Q, dP, dQ);

    P.print();

    gsl::vector Prow = P.row(0);
    Prow=2.*Prow;
    Prow.print();
    P.row(0)=Prow;
    P.print();

    cout << "Last entry of P is " << P.get((p+1)*(p+1)-1,N-1) << endl;
    cout << "Last entry of Q is " << Q.get((p+1)*(p+1)-1,N-1) << endl;
    cout << "Last entry of dP is " << dP.get((p+1)*(p+1)-1,N-1) << endl;
    cout << "Last entry of dQ is " << dQ.get((p+1)*(p+1)-1,N-1) << endl;

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

