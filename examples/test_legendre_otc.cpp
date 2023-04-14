#include <spheroidal/legendre_otc.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>  
#include <math.h> 
#include <gsl_wrapper/utils.hpp>
#include <gsl_wrapper/core.h>
using namespace std;


int main(int argc, char *argv[])
{
    int p=2, N=3;
    double ua=1.1, ub=2.1;
    
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

    // double h=(ub-ua)/(N-1);
    // vector<double> u;
    // // for (double ui = ua; ui <= ub; ui+=h)
    // //     {u.push_back(ui);}

    // double uj=ua;
    // for (int j=0; j<N; ++j)
    // {
    //     u.push_back(ua+j*h);
    //     // cout << "u[" << j << "]=" << uj << endl;
    // }

    gsl::vector u = gsl::linspace(ua,ub,N);

    cout << "length of u is " << u.size() << endl;
    cout << "Last u-value is " << u(u.size()-1) << endl;

    gsl::matrix P, dP;
    gsl::matrix Q, dQ;

    // legendre_otc(p, u, P, Q);
    // legendre_otc(p, u, P);
    Dlegendre_otc(p, u, P, Q, dP, dQ);

    cout << "Last entry of P is " << P.get((p+1)*(p+1)-1,N-1) << endl;
    cout << "Last entry of Q is " << Q.get((p+1)*(p+1)-1,N-1) << endl;
    cout << "Last entry of dP is " << dP.get((p+1)*(p+1)-1,N-1) << endl;
    cout << "Last entry of dQ is " << dQ.get((p+1)*(p+1)-1,N-1) << endl;

    int sp=(p+1)*(p+1);

    
    // Open a file for writing
    ofstream Pfile("../../data/P.txt");
    ofstream Qfile("../../data/Q.txt");
    ofstream dPfile("../../data/dP.txt");
    ofstream dQfile("../../data/dQ.txt");

    // Set precision to 16 decimal places
    Pfile << fixed << setprecision(16);
    Qfile << fixed << setprecision(16);
    dPfile << fixed << setprecision(16);
    dQfile << fixed << setprecision(16);

    // Check if files are open
    if (!Pfile.is_open()) {
        cerr << "Unable to open P file" << std::endl;
        return 1;
    }
    if (!Qfile.is_open()) {
        cerr << "Unable to open Q file" << std::endl;
        return 1;
    }
    if (!dPfile.is_open()) {
        cerr << "Unable to open dP file" << std::endl;
        return 1;
    }
    if (!dQfile.is_open()) {
        cerr << "Unable to open dQ file" << std::endl;
        return 1;
    }

    // Write the vector to the files
    for (int i=0; i<sp; ++i) {
        for (int j=0; j<N; ++j) {
            Pfile << scientific << P.get(i,j) << " ";
            Qfile << scientific << Q.get(i,j) << " ";
            dPfile << scientific << dP.get(i,j) << " ";
            dQfile << scientific << dQ.get(i,j) << " ";
        }
        Pfile << endl;
        Qfile << endl;
        dPfile << endl;
        dQfile << endl;
    }

    // Close the files
    Pfile.close();
    Qfile.close();
    dPfile.close();
    dQfile.close();

    return 0;
}

