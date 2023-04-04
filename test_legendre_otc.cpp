#include "legendre_otc.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>  
#include <math.h> 
using namespace std;


int main(int argc, char *argv[])
{
    int p=atoi(argv[1]);
    double ua=atof(argv[2]);
    double ub=atof(argv[3]);
    int N=atoi(argv[4]);

    double h=(ub-ua)/(N-1);
    vector<double> u;
    // for (double ui = ua; ui <= ub; ui+=h)
    //     {u.push_back(ui);}

    double uj=ua;
    for (int j=0; j<N; ++j)
    {
        u.push_back(ua+j*h);
        // cout << "u[" << j << "]=" << uj << endl;
    }

    cout << "length of u is " << u.size() << endl;
    cout << "Last u-value is " << u[u.size()-1] << endl;
    
    int sp=(p+1)*(p+1);

    // vector<double> h=cont_frac(1,1,u);

    // cout << endl;
    // for (int i=0; i< 3; i++) {cout << h[i] << endl;}

    
    vector<vector<vector<double> > > PQdPdQ;
    PQdPdQ=legendre_otc(p,u,1,1,1);

    vector<vector<double> > P = PQdPdQ[0];
    vector<vector<double> > Q = PQdPdQ[1];
    vector<vector<double> > dP = PQdPdQ[2];
    vector<vector<double> > dQ = PQdPdQ[3];

    // Open a file for writing
    ofstream Pfile("data/P.txt");
    ofstream Qfile("data/Q.txt");
    ofstream dPfile("data/dP.txt");
    ofstream dQfile("data/dQ.txt");

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
        cerr << "Unable to open P file" << std::endl;
        return 1;
    }
    if (!dPfile.is_open()) {
        cerr << "Unable to open P file" << std::endl;
        return 1;
    }
    if (!dQfile.is_open()) {
        cerr << "Unable to open P file" << std::endl;
        return 1;
    }

    // Write the vector to the files
    for (int i=0; i<sp; ++i) {
        for (int j=0; j<N; ++j) {
            Pfile << P[i][j] << " ";
            Qfile << Q[i][j] << " ";
            dPfile << dP[i][j] << " ";
            dQfile << dQ[i][j] << " ";
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

