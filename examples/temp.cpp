#include <spheroidal/spheroidal_coordinate_functions.h>
#include <spheroidal/spheroidal_double_layer.h>
#include <spheroidal/grid_functions.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>  
#include <math.h> 
#include <yawg/utils.hpp>
#include <yawg/core.h>
#include <chrono>
using namespace std;
using namespace std::chrono;


// Test surface density: sigma = exp(-cos(theta)^2)*sin(theta)^2*sin(phi)
gsl::matrix test_density(int p)
{
    gsl::matrix theta, phi;
    gl_grid(p, theta, phi);
    gsl::matrix sigma(2*p*(p+1),1);
    for (int i=0; i<p+1; i++)
    {
        for (int j=0; j<2*p; j++)
        {
            sigma(j*(p+1)+i,0) = exp(-cos(theta(i,j))*cos(theta(i,j)))*(sin(theta(i,j))*sin(theta(i,j)))*sin(phi(i,j));
        }
    }
    return sigma;
}

int main()
{
    int p=16;
    double u0=2/sqrt(3);
    int peval=8, npeval=2*peval*(peval+1); 
    
    //Set up evaluation points
    // ------------------------------------------------------------------------
    gsl::matrix THETA_eval, PHI_eval;
    gl_grid(peval, THETA_eval, PHI_eval);
    gsl::vector u_coincident(npeval), theta_eval(npeval), v_eval(npeval), phi_eval(npeval);
    u_coincident=gsl::linspace(u0,u0,npeval);

    int j1=0;
    for (int j=0; j<PHI_eval.ncols(); j++)
    {
        theta_eval.subvector(j*(peval+1),peval+1)=THETA_eval.column(j);
        phi_eval.subvector(j*(peval+1),peval+1)=PHI_eval.column(j);
    }
    for (int i=0; i<npeval; i++)
    {
        v_eval(i)=cos(theta_eval(i));
    }

    gsl::matrix  S_coincident(npeval,3), S_far(npeval,3), S_near(npeval,3);
    S_coincident.column(0)=u_coincident; S_coincident.column(1)=v_eval; S_coincident.column(2)=phi_eval;
    S_far=S_coincident; S_far.column(0)=3*u_coincident;
    S_near=S_coincident; S_near.column(0)=1.1*u_coincident;
    // ------------------------------------------------------------------------

    gsl::matrix sigma = test_density(p);
    cout << "p=" << p << endl;

    auto start = high_resolution_clock::now();
    gsl::cmatrix DL_coincident = spheroidal_double_layer(sigma,u0, S_coincident, 1);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by coincident evaluation: " << duration.count() << " microseconds" << endl;

    start = high_resolution_clock::now();
    gsl::cmatrix DL_far = spheroidal_double_layer(sigma,u0, S_far, 1);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by far evaluation: " << duration.count() << " microseconds" << endl;



    std::string density_filename="../../tests/data/std/Density_near_16.csv";
    FILE* density_file = fopen(density_filename.c_str(), "r");
    gsl::matrix sigma_pcp; sigma_pcp.load_csv(density_file);
    fclose(density_file);

    start = high_resolution_clock::now();
    gsl::cmatrix DL_pcp = spheroidal_double_layer(sigma_pcp, u0, S_near, 1);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by near evaluation: " << duration.count() << " microseconds" << endl;
    
    return 0;

}


