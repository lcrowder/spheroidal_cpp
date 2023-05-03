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
#include <string>
using namespace std;


gsl::vector PtChargePotential(gsl::vector ptch, gsl::matrix Xptch, gsl::matrix Y)
{
    int M=ptch.size();
    int np=Y.nrows();
    gsl::vector pcp(np);
    for (int i=0; i<np; i++)
    {
        gsl::vector Yi=Y.row(i);
        pcp(i)=0;
        for (int j=0; j<M; j++)
        {
            gsl::vector Xj=Xptch.row(j);
            double r=0;
            for (int k=0; k<3; k++)
            {
                r+=pow(Xj(k)-Yi(k),2);
            }
            r=sqrt(r);
            pcp(i)+=ptch(j)/(4*M_PI*r);
        }
    }
    return pcp;
}

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
    gsl::vector parr = gsl::arange(2, 17, 1);
    double u0=2/sqrt(3);
    int peval=8, speval=(peval+1)*(peval+1), npeval=2*peval*(peval+1); 
    
    //Evaluation points
    gsl::matrix THETA_eval, PHI_eval;
    gl_grid(peval, THETA_eval, PHI_eval);
    gsl::vector u_far(npeval), u_coincident(npeval), u_near(npeval), theta_eval(npeval), v_eval(npeval), phi_eval(npeval);
    u_coincident=gsl::linspace(u0,u0,npeval);
    u_far = 3*u_coincident;
    u_near = 1.1*u_coincident;
    for (int j=0; j<PHI_eval.ncols(); j++)
    {
        theta_eval.subvector(j*(peval+1),peval+1)=THETA_eval.column(j);
        phi_eval.subvector(j*(peval+1),peval+1)=PHI_eval.column(j);
    }
    for (int i=0; i<npeval; i++)
    {
        v_eval(i)=cos(theta_eval(i));
    }
    gsl::matrix S_far(npeval,3), S_near(npeval,3), S_coincident(npeval,3);
    S_far.column(0)=u_far; S_far.column(1)=v_eval; S_far.column(2)=phi_eval;
    S_near.column(0)=u_near; S_near.column(1)=v_eval; S_near.column(2)=phi_eval;
    S_coincident.column(0)=u_coincident; S_coincident.column(1)=v_eval; S_coincident.column(2)=phi_eval;

    // point charges potential to compare with solution to BIE in data files
    gsl::vector ptch(3); ptch(0)=1.; ptch(1)=2.; ptch(2)=-0.5;
    gsl::matrix Xptch(3,3);
    Xptch.row(0)=gsl::linspace(0,0,3);
    Xptch.row(1)=gsl::linspace(.1,.1,3);
    Xptch(2,0)=0.1; Xptch(2,1)=-0.1; Xptch(2,2)=-0.2;
    gsl::matrix pcp = PtChargePotential(ptch, Xptch, spheroidal_to_cart(S_near,1/u0));
    FILE* pcp_file = fopen("../../tests/data/PointChargePotential.csv", "w");
    pcp.print_csv(pcp_file);
    fclose(pcp_file);
    
    printf("wrote pcp file\n");

    //Completion terms for near singular test
    FILE* completion_file = fopen("../../tests/data/std/NearCompletionTerms.csv", "r");
    gsl::matrix C; 
    C.load_csv(completion_file);
    fclose(completion_file);

    printf("read completion terms\n");

    // Increase the order p to see convergence
    for (int l=0; l<parr.size(); l++)
    {
        int p = floor(parr(l));
        int sp=(p+1)*(p+1);
        int np=2*p*(p+1);

        gsl::matrix sigma = test_density(p);

        auto pstr = std::to_string(p);
        
        // Far evaluation
        gsl::cmatrix DL_far = spheroidal_double_layer(sigma,u0, S_far, 1);
        std::string far_filename="../tests/data/DL_far_"+pstr+".csv";
        FILE* far_file = fopen(far_filename.c_str(), "w");
        DL_far.print_csv(far_file);
        fclose(far_file);

        // Coincident evaluation
        gsl::cmatrix DL_coincident = spheroidal_double_layer(sigma,u0, S_coincident, 1);
        std::string coincident_filename="../tests/data/DL_coincident_"+pstr+".csv";
        FILE* coincident_file = fopen(coincident_filename.c_str(), "w");
        DL_coincident.print_csv(coincident_file);
        fclose(coincident_file);

        // NEAR evaluation
        std::string density_filename="../tests/data/std/Density_near_"+pstr+".csv";
        FILE* density_file = fopen(density_filename.c_str(), "r");
        gsl::matrix sigma_pcp; sigma_pcp.load_csv(density_file);
        gsl::cmatrix DL_pcp = spheroidal_double_layer(sigma_pcp, u0, S_near, 1);
        fclose(density_file);

        DL_pcp+=C.column(l); // D[sigma] + C is the solution to the BIE

        std::string near_filename="../tests/data/Solution_near_"+pstr+".csv";
        FILE* near_file = fopen(near_filename.c_str(), "w");
        DL_pcp.print_csv(near_file);
        fclose(near_file);

    }

    return 0;
}

