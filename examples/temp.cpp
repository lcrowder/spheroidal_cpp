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
#include <yawg/lls.h>
#include <string>
#include <unistd.h>
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
    gsl::vector p_array = gsl::arange(2, 16, 1);
    double u0=2/sqrt(3);
    int peval=8, npeval=2*peval*(peval+1), npeval_short=npeval-2*(peval+1); 
    
    //Set up evaluation points
    // ------------------------------------------------------------------------
    gsl::matrix THETA_eval, PHI_eval;
    gl_grid(peval, THETA_eval, PHI_eval);
    gsl::vector u_far(npeval_short), u_coincident(npeval_short), theta_eval(npeval_short), v_eval(npeval_short), phi_eval(npeval_short);
    gsl::vector u_near(npeval), theta_near(npeval), v_near(npeval), phi_near(npeval);
    u_coincident=gsl::linspace(u0,u0,npeval_short);
    u_far = 3*u_coincident;
    u_near = 1.1*gsl::linspace(u0,u0,npeval);

    int j1=0;
    for (int j=0; j<PHI_eval.ncols(); j++)
    {
        theta_near.subvector(j*(peval+1),peval+1)=THETA_eval.column(j);
        phi_near.subvector(j*(peval+1),peval+1)=PHI_eval.column(j);
        if (j!=0 and j!=peval)
        {
            theta_eval.subvector(j1*(peval+1),peval+1)=THETA_eval.column(j);
            phi_eval.subvector(j1*(peval+1),peval+1)=PHI_eval.column(j);
            j1+=1;
        }
    }
    for (int i=0; i<npeval; i++)
    {
        v_near(i)=cos(theta_near(i));
        if (i < npeval_short)
        {
            v_eval(i)=cos(theta_eval(i));
        }
    }

    gsl::matrix S_far(npeval_short,3), S_near(npeval,3), S_coincident(npeval_short,3);
    S_far.column(0)=u_far; S_far.column(1)=v_eval; S_far.column(2)=phi_eval;
    S_near.column(0)=u_near; S_near.column(1)=v_near; S_near.column(2)=phi_near;
    S_coincident.column(0)=u_coincident; S_coincident.column(1)=v_eval; S_coincident.column(2)=phi_eval;
    // ------------------------------------------------------------------------


    // Load far evaluation data
    // ------------------------------------------------------------------------
    FILE* far_std_file = fopen("../../tests/data/std/DL_far_16_std.csv", "r");
    gsl::matrix DL_far_std;
    DL_far_std.load_csv(far_std_file);
    DL_far_std.print();
    // ------------------------------------------------------------------------

    //Load coincident evaluation data
    // ------------------------------------------------------------------------
    FILE* coincident_std_file = fopen("../../tests/data/std/DL_coincident_16_std.csv", "r");
    gsl::cmatrix DL_coincident_std;
    DL_coincident_std.load_csv(coincident_std_file);
    // ------------------------------------------------------------------------


    // Compute point charges potential to compare with solution to BIE in data files
    // ------------------------------------------------------------------------
    gsl::vector ptch(3); ptch(0)=1.; ptch(1)=2.; ptch(2)=-0.5;
    gsl::matrix Xptch(3,3);
    Xptch.row(0)=gsl::linspace(0,0,3);
    Xptch.row(1)=gsl::linspace(.1,.1,3);
    Xptch(2,0)=0.1; Xptch(2,1)=-0.1; Xptch(2,2)=-0.2;
    gsl::matrix pcp = PtChargePotential(ptch, Xptch, spheroidal_to_cart(S_near,1/u0));
    FILE* pcp_file = fopen("../../tests/data/PointChargePotential.csv", "w");
    pcp.print_csv(pcp_file);
    fclose(pcp_file);
    // ------------------------------------------------------------------------

    //Completion terms for near singular test
    // ------------------------------------------------------------------------
    FILE* completion_file = fopen("../../tests/data/std/NearCompletionTerms.csv", "r");
    gsl::matrix C; 
    C.load_csv(completion_file);
    fclose(completion_file);
    // ------------------------------------------------------------------------


    // Calculate the errors for each regime and order p
    // ------------------------------------------------------------------------
    gsl::matrix DL_far_err(npeval_short, p_array.size());
    gsl::matrix DL_coincident_err(npeval_short, p_array.size());
    gsl::matrix DL_near_err(npeval, p_array.size());

    // Increase the order p to see convergence
    for (int l=0; l<p_array.size(); l++)
    {
        int p = floor(p_array(l));
        int sp=(p+1)*(p+1);
        int np=2*p*(p+1);

        gsl::matrix sigma = test_density(p);

        auto pstr = std::to_string(p);
        
        // Write Far evaluation data. Not needed for this but Leo likes to have it for Matlab tests
        gsl::cmatrix DL_far = spheroidal_double_layer(sigma,u0, S_far, 1);
        std::string far_filename="../../tests/data/DL_far_"+pstr+".csv";
        FILE* far_file = fopen(far_filename.c_str(), "w");
        DL_far.print_csv(far_file);
        fclose(far_file);

        // Write coincident evaluation data. Not needed for this but Leo likes to have it for Matlab tests.
        gsl::cmatrix DL_coincident = spheroidal_double_layer(sigma,u0, S_coincident, 1);
        std::string coincident_filename="../../tests/data/DL_coincident_"+pstr+".csv";
        FILE* coincident_file = fopen(coincident_filename.c_str(), "w");
        DL_coincident.print_csv(coincident_file);
        fclose(coincident_file);

        // Compute log10 error for far evaluation
        for (size_t i = 0; i < npeval_short; ++i)
            DL_far_err(i,l)=log10((DL_far(i,0)-DL_far_std(i,0)).abs()); 

        for (size_t i = 0; i < npeval_short; ++i)
            DL_coincident_err(i,l)=log10((DL_coincident(i,0)-DL_coincident_std(i,0)).abs());

        // Near evaluation
        // read in densities (sigma) obtained from solving the BIE on the spheroid with the singluar quadrature matrix operator
        std::string density_filename="../../tests/data/std/Density_near_"+pstr+".csv";
        FILE* density_file = fopen(density_filename.c_str(), "r");
        gsl::matrix sigma_pcp; sigma_pcp.load_csv(density_file);
        gsl::cmatrix DL_pcp = spheroidal_double_layer(sigma_pcp, u0, S_near, 1);
        fclose(density_file);

        DL_pcp+=C.column(l); // D[sigma] + C is the solution to the BIE

        //Not needed for this but Leo likes to have it for Matlab tests.
        std::string near_filename="../../tests/data/Solution_near_"+pstr+".csv";
        FILE* near_file = fopen(near_filename.c_str(), "w");
        DL_pcp.print_csv(near_file);
        fclose(near_file);

        // Compute log10 error for near singular evaluation
        for (int i=0; i<npeval; i++)
        {
            DL_near_err(i,l)=log10((DL_pcp(i,0)-pcp(i,0)).abs()) ;
        }
    }

    // Linear fit to the log of the error.
    // ------------------------------------------------------------------------
    // From Matlab tests, slope should be between -0.4 and -0.2 for all regimes
    // From Matlab tests, RSS should be less than 12 for all regimes 
    double min_tol_slope=-0.45, max_tol_slope=-0.2, max_tol_RSS=12.;
    double far_min_slope=100., far_max_slope=-100., coincident_min_slope=100., coincident_max_slope=-100., near_min_slope=100., near_max_slope=-100.;
    double far_max_RSS=0., near_max_RSS=0., coincident_max_RSS=0.;
    
    // Loop over target points. We will do a linear least squares fit for each one.
    for (int i=0; i<npeval; i++)
    {
        gsl::vector near_y = DL_near_err.row(i);

        double lsq_near_beta0, lsq_near_beta1, near_RSS;
        gsl::fit_linear( p_array, near_y, lsq_near_beta0, lsq_near_beta1, near_RSS );

        // Get slope of current linear fit
        double near_slope=lsq_near_beta1;
        near_min_slope = min(near_min_slope, near_slope);
        near_max_slope = max(near_max_slope, near_slope);

        near_max_RSS=max(near_max_RSS, near_RSS);
        
        if (i < npeval_short)
        {
            gsl::vector far_y=DL_far_err.row(i);
            gsl::vector coincident_y=DL_coincident_err.row(i);
            
            double lsq_far_beta0, lsq_far_beta1, far_RSS;
            gsl::fit_linear( p_array, far_y, lsq_far_beta0, lsq_far_beta1, far_RSS );

            double lsq_coincident_beta0, lsq_coincident_beta1, coincident_RSS;
            gsl::fit_linear( p_array, coincident_y, lsq_coincident_beta0, lsq_coincident_beta1, coincident_RSS );

            // Get slope of current linear fit
            double far_slope=lsq_far_beta1;
            far_min_slope = min(far_min_slope, far_slope);
            far_max_slope = max(far_max_slope, far_slope);

            double coincident_slope=lsq_coincident_beta1;
            coincident_min_slope = min(coincident_min_slope, coincident_slope);
            coincident_max_slope = max(coincident_max_slope, coincident_slope);

            far_max_RSS=max(far_max_RSS, far_RSS);
            coincident_max_RSS=max(coincident_max_RSS, coincident_RSS);
        }

    }

    printf("far_min_slope=%f\n",far_min_slope);
    printf("far_max_slope=%f\n",far_max_slope);
    printf("coincident_min_slope=%f\n",coincident_min_slope);
    printf("coincident_max_slope=%f\n",coincident_max_slope);
    printf("near_min_slope=%f\n",near_min_slope);
    printf("near_max_slope=%f\n",near_max_slope);
    printf("far_max_RSS=%f\n",far_max_RSS);
    printf("coincident_max_RSS=%f\n",coincident_max_RSS);
    printf("near_max_RSS=%f\n",near_max_RSS);

    return 0;

}

