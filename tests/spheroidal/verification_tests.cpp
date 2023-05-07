#define CATCH_CONFIG_MAIN
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
#include <catch2/catch_test_macros.hpp>
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


TEST_CASE("Convergence Testing", "[spheroidal_double_layer]")
{
    gsl::vector p_array = gsl::arange(2, 16, 1);
    double u0=2/sqrt(3);
    int peval=8, speval=(peval+1)*(peval+1), npeval=2*peval*(peval+1); 
    
    //Set up evaluation points
    // ------------------------------------------------------------------------
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
    // ------------------------------------------------------------------------

    // Load far evaluation data
    // ------------------------------------------------------------------------
    FILE* far_std_file = fopen("../../tests/data/std/DL_far_16_std.csv", "r");
    gsl::cmatrix DL_far_std;
    DL_far_std.load_csv(far_std_file);
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
    gsl::matrix DL_far_err(npeval-2*(peval+1), p_array.size());
    gsl::matrix DL_coincident_err(npeval-2*(peval+1), p_array.size());
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
        for (size_t i = 0; i < npeval - 2*(peval+1); ++i)
            DL_far_err(i,l)=log10((DL_far(i,0)-DL_far_std(i,0)).abs()); 

        for (size_t i = 0; i < npeval - 2*(peval+1); ++i)
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
    double min_tol_slope=-0.4, max_tol_slope=-0.2, max_tol_RSS=12.;
    double far_min_slope=0., far_max_slope=0., coincident_min_slope=0., coincident_max_slope=0., near_min_slope=0., near_max_slope=0.;
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
        
        if (i < npeval-2*(peval+1))
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

    // Check that the slopes are within the tolerance
    REQUIRE(far_min_slope > min_tol_slope) ;
    REQUIRE(far_max_slope < max_tol_slope) ;
    REQUIRE(coincident_min_slope > min_tol_slope) ;
    REQUIRE(coincident_max_slope < max_tol_slope) ;
    REQUIRE(near_min_slope > min_tol_slope) ;
    REQUIRE(near_max_slope < max_tol_slope) ;

    // Check that the RSS is within the tolerance
    REQUIRE(far_max_RSS < max_tol_RSS) ;
    REQUIRE(coincident_max_RSS < max_tol_RSS) ;
    REQUIRE(near_max_RSS < max_tol_RSS) ;

}




// OLD MAIN ROUTINE:

// int main()
// {
//     gsl::vector parr = gsl::arange(2, 17, 1);
//     double u0=2/sqrt(3);
//     int peval=8, speval=(peval+1)*(peval+1), npeval=2*peval*(peval+1); 
    
//     //Evaluation points
//     gsl::matrix THETA_eval, PHI_eval;
//     gl_grid(peval, THETA_eval, PHI_eval);
//     gsl::vector u_far(npeval), u_coincident(npeval), u_near(npeval), theta_eval(npeval), v_eval(npeval), phi_eval(npeval);
//     u_coincident=gsl::linspace(u0,u0,npeval);
//     u_far = 3*u_coincident;
//     u_near = 1.1*u_coincident;
//     for (int j=0; j<PHI_eval.ncols(); j++)
//     {
//         theta_eval.subvector(j*(peval+1),peval+1)=THETA_eval.column(j);
//         phi_eval.subvector(j*(peval+1),peval+1)=PHI_eval.column(j);
//     }
//     for (int i=0; i<npeval; i++)
//     {
//         v_eval(i)=cos(theta_eval(i));
//     }
//     gsl::matrix S_far(npeval,3), S_near(npeval,3), S_coincident(npeval,3);
//     S_far.column(0)=u_far; S_far.column(1)=v_eval; S_far.column(2)=phi_eval;
//     S_near.column(0)=u_near; S_near.column(1)=v_eval; S_near.column(2)=phi_eval;
//     S_coincident.column(0)=u_coincident; S_coincident.column(1)=v_eval; S_coincident.column(2)=phi_eval;

//     // Load far evaluation data
//     FILE* far_std_file = fopen("../../tests/data/std/DL_far_16_std.csv", "r");
//     gsl::cmatrix DL_far_std;
//     DL_far_std.load_csv(far_std_file);

//     //Load coincident evaluation data
//     FILE* coincident_std_file = fopen("../../tests/data/std/DL_coincident_16_std.csv", "r");
//     gsl::cmatrix DL_coincident_std;
//     DL_coincident_std.load_csv(coincident_std_file);


//     // point charges potential to compare with solution to BIE in data files
//     gsl::vector ptch(3); ptch(0)=1.; ptch(1)=2.; ptch(2)=-0.5;
//     gsl::matrix Xptch(3,3);
//     Xptch.row(0)=gsl::linspace(0,0,3);
//     Xptch.row(1)=gsl::linspace(.1,.1,3);
//     Xptch(2,0)=0.1; Xptch(2,1)=-0.1; Xptch(2,2)=-0.2;
//     gsl::matrix pcp = PtChargePotential(ptch, Xptch, spheroidal_to_cart(S_near,1/u0));
//     FILE* pcp_file = fopen("../../tests/data/PointChargePotential.csv", "w");
//     pcp.print_csv(pcp_file);
//     fclose(pcp_file);
    
//     printf("wrote pcp file\n");

//     //Completion terms for near singular test
//     FILE* completion_file = fopen("../../tests/data/std/NearCompletionTerms.csv", "r");
//     gsl::matrix C; 
//     C.load_csv(completion_file);
//     fclose(completion_file);

//     printf("read completion terms\n");

//     gsl::matrix DL_far_err(npeval-2*(peval+1), parr.size()-1);
//     gsl::matrix DL_coincident_err(npeval-2*(peval+1), parr.size());
//     gsl::matrix DL_near_err(npeval, parr.size());

//     // Increase the order p to see convergence
//     for (int l=0; l<parr.size(); l++)
//     {
//         int p = floor(parr(l));
//         int sp=(p+1)*(p+1);
//         int np=2*p*(p+1);

//         gsl::matrix sigma = test_density(p);

//         auto pstr = std::to_string(p);
        
//         // Far evaluation
//         gsl::cmatrix DL_far = spheroidal_double_layer(sigma,u0, S_far, 1);
//         std::string far_filename="../../tests/data/DL_far_"+pstr+".csv";
//         FILE* far_file = fopen(far_filename.c_str(), "w");
//         DL_far.print_csv(far_file);
//         fclose(far_file);

//         // Coincident evaluation
//         gsl::cmatrix DL_coincident = spheroidal_double_layer(sigma,u0, S_coincident, 1);
//         std::string coincident_filename="../../tests/data/DL_coincident_"+pstr+".csv";
//         FILE* coincident_file = fopen(coincident_filename.c_str(), "w");
//         DL_coincident.print_csv(coincident_file);
//         fclose(coincident_file);

//         for (int i=0; i<npeval; i++)
//         {
//             if (phi_eval(i)!=0. and phi_eval(i)!=M_PI)
//             {
//                 if (p<16)
//                 {
//                     DL_far_err(i,l)=log10((DL_far(i,0)-DL_far_std(i,0)).abs()); 
//                 }
//                 DL_coincident_err(i,l)=log10((DL_coincident(i,0)-DL_coincident_std(i,0)).abs()); 
//             }
//         }


//         // NEAR evaluation
//         std::string density_filename="../../tests/data/std/Density_near_"+pstr+".csv";
//         FILE* density_file = fopen(density_filename.c_str(), "r");
//         gsl::matrix sigma_pcp; sigma_pcp.load_csv(density_file);
//         gsl::cmatrix DL_pcp = spheroidal_double_layer(sigma_pcp, u0, S_near, 1);
//         fclose(density_file);

//         DL_pcp+=C.column(l); // D[sigma] + C is the solution to the BIE

//         std::string near_filename="../../tests/data/Solution_near_"+pstr+".csv";
//         FILE* near_file = fopen(near_filename.c_str(), "w");
//         DL_pcp.print_csv(near_file);
//         fclose(near_file);

//         for (int i=0; i<npeval; i++)
//         {
//             DL_near_err(i,l)=log10((DL_pcp(i,0)-pcp(i,0)).abs()) ;
//         }
//     }

//     // Linear fit to the log of the error
//     double min_tol_slope=-0.4, max_tol_slope=-0.2, max_tol_RSS=12.;
//     double far_min_slope=0., far_max_slope=0., coincident_min_slope=0., coincident_max_slope=0., near_min_slope=0., near_max_slope=0.;
//     double far_max_RSS=0., near_max_RSS=0., coincident_max_RSS=0.;

//     for (int i=0; i<npeval; i++)
//     {
//         gsl::matrix lsq_X(parr.size(), 2); //For coincident and near
//         lsq_X.column(0)=gsl::linspace(1,1,parr.size());
//         lsq_X.column(1)=gsl::linspace(2,16,parr.size());

//         gsl::vector near_y=DL_near_err.row(i);

//         // ------------------------------------------------------------------------------------------
//         // Do a least squares linear fit on lsq_X with near_y.
//         // Assume the coefficients are output in a gsl::vector [intercept, slope] lsq_near_beta.
//         // ------------------------------------------------------------------------------------------
        
//         // Get slope of current linear fit
//         double near_slope=lsq_near_beta(1);
//         near_min_slope = min(near_min_slope, near_slope);
//         near_max_slope = max(near_max_slope, near_slope);

//         // Get the y values of the linear fit
//         lsq_near_y = lsq_X * lsq_near_beta;

//         // Compute the residual sum of squares ror goodness of fit
//         double near_RSS=0.;
//         for (int j=0; j<near_y.size(); j++)
//         {
//             near_RSS+=(far_y(j)-lsq_near_y(j))*(near_y(j)-lsq_near_y(j));
//         }
//         near_max_RSS=max(near_max_RSS, near_RSS);

        
//         if (i < npeval-2*(peval+1))
//         {
//             gsl::vector far_y=DL_far_err.row(i);
//             gsl::vector coincident_y=DL_coincident_err.row(i);

//             // Least squares variables
//             gsl::matrix lsq_far_X(far_y.size(), 2);
//             lsq_far_X.column(0)=gsl::linspace(1,1,far_y.size());
//             lsq_far_X.column(1)=gsl::linspace(2,15,far_y.size());

//             // ------------------------------------------------------------------------------------------
//             // Do a least squares linear fit on lsq_far_X with far_y and lsq_coincident_X with coincident_y.
//             // Assume the coefficients are output in a gsl::vector [intercept, slope] lsq_far_beta and lsq_coincident_beta.
//             // ------------------------------------------------------------------------------------------
            
//             // Get slope of current linear fit
//             double far_slope=lsq_far_beta(1);
//             far_min_slope = min(far_min_slope, far_slope);
//             far_max_slope = max(far_max_slope, far_slope);
//             double coincident_slope=lsq_coincident_beta(1);
//             coincident_min_slope = min(coincident_min_slope, coincident_slope);
//             coincident_max_slope = max(coincident_max_slope, coincident_slope);
            
//             lsq_far_y = lsq_far_X * lsq_far_beta;
//             lsq_coincident_y = lsq_X * lsq_coincident_beta;

//             // Calculate the residual sum of squares for goodness of fit
//             double far_RSS=0., coincident_RSS=0.;
//             for (int j=0; j<far_y.size(); j++)
//             {
//                 far_RSS+=(far_y(j)-lsq_far_y(j))*(far_y(j)-lsq_far_y(j));
//             }
//             for (int j=0; j<coincident_y.size(); j++)
//             {
//                 coincident_RSS+=(coincident_y(j)-lsq_coincident_y(j))*(coincident_y(j)-lsq_coincident_y(j));
//             }

//             far_max_RSS=max(far_max_RSS, far_RSS);
//             coincident_max_RSS=max(coincident_max_RSS, coincident_RSS);
//         }
        


//     }

//     // Check that the slopes are within the tolerance
//     if (far_min_slope < min_tol_slope) ;
//     if (far_max_slope > max_tol_slope) ;
//     if (coincident_min_slope < min_tol_slope) ;
//     if (coincident_max_slope > max_tol_slope) ;
//     if (near_min_slope < min_tol_slope) ;
//     if (near_max_slope > max_tol_slope) ;

//     // Check that the RSS is within the tolerance
//     if (far_max_RSS > max_tol_RSS) ;
//     if (coincident_max_RSS > max_tol_RSS) ;
//     if (near_max_RSS > max_tol_RSS) ;

//     return 0;
// }

