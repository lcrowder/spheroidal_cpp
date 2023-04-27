#include <spheroidal/legendre_otc.h>
#include <spheroidal/grid_functions.h>
#include <spheroidal/spheroidal_harmonic_transforms.h>
#include <spheroidal/spheroidal_coordinate_functions.h>
#include <spheroidal/spheroidal_double_layer.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <yawg/utils.hpp>
#include <yawg/core.h>
#include <yawg/fft.h>
#include <yawg/legendre.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_sf_gamma.h>
using namespace std;

gsl::matrix solid_harmonic(int p, gsl::vector u_x, int region)
{
    int N=u_x.size();
    int sp=(p+1)*(p+1);
    gsl::matrix F(N,sp);
    gsl::matrix P;

    if (region==-1)
    {
        legendre_otc(p,u_x,P);
        F=P.T();
    }
    else if (region==0)
    {
        for (int i=0; i<N; i++)
            for (int j=0; j<sp; j++)
                F(i,j)=1.;
    }
    else if (region==1)
    {
        gsl::matrix Q;
        legendre_otc(p,u_x,P,Q);
        F=Q.T();
    }
    return F;
}

void DLspectrum(int p, double u0, gsl::vector &lambda_int, gsl::vector &lambda_surf, gsl::vector &lambda_ext)
{
    int sp=(p+1)*(p+1);
    gsl::matrix P, Q, dP, dQ;
    gsl::vector uvec(1); uvec(0)=u0;
    Dlegendre_otc(p,uvec,P,Q,dP,dQ);
    lambda_int.resize(sp);
    lambda_surf.resize(sp);
    lambda_ext.resize(sp);

    for (int i=0; i<sp; i++)
    {
        int n=floor(sqrt(i));
        int m=i-n*n-n;
        double anm =gsl_sf_fact(n-m)/(1.*gsl_sf_fact(n+m))*pow(-1,m)*(u0*u0-1.);

        lambda_int(i) =  anm * dQ(i,0);
        lambda_surf(i) = anm/2. * (P(i,0)*dQ(i,0)+dP(i,0)*Q(i,0));
        lambda_ext(i) =  anm * dP(i,0);
    }
}

gsl::cmatrix Ynm_matrix(int p, gsl::vector v, gsl::vector phi)
{
    int N=v.size();
    int sp=(p+1)*(p+1);
    gsl::cmatrix Y(N,sp);
    
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<sp; j++)
        {
            int n=floor(sqrt(j));
            int m=j-n*n-n;
            Y(i,j)=gsl::spherical_harmonic(n,m,acos(v(i)),phi(i));
        }
    }
    return Y;
}


gsl::cmatrix spheroidal_double_layer(gsl::matrix sigma, double u0, gsl::matrix X, int target_coords)
{
    gsl::cmatrix shc = spheroidal_analysis(sigma);
    int sp = shc.nrows();
    int num_surfs = shc.ncols(); //Number of surfaces
    int p = sqrt(sp)-1;
    int nt = X.nrows(); //Number of targets

    if (X.ncols()!=3*num_surfs){ cerr << "number of columns in X must be 3*num_surfs" << endl; }
    gsl::matrix S = X;

    if (target_coords==0)
    {
        S = cart_to_spheroidal(X,1/u0);
    }
    else if (target_coords > 1)
    {
        cerr << "invalid coordinates given. Select 0 for cartesian or 1 for spheroidal coordinates." << endl;
    }

    // Calculate spectra for interior, exterior, surface
    gsl::vector spectra_int, spectra_surf, spectra_ext;
    DLspectrum(p,u0,spectra_int,spectra_surf,spectra_ext);

    vector<gsl::vector> spectra_regions={spectra_int,spectra_surf,spectra_ext};

    gsl::cmatrix DL(nt,num_surfs);

    // Loop over all the different surfaces we want to evaluate
    for (int k=0; k<num_surfs; k++)
    {
        gsl::vector u_x = S.column(3*k);
        gsl::vector v_x = S.column(3*k+1);
        gsl::vector phi_x = S.column(3*k+2);

        // Count number of target which are in each region. 
        int int_count, surf_count, ext_count;
        int_count=0; surf_count=0; ext_count=0;
        for (int j=0; j<nt; j++)
        {
            if (u_x(j)<u0) int_count++;
            else if (u_x(j)==u0) surf_count++;
            else ext_count++;
        }
        
        //Sort targets into different regimes and get global indices of targets in each regime
        vector<int> regions = {int_count , surf_count , ext_count};
        // if (regions[0]>0) {gsl::matrix S_int(int_count,3); vector<int> j_int;}
        // if (regions[1]>0) {gsl::matrix S_surf(surf_count,3); vector<int> j_surf;}
        // if (regions[2]>0) {gsl::matrix S_ext(ext_count,3); vector<int> j_ext;}
        gsl::matrix S_int(int_count,3); vector<int> j_int;
        gsl::matrix S_surf(surf_count,3); vector<int> j_surf;
        gsl::matrix S_ext(ext_count,3); vector<int> j_ext;
        
        for (int j=0; j<nt; j++)
        {
            if (u_x(j)<u0)
            {
                S_int.row(int_count)=S.row(j);
                j_int.push_back(j);
            }
            else if (u_x(j)==u0)
            {
                S_surf.row(surf_count)=S.row(j);
                j_surf.push_back(j);
            }
            else
            {   
                S_ext.row(ext_count)=S.row(j);
                j_ext.push_back(j);
            }
        }
        
        vector<gsl::matrix> Sregions = {S_int, S_surf, S_ext};
        vector<gsl::cmatrix> DLregions;
        
        // Loop over all targets
        for (int r=0; r<3; r++)
        {
            // If there are more than 0 targets in this region, then we compute double layer potential
            int nt_r = regions[r];
            if (nt_r > 0)
            {
                // spheroidal coordinates of targets in this region
                gsl::matrix Sr=Sregions[r];
                gsl::vector u_x_r=Sr.column(0);
                gsl::vector v_x_r=Sr.column(1);
                gsl::vector phi_x_r=Sr.column(2);

                gsl::matrix Fr=solid_harmonic(p, u_x_r, r-1);
                gsl::cmatrix Yr=Ynm_matrix(p, v_x_r, phi_x_r);

                // Loop over targets in this region
                // for (int j=0; j<nt_r; j++)
                // {
                    
                // }
            }
        }

    }
    return DL;
}




// MATLAB CODE
//-----------------------------------------------------------------------------
// ...

    
//             FYr=Fr.*Yr;
//             DLcoefs_r=spectra_regions{r}.*shc(:,k);
//             DLregions{r}=FYr*DLcoefs_r;
//         end
//     end
    
//     % For debugging:
// %---------------------------------------------------------------
// %         sprintf('size Y: ')
// %         disp(size(Y))
// %         sprintf('size spectra: ')
// %         disp(size(spectra))
// %         sprintf('size shc: ')
// %         disp(size(shc))
// %         sprintf('size FY: ')
// %         disp(size(FY))
// %---------------------------------------------------------------

//     % Recombine all DLs from interior/exterior/surface and resort to
//     % original order of targets
//     DLk_unsorted=cat(1,DLregions{1},DLregions{2},DLregions{3});
//     DL(:,k)=DLk_unsorted(j_unsort,:);

// %     fprintf('completed evaluation of surface density %d\n',k)
// end


