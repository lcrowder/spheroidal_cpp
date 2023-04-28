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

    //

    gsl::matrix P,Q;
    legendre_otc(p,u_x,P,Q);

    for (int i=0; i<N; i++)
            for (int j=0; j<sp; j++)
                if (region==-1) F(i,j)=P(j,i);
                else if (region==1) F(i,j)=Q(j,i);
                else if (region==0) F(i,j)=1.;
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


gsl::cmatrix spheroidal_double_layer(gsl::cmatrix sigma, double u0, gsl::matrix X, int target_coords)
{
    printf("sigma after input:\n");
    sigma.print();
    
    gsl::cmatrix shc = spheroidal_analysis(sigma);
    printf("Computed spheroidal harmonic coefficients\n");
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
        printf("Evaluating surface %d\n",k);

        gsl::matrix Sk = S.submatrix(0,3*k,nt,3); // Get coordinates of targets coresponding to surface k
        gsl::vector u_x = Sk.column(0);
        gsl::vector v_x = Sk.column(1);
        gsl::vector phi_x = Sk.column(2);

        // Count number of target which are in each region. 
        //Sort targets into different regimes and get global indices of targets in each regime
        
        int int_count=0, surf_count=0, ext_count=0;
        vector<int> j_int, j_surf, j_ext;
        gsl::matrix S_int_padded(nt,3), S_surf_padded(nt,3), S_ext_padded(nt,3);    

        for (int j=0; j<nt; j++)
        {
            if (u_x(j)<u0)
            {
                S_int_padded.row(int_count) = Sk.row(j);
                j_int.push_back(j);
                int_count++;
            }
            else if (u_x(j)==u0)
            {
                S_surf_padded.row(surf_count)  = Sk.row(j);
                j_surf.push_back(j);
                surf_count++;
            }
            else
            {   
                S_ext_padded.row(ext_count)  = Sk.row(j);
                j_ext.push_back(j);
                ext_count++;
            }
        }
        
        printf("Number of targets in int/surf/ext regions: %d, %d, %d\n",int_count,surf_count,ext_count);

        gsl::matrix S_int, S_surf, S_ext;

        if (int_count > 0)
        {
            S_int.resize(int_count,3);
            S_int = S_int_padded.submatrix(0,0, int_count,3);
        }
        if (surf_count > 0)
        {
            S_surf.resize(surf_count,3);
            S_surf = S_surf_padded.submatrix(0,0, surf_count,3);
        }
        if (ext_count > 0)
        {
            S_ext.resize(ext_count,3);
            S_ext = S_ext_padded.submatrix(0,0, ext_count,3);
        }

        printf("Chopped up S matrices.\n");

        vector<gsl::matrix> Sregions = {S_int, S_surf, S_ext};
        vector<int> region_counts = {int_count , surf_count , ext_count};
        vector<vector<int> > jregions= {j_int, j_surf, j_ext};
        
        // Loop over each region
        for (int r=0; r<3; r++)
        {
            
            // If there are more than 0 targets in this region, then we compute double layer potential
            int nt_r = region_counts[r];

            if (nt_r > 0)
            {
                printf("Evaluating region %d\n",r-1);

                // spheroidal coordinates of targets in this region
                gsl::matrix Sr=Sregions[r];
                gsl::vector u_x_r=Sr.column(0);
                gsl::vector v_x_r=Sr.column(1);
                gsl::vector phi_x_r=Sr.column(2);

                gsl::matrix Fr=solid_harmonic(p, u_x_r, r-1);

                printf("Solid harmonics:\n");
                Fr.print();

                gsl::cmatrix Yr=Ynm_matrix(p, v_x_r, phi_x_r);

                printf("Spherical harmonics:\n");
                Yr.print();

                gsl::cvector DLcoefs_r = spectra_regions[r];
                gsl::cmatrix FYr (nt_r,sp);

                // Construct a matrix of solid harmonics and a vector of coefficients
                for (int i=0; i<sp; i++)
                {
                    DLcoefs_r(i)*=shc(i,k);
                    for (int j=0; j<nt_r; j++)
                    {
                        FYr(j,i)=Fr(j,i)*Yr(j,i);
                    }  
                }

                printf("Eigenfunctions:\n");
                FYr.print();

                printf("DL coefficients:\n");
                DLcoefs_r.print();

                // Multiply the two matrices to get the double layer potential
                gsl::cvector DLregions_r=FYr*DLcoefs_r;



                // Recombine all DLs from interior/exterior/surface
                for (int j=0; j<nt_r; j++)
                {
                    int j_global=jregions[r][j]; 
                    DL(j_global,k)=DLregions_r(j);
                }
            }
        }
    }
    return DL;
}


//Evaluate double layer potential at sigma, which is on surface.
gsl::cmatrix spheroidal_double_layer(gsl::cmatrix sigma, double u0)
{
    gsl::cmatrix shc = spheroidal_analysis(sigma);
    printf("Computed spheroidal harmonic coefficients\n");

    int sp = shc.nrows();
    int num_surfs = shc.ncols(); //Number of surfaces
    int p = sqrt(sp)-1;

    int nt=sigma.nrows(); //Number of targets

    // Calculate spectra for interior, exterior, surface
    gsl::vector spectra_int, spectra_surf, spectra_ext;
    DLspectrum(p,u0,spectra_int,spectra_surf,spectra_ext);

    gsl::cmatrix DL(nt,num_surfs);

    gsl::matrix V, PHI;
    gl_grid(p, V,PHI);
    gsl::vector u_x(nt), v_x(nt), phi_x(nt);
    u_x=gsl::linspace(u0,u0,nt);
    for (int j=0; j<V.ncols(); j++)
    {
        v_x.subvector(j*(p+1),p+1)=V.column(j);
        phi_x.subvector(j*(p+1),p+1)=PHI.column(j);
    }

    
    gsl::matrix F=solid_harmonic(p, u_x, 0);

    printf("Solid harmonics:\n");
    F.print();

    gsl::cmatrix Y=Ynm_matrix(p, v_x, phi_x);

    printf("Spherical harmonics:\n");
    Y.print();

    gsl::cvector DLcoefs = spectra_surf;
    gsl::cmatrix FY(nt,sp);

    // Loop over all the different surfaces we want to evaluate
    for (int k=0; k<num_surfs; k++)
    {
        printf("Evaluating surface %d\n",k);

        // Construct a matrix of solid harmonics and a vector of coefficients
        for (int i=0; i<sp; i++)
        {
            DLcoefs(i)*=shc(i,k);
            for (int j=0; j<nt; j++)
            {
                FY(j,i)=F(j,i)*Y(j,i);
            }  
        }

        printf("Eigenfunctions:\n");
        FY.print();

        printf("DL coefficients:\n");
        DLcoefs.print();

        // Multiply the two matrices to get the double layer potential
        DL.column(k) = FY*DLcoefs;
        
    }
    return DL;
}


//Allow vector input for sigma
gsl::cvector spheroidal_double_layer(gsl::cvector sigma, double u0, gsl::matrix X, int target_coords)
{
    gsl::cmatrix sigma_matrix(sigma.size(),1);
    sigma_matrix.column(0)=sigma;
    
    printf("Sigma vector:\n");
    sigma.print();
    printf("Sigma matrix:\n");
    sigma_matrix.print();

    gsl::cmatrix DL = spheroidal_double_layer(sigma_matrix,u0,X,target_coords);
    gsl::cvector DLvec=DL.column(0);

    printf("DL matrix:\n");
    DL.print();
    printf("DL vector:\n");
    DLvec.print();

    return DLvec;
}

gsl::cvector spheroidal_double_layer(gsl::cvector sigma, double u0)
{
    gsl::cmatrix sigma_matrix(sigma.size(),1);
    sigma_matrix.column(0)=sigma;
    gsl::cmatrix DL = spheroidal_double_layer(sigma_matrix,u0);
    gsl::cvector DLvec = DL.column(0);
    return DLvec;
}


