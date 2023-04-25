// #include <spheroidal/legendre_otc.h>
#include <spheroidal/grid_functions.h>
#include <spheroidal/spheroidal_harmonic_transforms.h>
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


gsl::matrix get_legendre_matrix(int p, int m)
{
    //Get the Gauss-Legendre quadrature points and weights
    gsl::vector v, w;
    gsl::leggauss(p+1, v, w);
    gsl::vector wP(p+1);
    gsl::matrix L(p+1,p+1);

    for (int n=fabs(m); n<p+1; n++)
    {
        for (int j=0; j<p+1; j++)
        {
            wP(j) = w(j) * gsl::spherical_harmonic(n, m, v(j));
        }
        L.row(n) = wP; 
    }
    return L;
}

gsl::matrix get_legendre_matrix_inv(int p, int m)
{
    //Get the Gauss-Legendre quadrature points and weights
    gsl::vector v, w;
    gsl::leggauss(p+1, v, w);
    gsl::vector P(p+1);
    gsl::matrix L_inv(p+1,p+1);

    for (int n=fabs(m); n<p+1; n++)
    {
        for (int j=0; j<p+1; j++)
        {
            P(j) = gsl::spherical_harmonic(n, m, v(j));
        }
        L_inv.column(n) = P; 
    }
    return L_inv;
}

// Get index in spheroidal harmonic ordering: n=0,1,2,...,p; m=-n,-n+1,...,n
int geti(int n, int m)
{
    return m + n * (n + 1);
}


gsl::cmatrix spheroidal_analysis(gsl::matrix f)
{
    int np=f.nrows(); //number of points on a spheroid surface
    int ns=f.ncols(); //number of surfaces to evaluate
    int p = floor((sqrt(2.*np+1)-1.)/2.);
    int sp = (p+1)*(p+1);

    if (p-(sqrt(2.*np+1)-1.)/2. != 0 )
    { 
        cerr << "The input should be defined on a Gauss-Legendre--uniform grid. Check the size of input.\n";
    }
    gsl::cmatrix shc(sp,ns);

    // f.print();

    for (int i=0; i<ns; i++)
    {
        // Pull out the ith column of f, corrresponding to the ith surface
        gsl::vector fi = f.column(i);
        // fi.print();

        // Reshape the vector into a matrix with 2p rows and p+1 columns
        gsl::matrix f_matrix(p+1, 2*p);
        for (int j=0; j<2*p; j++)
        {
            // gsl::vector fi_col = fi.subvector(j*(p+1),p+1);
            // f_matrix.column(j) = fi_col;
            f_matrix.column(j) = fi.subvector(j*(p+1),p+1);
        }
        // f_matrix.print();

        // Take the Fourier transform along the rows of the matrix
        gsl::cmatrix fi_fourier = gsl::fft(f_matrix, 2);
        // fi_fourier.print();

        // Reshape the matrix of fourier coefficients into a vector (column major). 
        gsl::cvector fi_fourier_vec((2*p)*(p+1));
        for (int j=0; j<2*p; j++)
        {
            gsl::cvector fi_fourier_col = fi_fourier.column(j);
            fi_fourier_vec.subvector(j*(p+1),p+1) = fi_fourier_col;
        }
        fi_fourier_vec *= (M_PI/p);
        // fi_fourier_vec.print();

        // Shift the vector so that m=0 frequency is in the center
        gsl::cvector fi_fourier_vec1 = gsl::circshift(fi_fourier_vec, p*(p+1));
        // fi_fourier_vec1.print();
        
        //Set up for trapezoidal rule.
        gsl::cvector fi_fourier_vec2((2*p+1)*(p+1));
        for (int j=0; j<2*p; j++)
        {
            gsl::cvector fi_fourier_col = fi_fourier_vec1.subvector(j*(p+1),p+1);
            if (j==0)
            {
                fi_fourier_col /= 2.;
                fi_fourier_vec2.subvector(0,p+1) = fi_fourier_col;
                fi_fourier_vec2.subvector(2*p*(p+1),p+1) = fi_fourier_col;
            }
            else
            {
                fi_fourier_vec2.subvector(j*(p+1),p+1) = fi_fourier_col;
            }
        }
        // fi_fourier_vec2.print();

        // printf("Starting legendre basis transformation now.\n");

        //Legendre basis transform
        gsl::cmatrix shc_i(p+1,2*p+1);
        for (int m=-p; m<p+1; m++)
        {
            gsl::cvector fi_fourier_col = fi_fourier_vec2.subvector((p+m)*(p+1),p+1);
            // fi_fourier_col.print();
            gsl::cmatrix L = get_legendre_matrix(p,m);
            // L.print();

            gsl::cvector shc_col = L*fi_fourier_col;
            // shc_col.print();
            
            shc_i.column((m+p)) = shc_col;
        }
        // shc_i.print();
        
        
        // SHRINK shc_i
        gsl::cvector shrink_shc_i(sp);
        int count=0;
        for (int n=0; n<p+1; n++)
        {
            for (int m=-n; m<n+1; m++)
            {
                int inm=n;
                int jnm=m+p;
                if (fabs(m)<=n)
                {
                    shrink_shc_i(count) = shc_i(inm,jnm);
                    count++;
                }
                
            }
        }
        // shrink_shc_i.print();

        shc.column(i)=shrink_shc_i;
    }
    return shc;
}



gsl::cmatrix spheroidal_snythesis(gsl::cmatrix shc)
{
    int sp=shc.nrows(); //number of points on a spheroid surface
    int ns=shc.ncols(); //number of surfaces to evaluate
    int p = floor(sqrt(sp)-1.);
    int np = 2*p*(p+1);

    if (p-(sqrt(sp)-1.) != 0 )
    { 
        cerr << "The input should be defined on a Gauss-Legendre--uniform grid. Check the size of input.\n";
    }

    gsl::cmatrix f(np,ns);

    // Loop over surfaces
    for (int i=0; i<ns; i++)
    {
        gsl::cvector shc_i = shc.column(i);
        shc_i.print();

        //EXPAND shc
        gsl::cmatrix shc_i_matrix(p+1,2*p+1);
        for (int n=0; n<p+1; n++)
        {
            for (int m=-n; m<n+1; m++)
            {
                int idx=geti(n,m);
                int inm=n;
                int jnm=m+p;

                if (fabs(m)<=n)
                {
                    shc_i_matrix(inm,jnm) = shc_i(idx);
                }
                else
                {
                    shc_i_matrix(inm,jnm) = 0.;
                }
            }
        }
        shc_i_matrix.print();


        //Inverse Legendre basis transform. Since we chop the last frequency, this loop is shorter than the one in synthesis.
        // Also we need to multiply the m=-p frequency by 2.
        gsl::cvector fi_fourier_vec(np);
        for (int m=-p; m<p; m++)
        {
            gsl::cvector shc_i_col = shc_i_matrix.column((p+m));
            shc_i_col.print();
            
            gsl::cmatrix L_inv = get_legendre_matrix_inv(p,m);
            L_inv.print();

            gsl::cvector fi_fourier_m = L_inv*shc_i_col;
            fi_fourier_m.print();
            
            if (m==-p) {fi_fourier_m *= 2.;}
            fi_fourier_vec.subvector((m+p)*(p+1),p+1) = fi_fourier_m;
        }
        fi_fourier_vec.print();
        

        //Undo the steps we took rearranging f for the fft:
        // Reshape fourier coefficients into a matrix
        gsl::cvector fi_fourier_vec1 = gsl::circshift(fi_fourier_vec, -p*(p+1));
        fi_fourier_vec1.print();

        gsl::cmatrix fi_fourier(p+1,2*p);
        for (int j=0; j<2*p; j++)
        {
            fi_fourier.column(j) = fi_fourier_vec1.subvector(j*(p+1),p+1);
        }
        fi_fourier.print();
        
        //Inverse Fourier transform along the rows of the matrix
        gsl::cmatrix f_matrix = gsl::ifft(fi_fourier, 2);
        f_matrix *= 2*p;
        f_matrix.print();

        // Reshape the matrix back into a vector
        gsl::cvector fi(np);
        for (int j=0; j<2*p; j++)
        {
            fi.subvector(j*(p+1),p+1) = f_matrix.column(j) ;
        }
        fi.print();

        //Store in big matrix f
        f.column(i) = fi;
    }
    return f;
}

