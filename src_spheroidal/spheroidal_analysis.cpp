// #include <spheroidal/legendre_otc.h>
#include <spheroidal/grid_functions.h>
#include <spheroidal/spheroidal_analysis.h>
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
        L.row(n+fabs(m)) = wP; 
    }
    return L;
}


gsl::matrix spheroidal_analysis(gsl::matrix f)
{
    int np=f.nrows(); //number of points on a spheroid surface
    int ns=f.ncols(); //number of surfaces to evaluate
    int p = floor((sqrt(2.*np+1)-1.)/2.);
    int sp = (p+1)*(p+1);

    if (p-(sqrt(2.*np+1)-1.)/2. != 0 )
    { 
        cerr << "The input should be defined on a Gauss-Legendre--uniform grid. Check the size of input.\n";
    }
    gsl::matrix shc(sp,ns);

    f.print();

    for (int i=0; i<ns; i++)
    {
        // Pull out the ith column of f, corrresponding to the ith surface
        gsl::vector fi = f.column(i);
        fi.print();

        // Reshape the vector into a matrix with 2p rows and p+1 columns
        gsl::matrix f_matrix(p+1, 2*p);
        for (int j=0; j<2*p; j++)
        {
            // gsl::vector fi_col = fi.subvector(j*(p+1),p+1);
            // f_matrix.column(j) = fi_col;
            f_matrix.column(j) = fi.subvector(j*(p+1),p+1);
        }
        f_matrix.print();

        // Take the Fourier transform along the rows of the matrix
        gsl::cmatrix fi_fourier = gsl::fft(f_matrix, 2);
        fi_fourier.print();

        // Reshape the matrix of fourier coefficients into a vector (column major). 
        gsl::cvector fi_fourier_vec((2*p)*(p+1));
        for (int j=0; j<2*p; j++)
        {
            gsl::cvector fi_fourier_col = fi_fourier.column(j);
            fi_fourier_vec.subvector(j*(p+1),p+1) = fi_fourier_col;
        }
        fi_fourier_vec *= (M_PI/p);
        fi_fourier_vec.print();

        // Shift the vector so that m=0 frequency is in the center
        gsl::cvector fi_fourier_vec1 = gsl::circshift(fi_fourier_vec, p*(p+1));
        fi_fourier_vec1.print();
        
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
        fi_fourier_vec2.print();


        //Legendre basis transform
        gsl::cvector shc_i(fi_fourier_vec2.size());
        for (int m=-p; m<p+1; m++)
        {
            gsl::cvector fi_fourier_col = fi_fourier_vec2.subvector((p+m)*(p+1),p+1);
            gsl::cmatrix L = get_legendre_matrix(p,m);
            
            //// FIX THIS
            gsl::cvector shc_col = L*fi_fourier_col;
            //// FIX THIS

            // shc_col.print();
            
            shc_i.subvector((m+p)*(p+1),p+1) = shc_col;
        }
        shc_i.print();
        
        
        
        // SHRINK shc_i

        
        



    }

    return shc;

}

