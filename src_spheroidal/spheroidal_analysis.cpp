// #include <spheroidal/legendre_otc.h>
#include <spheroidal/grid_functions.h>
#include <spheroidal/spheroidal_analysis.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <yawg/utils.hpp>
#include <yawg/core.h>
#include <yawg/fft.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_sf_gamma.h>
using namespace std;


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
    gsl::matrix shc; 
    shc.resize(sp,ns);

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
            gsl::vector fi_col = fi.subvector(j*(p+1),p+1);
            f_matrix.column(j) = fi_col;
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
            gsl::cvector fi_fourier_col = fi_fourier_vec1.subvector(j,p+1);
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


        

        // gsl::cvector endpts = fi_fourier_vec.subvector(0,p+1);
        // endpts/=2.;
        
        // gsl::cvector fi_fourier_vec2(2*p*(p+1)+p+1);

        // fi_fourier_vec2.subvector(0,2*p*(p+1)) = fi_fourier_vec;
        // fi_fourier_vec2.subvector(2*p*(p+1), p+1) = endpts;

        // fi_fourier_vec2.print();


    }

    return shc;

}



//  MATLAB CODE:
// -----------------------------------------------------------------
// % SPHEROIDALANA(F) - Calculate the spheroidal harmonics transform of the function
// % set F. Each column of F is a function defined on the parameter domain
// % defined by parDomain.

// %d1 = number of points on a single particle surface
// %d2 = number of particles (columns of F)
// [d1,d2] = size(f);

// %d1 should be 2p(p+1)
// p = (sqrt(2*d1+1)-1)/2;

// if(p~=fix(p))
// error(['The input should be defined on a Gauss-Legendre--uniform'...
//        ' grid. Check the size of input.']);
// end

// %-- Allocating memory
// shc = zeros((2*p+1)*(p+1),d2);

// %-- Reshaping to match the parametrization grid structure
// f = reshape(f,p+1,2*p,d2);

// %-- DFT (identical to spherical harmonics)
// % Fourier transform along each row
// f = fft(f,[],2)*pi/p;
// f = circshift(reshape(f,2*p*(p+1),d2),p*(p+1));
// f(1:p+1,:) = f(1:p+1,:)/2;
// f = [f; f(1:p+1,:)];

// LT = []; LT_Freqs = [];
// LT{end+1} = getLegMat(p);
// LT_Freqs(end+1) = p;
// idx = length(LT);

// for m = -p:p
//     ind = (m+p)*(p+1)+1;
//     shc(ind:ind+p,:) = LT{idx}{p+1+m}*f(ind:ind+p,:);
// end
// % sprintf("before shrinkvec\n")
// % disp(shc)
// % disp(length(shc))

// shc = shrinkShVec(shc);

// % sprintf("after shrinkvec\n")
// % disp(shc)
// % disp(length(shc))

// end