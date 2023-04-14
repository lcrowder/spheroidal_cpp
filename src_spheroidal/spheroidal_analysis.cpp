#include <spheroidal/legendre_otc.h>
#include <spheroidal/grid_functions.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <gsl_wrapper/utils.hpp>
#include <gsl_wrapper/core.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_sf_gamma.h>
using namespace std;


gsl::vector spheroidal_analysis(gsl::matrix f)
{
    int np=f.nrows(); //number of points on a spheroid surface
    int ns=f.ncols(); //number of surfaces to evaluate
    int p = floor((sqrt(2.*d1+1)-1.)/2.);
    int sp = (p+1)*(p+1);

    if (p-(sqrt(2.*d1+1)-1.)/2. != 0 )
    { 
        cerr << "The input should be defined on a Gauss-Legendre--uniform grid. Check the size of input.";
    }
    
    gsl::matrix shc; 
    shc.resize(sp,ns);

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