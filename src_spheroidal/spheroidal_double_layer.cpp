#include <spheroidal/legendre_otc.h>
#include <spheroidal/grid_functions.h>
#include <spheroidal/spheroidal_harmonic_transforms.h>
#include <spheroidal/spheroidal_coordinate_functions.h>
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

gsl::matrix solid_harmonic(int p, gsl::vector u_x, int regime=0)
{
    int N=u_x.size();
    int sp=(p+1)*(p+1);
    gsl::matrix F(N,sp);
    gsl::matrix P;

    if (regime==-1)
    {
        legendre_otc(p,u_x,P);
        F=P.T();
    }
    else if (regime==0)
    {
        for (int i=0; i<N; i++)
            for (int j=0; j<sp; j++)
                F(i,j)=1.;
    }
    else if (regime==1)
    {
        gsl::matrix Q;
        legendre_otc(p,u_x,P,Q);
        F=Q.T();
    }
    return F;
}


// function [lambda_int,lambda_surf,lambda_ext]=DLspectrum(p,u0)
//     sp=(p+1)^2;
//     ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
//     anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^mm .*(u0.^2-1);
//     L=legendre_otc(p,u0,1,1,1);
//     P=L{1}; Q=L{2}; dP=L{3}; dQ=L{4};

//     lambda_int=anm.*dQ;
//     lambda_surf=anm./2.*(P.*dQ+dP.*Q);
//     lambda_ext=anm.*dP;   
// end


