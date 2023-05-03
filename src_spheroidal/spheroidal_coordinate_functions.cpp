#include <spheroidal/spheroidal_coordinate_functions.h>
#include <spheroidal/grid_functions.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>  
#include <math.h> 
#include <yawg/utils.hpp>
#include <yawg/core.h>
using namespace std;


gsl::matrix spheroidal_to_cart(gsl::matrix S, double a)
{
    int s1 = S.nrows();
    int s2 = S.ncols();
    if (floor(s2/3) != s2/3.)
    {
        cerr << "Given points are not correct dimensions. #columns must be a multiple of 3." << endl;
    }

    gsl::matrix X(s1,s2);

    for (int j=0; j<s2/3; j++)
    {
        gsl::vector u = S.column(3*j);
        gsl::vector v = S.column(3*j+1);
        gsl::vector phi= S.column(3*j+2);

        // printf("extracted u,v,phi\n");

        gsl::vector x(s1),y(s1),z(s1);
        for (int k=0; k<s1; k++)
        {
            x(k) = a*sqrt(pow(u(k),2)-1)*sqrt(1-pow(v(k),2))*cos(phi(k));
            y(k) = a*sqrt(pow(u(k),2)-1)*sqrt(1-pow(v(k),2))*sin(phi(k));
            z(k) = a*u(k)*v(k);
        }

        // printf("computed cos, sin, sqrt\n");
                                                           
        X.column(3*j) = x;
        X.column(3*j+1) = y;
        X.column(3*j+2) = z;

        // printf("computed X\n");
    }
    return X;
}

gsl::matrix cart_to_spheroidal(gsl::matrix X, double a)
{
    int s1 = X.nrows();
    int s2 = X.ncols();
    if (floor(s2/3) != s2/3.)
    {
        cerr << "Given points are not correct dimensions. #columns must be a multiple of 3." << endl;
    }

    gsl::matrix S(s1,s2);

    for (int j=0; j<s2/3; j++)
    {
        gsl::vector x = X.column(3*j);
        gsl::vector y = X.column(3*j+1);
        gsl::vector z= X.column(3*j+2);

        printf("extracted x,y,z\n");

        gsl::vector u(s1),v(s1),phi(s1);
        for (int k=0; k<s1; k++)
        {
            u(k) = (sqrt(pow(x(k),2) + pow(y(k),2) + pow((z(k)+a),2)) + sqrt(pow(x(k),2) + pow(y(k),2) + pow((z(k)-a),2)))/(2.*a);
            v(k) = (sqrt(pow(x(k),2) + pow(y(k),2) + pow(z(k)+a,2)) - sqrt(pow(x(k),2) + pow(y(k),2) + pow(z(k)-a,2)))/(2.*a);
            phi(k) = atan2(y(k),x(k));
        }

        // printf("computed sqrt, atan2\n");
                                                           
        S.column(3*j) = u;
        S.column(3*j+1) = v;
        S.column(3*j+2) = phi;

        // printf("computed S\n");
    }
    return S;
}

