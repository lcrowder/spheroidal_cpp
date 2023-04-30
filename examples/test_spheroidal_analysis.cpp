#include <spheroidal/spheroidal_harmonic_transforms.h>
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
using namespace std;


int main(int argc, char *argv[])
{
    // int p=1;
    // int sp=(p+1)*(p+1);
    
    // if (argc==2){p=atoi(argv[1]);}
    
    // printf("p=%d\n", p);

    // gsl::matrix theta, phi;
    // gl_grid(p, theta, phi);
    // theta.print();
    // phi.print();

    // gsl::matrix f;
    // f.resize(2*p*(p+1),2);
    // f.column(0) = gsl::linspace(1. , 2.*p*(p+1), 2*p*(p+1));
    // f.column(1) = gsl::linspace(1. , 2.*p*(p+1), 2*p*(p+1));

    // gsl::cmatrix shc=spheroidal_analysis(f);

    // shc.print();

    // printf("Starting spheroidal synthesis\n");

    // gsl::cmatrix f2=spheroidal_snythesis(shc);

    // f2.print();

    
    // ------------------------------------------------------------
    //coordinate transform tests
    // ------------------------------------------------------------

    // gsl::vector u=gsl::linspace(1.5,1.5,10);
    // gsl::vector v=gsl::linspace(-1,1,10);
    // gsl::vector phi=gsl::linspace(0,2*M_PI,10);

    // gsl::matrix S(10,3);
    // S.column(0)=u;
    // S.column(1)=v;
    // S.column(2)=phi;
    // S.print();

    // gsl::matrix X=spheroidal_to_cart(S,1/1.5);
    // X.print();

    // gsl::matrix S2=cart_to_spheroidal(X,1/1.5);
    // S2.print();

    // gsl::matrix X2 = spheroidal_to_cart(S2,1/1.5);
    // X2-=X;
    // X2.print();

    // ------------------------------------------------------------




    // ------------------------------------------------------------
    // double layer tests
    // ------------------------------------------------------------

    //Target points
    gsl::vector ut=gsl::linspace(1.5,3.5,4);
    gsl::vector vt=gsl::linspace(-1,1,4);
    gsl::vector phit=gsl::linspace(0,2*M_PI,4);
    
    gsl::matrix S(4,3);
    // gsl::matrix S(4,6);
    S.column(0)=ut; S.column(1)=vt; S.column(2)=phit;
    // S.column(3)=2*ut; S.column(4)=vt; S.column(5)=phit;

    S.print();

    int p=1;
    int np = 2*p*(p+1);
    gsl::vector sigma0=gsl::linspace(1., 2.*p*(p+1), np);
    gsl::cmatrix sigma(np,1); sigma.column(0) = sigma0; 
    // sigma.column(1)=sigma0;
    sigma.print();

    double u0=2./sqrt(3);
    int coords=1;

    // gsl::cmatrix DL = spheroidal_double_layer(sigma, u0, S, coords);

    printf("sigma0:\n");
    sigma0.print();
    // gsl::cvector DLvec = spheroidal_double_layer(sigma0, u0, S, coords);
    gsl::cmatrix DL= spheroidal_double_layer(sigma, u0);

    DL.print();
    // DLvec.print();


    // int regime =1;
    // gsl::matrix F=solid_harmonic(p,u,regime);
    // F.print();

    // gsl::vector li,ls,le;
    // DLspectrum(p,1.5,li,ls,le);
    // li.print();
    // ls.print();
    // le.print();




    // ------------------------------------------------------------
 

    // cout << "Last entry of P is " << P.get((p+1)*(p+1)-1,N-1) << endl;
    // cout << "Last entry of Q is " << Q.get((p+1)*(p+1)-1,N-1) << endl;
    // cout << "Last entry of dP is " << dP.get((p+1)*(p+1)-1,N-1) << endl;
    // cout << "Last entry of dQ is " << dQ.get((p+1)*(p+1)-1,N-1) << endl;

    // FILE* Pfile = fopen("../../data/P.csv","w");
    // FILE* Qfile = fopen("../../data/Q.csv","w");
    // FILE* dPfile = fopen("../../data/dP.csv","w");
    // FILE* dQfile = fopen("../../data/dQ.csv","w");

    // P.print_csv(Pfile); 
    // Q.print_csv(Qfile); 
    // dP.print_csv(dPfile); 
    // dQ.print_csv(dQfile);

    // fclose(Pfile);
    // fclose(Qfile);
    // fclose(dPfile);
    // fclose(dQfile);

    return 0;
}

