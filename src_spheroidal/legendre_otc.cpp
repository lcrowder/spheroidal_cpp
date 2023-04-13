#include <spheroidal/legendre_otc.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>  
#include <gsl_wrapper/utils.hpp>
#include <gsl_wrapper/core.h>
#include <fmt/core.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_sf_gamma.h>
using namespace std;

// Lorentz's Algorithm to compute continued fraction
gsl::vector cont_frac(int n, int m, gsl::vector u)
{
    ofstream logfile("cont_frac.log");
    
    double tol=1.e-15;
    int max_iters=100000;
    double tiny=1.e-300;
    double a, delta, b, f0, f1, c0, c1, d0, d1;

    int N=u.size();
    gsl::vector H(N);

    // logfile << "Length of u array is " << N << endl;
    
    // Loop over each u value, since each may need different number of iterations to converge
    for(int i=0 ; i<N; ++i)
    {
        double ui=u(i);

        // logfile << "u[" << i << "]=" << ui << endl; 

        c0=tiny;
        f0=tiny;
        d0=0.;

        for(int k=1 ; k<=max_iters; ++k)
        {
            a=-(1.*(n+k-1+m))/(n+k-m);
            b=-(2*(n+k-1)+1)*ui/(n+k-m);
            d1=b+a*d0;
            if (d1==0.) {d1=tiny;}
            c1=b+a/c0;
            if (c1==0.) {c1=tiny;}
            d1=1./d1;
            delta=c1*d1;
            f1=delta*f0;

            if (fabs(delta-1.)<tol)
            {
                // logfile << "continued fraction algorithm converged in " << k << " iterations" <<endl;
                break;
            }
            
            if (k==max_iters)
            {
                logfile << "continued fraction algorithm reached max iterations." <<endl;
            }

            // Update c,d,f values
            c0=c1;
            d0=d1;
            f0=f1;

        }
        H(i)=f1;
        // logfile << "H["<< i << "]=" << H[i] << endl;
    }
    logfile.close();
    return H;
}

int geti(int n,int m)
{
    return m+n*(n+1);
}


void legendre_otc(int p, gsl::vector u, gsl::matrix P)
{
    ofstream logfile("legendre_otc.log");
    
    // Initialize some stuff
    int mm_index, neg_mm_index, mp1m_index, neg_mp1m_index;  
    int nm_index, nm1m_index, nm2m_index, neg_nm_index;
    int np1m_index, np2m_index, n, m;
    double mm_coef, neg_mm_coef, neg_nm_coef; 
    
    int N=u.size();
    
    // Check that all u-values are valid: u>1
    for (int j=0; j<N; ++j)
    {
        if (u(j)<=1.) {cerr << "all u-values must be strictly greater than 1."; break; }
    }
    
    int sp=(p+1)*(p+1); //Number of P and/or Q functions to compute

    // Initialize P, Q, dP, dQ as 2D vectors
    P.resize(sp, N);

    // Vectors of n and m values corresponding to each index.
    vector<int> nn(sp);
    vector<int> mm(sp);
    for (int i=0; i<sp; ++i)
    {
        nn[i]=sqrt(floor(i));
        mm[i]=i-nn[i]*(nn[i]+1);
        // logfile << "i=" << i << endl;
        // logfile << "n=" << nn[i] << ", m=" << mm[i] << endl;
    }

    // P_0^0(u)=1 
    for (int j=0; j< N; ++j){P.set(0, j, 1.);} 
    logfile << "P_0^0 computed." << endl;

    //populate with starting values, P_m^m (m=1:p) and P_{m+1}^m (m=1:p-1)
    if (p > 0)
    {
        for (int m=0; m<=p; ++m)
        {
            // logfile << "this m=" << m << endl;

            if (m>0)
            {
                mm_index = geti(m,m);
                // logfile << "found mm index: " << mm_index << endl;

                mm_coef = gsl_sf_doublefact(2.*m-1.);
                // logfile << "found (2m-1)!!" << endl;  

                neg_mm_index=geti(m,-m);
                // logfile << "found negative mm index: " << neg_mm_index << endl;

                neg_mm_coef=1./gsl_sf_fact(2.*m);
                // logfile << "found negative mm coefficient" << endl;

                for (int j=0; j< N; ++j)
                {
                    // logfile << "u^2-1=" << u[j]*u[j]-1. << endl;
                    // logfile << "exponent = " << 0.5*m << endl;
                    // logfile << "(u^2-1)^(m/2)=" << pow(u[j]*u[j]-1. , 0.5*m) << endl;
                    // logfile << "C(u^2-1)^(m/2)=" << mm_coef*pow(u[j]*u[j]-1. , 0.5*m) << endl;

                    P.set(mm_index, j, mm_coef*pow(u(j)*u(j)-1. , 0.5*m));     // P_m^m
                    // logfile << "found P_m^m[" << j <<"]" << endl;

                    P.set(neg_mm_index, j, neg_mm_coef * P.get(mm_index,j));      //P_m^{-m}
                    // logfile << "found P_m^{-m}[" << j <<"]" << endl;
                }
                
                logfile << "P_{" << m << "}^{" <<  m << "} computed." << endl;
                logfile << "P_{" << m << "}^{" << -m << "} computed." << endl;
                
                
            }
            
            if (m < p)
            {
                mm_index = geti(m,m);
                mp1m_index=geti(m+1,m);
                // logfile << "mp1m_index=" << mp1m_index << endl;
            
                neg_mp1m_index=geti(m+1,-m);
                // logfile << "neg_mp1m_index=" << neg_mp1m_index << endl;

                for (int j=0; j<N; ++j)
                {
                    // logfile << "P[2][" << j << "]=" << P[2][j] << endl;

                    P.set(mp1m_index, j, (2*m+1)* u(j)* P.get(mm_index,j));       //P_{m+1}^m
                    // logfile << "P[(m+1,m)[" << j << "]=" << P[mp1m_index][j] << endl;
                    
                    P.set(neg_mp1m_index, j, 1./gsl_sf_fact(2.*m+1.) *P.get(mp1m_index, j));    //P_{m+1}^{-m}
                    // logfile << "P[(m+1,-m][" << j << "]=" << P[neg_mp1m_index][j] << endl;
                }

                logfile << "P_{" << m+1 << "}^{" << m << "} computed." << endl;
                logfile << "P_{" << m+1 << "}^{" << -m << "} computed." << endl;
            }
        }
    }

    //Use recursion if higher order is needed.
    if (p > 1)
    {   
        for (int m=0; m<=p-2; ++m)
        {
            //Use forward recursion to compute more P_n^m.
            for (int n=m+2; n<=p; ++n)
            {
                nm_index=geti(n,m);
                nm1m_index=geti(n-1,m);
                nm2m_index=geti(n-2,m);
                neg_nm_index=geti(n,-m);
                neg_nm_coef=(gsl_sf_fact(1.*(n-m)))/gsl_sf_fact(1.*(n+m));
                
                for (int j=0; j< N; ++j)
                {
                    P.set(nm_index, j, ( (2*n-1)*u(j)*P.get(nm1m_index,j)-(n+m-1)*P.get(nm2m_index,j) )/(n-m));
                    P.set(neg_nm_index, j, neg_nm_coef*P.get(nm_index,j));
                }
                logfile << "P_{" << n << "}^{" << m << "} computed." << endl;
                logfile << "P_{" << n << "}^{" << -m << "} computed." << endl;
            }
        }
    }

    logfile << "P complete. " << endl;
    logfile.close();
}








// void legendre_otc(int p, vector<double> u, vector<vector<double> > &P, vector<vector<double> > &Q)
// {
//     ofstream logfile("legendre_otc.log");
    
//     // Initialize some stuff
//     int mm_index, neg_mm_index, mp1m_index, neg_mp1m_index;  
//     int nm_index, nm1m_index, nm2m_index, neg_nm_index;
//     int pp_index, neg_pp_index, pm_index, pm1m_index, neg_pm_index, neg_pm1m_index;
//     int np1m_index, np2m_index, n, m;
//     double mm_coef, Qpm1m_coef, neg_pm_coef, neg_pm1m_coef, neg_mm_coef, neg_nm_coef, neg_pp_coef, Qpp1p_coef; 
//     vector<double> H;
    
//     int N=u.size();
    
//     // Check that all u-values are valid: u>1
//     for (int j=0; j<N; ++j)
//     {
//         if (u[j]<=1.) {cerr << "all u-values must be strictly greater than 1."; break; }
//     }
    
//     int sp=(p+1)*(p+1); //Number of P and/or Q functions to compute

//     // Initialize P, Q, dP, dQ as 2D vectors
//     P.resize( sp, vector<double>(N));
//     Q.resize( sp, vector<double>(N));

//     // Vectors of n and m values corresponding to each index.
//     vector<int> nn(sp);
//     vector<int> mm(sp);
//     for (int i=0; i<sp; ++i)
//     {
//         nn[i]=sqrt(floor(i));
//         mm[i]=i-nn[i]*(nn[i]+1);
//         // logfile << "i=" << i << endl;
//         // logfile << "n=" << nn[i] << ", m=" << mm[i] << endl;
//     }

//     // P_0^0(u)=1 
//     for (int j=0; j< N; ++j){P[0][j]=1.;} 
//     logfile << "P_0^0 computed." << endl;

//     //populate with starting values, P_m^m (m=1:p) and P_{m+1}^m (m=1:p-1)
//     if (p > 0)
//     {
//         for (int m=0; m<=p; ++m)
//         {
//             // logfile << "this m=" << m << endl;

//             if (m>0)
//             {
//                 mm_index = geti(m,m);
//                 // logfile << "found mm index: " << mm_index << endl;

//                 mm_coef = gsl_sf_doublefact(2.*m-1.);
//                 // logfile << "found (2m-1)!!" << endl;  

//                 neg_mm_index=geti(m,-m);
//                 // logfile << "found negative mm index: " << neg_mm_index << endl;

//                 neg_mm_coef=1./gsl_sf_fact(2.*m);
//                 // logfile << "found negative mm coefficient" << endl;

//                 for (int j=0; j< N; ++j)
//                 {
//                     // logfile << "u^2-1=" << u[j]*u[j]-1. << endl;
//                     // logfile << "exponent = " << 0.5*m << endl;
//                     // logfile << "(u^2-1)^(m/2)=" << pow(u[j]*u[j]-1. , 0.5*m) << endl;
//                     // logfile << "C(u^2-1)^(m/2)=" << mm_coef*pow(u[j]*u[j]-1. , 0.5*m) << endl;

//                     P[mm_index][j]=mm_coef*pow(u[j]*u[j]-1. , 0.5*m);     // P_m^m
//                     // logfile << "found P_m^m[" << j <<"]" << endl;

//                     P[neg_mm_index][j]=neg_mm_coef*P[mm_index][j];      //P_m^{-m}
//                     // logfile << "found P_m^{-m}[" << j <<"]" << endl;
//                 }
                
//                 logfile << "P_{" << m << "}^{" <<  m << "} computed." << endl;
//                 logfile << "P_{" << m << "}^{" << -m << "} computed." << endl;
                
                
//             }
            
//             if (m < p)
//             {
//                 mm_index = geti(m,m);
//                 mp1m_index=geti(m+1,m);
//                 // logfile << "mp1m_index=" << mp1m_index << endl;
            
//                 neg_mp1m_index=geti(m+1,-m);
//                 // logfile << "neg_mp1m_index=" << neg_mp1m_index << endl;

//                 for (int j=0; j<N; ++j)
//                 {
//                     // logfile << "P[2][" << j << "]=" << P[2][j] << endl;

//                     P[mp1m_index][j] = (2*m+1)* u[j]* P[mm_index][j];       //P_{m+1}^m
//                     // logfile << "P[(m+1,m)[" << j << "]=" << P[mp1m_index][j] << endl;
                    
//                     P[neg_mp1m_index][j] =1./gsl_sf_fact(2.*m+1.) *P[mp1m_index][j];    //P_{m+1}^{-m}
//                     // logfile << "P[(m+1,-m][" << j << "]=" << P[neg_mp1m_index][j] << endl;
//                 }

//                 logfile << "P_{" << m+1 << "}^{" << m << "} computed." << endl;
//                 logfile << "P_{" << m+1 << "}^{" << -m << "} computed." << endl;
//             }
//         }
//     }

//     //Use recursion if higher order is needed.
//     if (p > 1)
//     {   
//         for (int m=0; m<=p-2; ++m)
//         {
//             //Use forward recursion to compute more P_n^m.
//             for (int n=m+2; n<=p; ++n)
//             {
//                 nm_index=geti(n,m);
//                 nm1m_index=geti(n-1,m);
//                 nm2m_index=geti(n-2,m);
//                 neg_nm_index=geti(n,-m);
//                 neg_nm_coef=(gsl_sf_fact(1.*(n-m)))/gsl_sf_fact(1.*(n+m));
                
//                 for (int j=0; j< N; ++j)
//                 {
//                     P[nm_index][j]=( (2*n-1)*u[j]*P[nm1m_index][j]-(n+m-1)*P[nm2m_index][j] )/(n-m);
//                     P[neg_nm_index][j]=neg_nm_coef*P[nm_index][j];
//                 }
//                 logfile << "P_{" << n << "}^{" << m << "} computed." << endl;
//                 logfile << "P_{" << n << "}^{" << -m << "} computed." << endl;
//             }
//         }
//     }

//     logfile << "P complete. " << endl;

//     //--------------------------------------------------------------------------------------
//     // calculate Q:

//     // Calculate highest order separately
//     H = cont_frac(p+1,p,u);
//     pp_index=geti(p,p);
//     neg_pp_index=geti(p,-p);
//     neg_pp_coef=1./gsl_sf_fact(2.*p);
//     Qpp1p_coef=gsl_sf_fact(2.*p)*pow(-1,p);
    
//     for (int j=0; j<N; ++j)
//     {
//         Q[pp_index][j]= Qpp1p_coef/ (((2*p+1)*u[j]-H[j])*P[pp_index][j]);     //Q_p^p
//         Q[neg_pp_index][j]=neg_pp_coef * Q[pp_index][j];                           //Q_p^{-p}
//     }

//     logfile << "Q_{" << p << "}^{" <<  p << "} computed." << endl;
//     logfile << "Q_{" << p << "}^{" << -p << "} computed." << endl;

//     if (p > 0)
//     {

//         for (int m=0; m<=p-1; ++m)
//         {

//             logfile << "before cont_frac: m=" << m <<endl;
//             H=cont_frac(p,m,u);
//             logfile << "after cont_frac: m=" << m <<endl;
//             pm_index=geti(p,m);
//             logfile << "("<<p << "," << m <<") index = " <<pm_index << endl;
//             pm1m_index=geti(p-1,m);
//             logfile << "("<<p-1 << "," << m <<") index = " <<pm1m_index << endl;
//             neg_pm_index=geti(p,-m);
//             logfile << "("<<p << "," << -m <<") index = " << neg_pm_index << endl;
//             neg_pm1m_index=geti(p-1,-m);
//             logfile << "("<<p-1 << "," << -m <<") index = " <<neg_pm1m_index << endl;

//             Qpm1m_coef=(1.*gsl_sf_fact(p+m-1.))/gsl_sf_fact(p-m+0.)*pow(-1,m);
//             neg_pm_coef=(1.*gsl_sf_fact(p-m+0.))/gsl_sf_fact(p+m+0.);
//             neg_pm1m_coef=(1.*gsl_sf_fact(p-1.-m))/gsl_sf_fact(p-1.+m);

//             logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;

//             for (int j=0; j<N; ++j)
//             {
//                 logfile << "j=" << j << endl;

//                 logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
                
//                 Q[pm1m_index][j]=Qpm1m_coef / (P[pm_index][j]-H[j]*P[pm1m_index][j]);

//                 logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
                
//                 Q[pm_index][j]=H[j]*Q[pm1m_index][j];

//                 logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;

//                 if (m>0)
//                 {
//                     Q[neg_pm1m_index][j]=neg_pm1m_coef*Q[pm1m_index][j];
//                     Q[neg_pm_index][j]=neg_pm_coef*Q[pm_index][j];
//                 }

//                 logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
//             }

//             logfile << "Q_{" << p << "}^{" <<  m << "} computed." << endl;
//             logfile << "Q_{" << p-1 << "}^{" <<  m << "} computed." << endl;
//             logfile << "Q_{" << p << "}^{" <<  -m << "} computed." << endl;
//             logfile << "Q_{" << p-1 << "}^{" <<  -m << "} computed." << endl;

//             logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;

//             // Use recursion to compute remaining Qnm's
//             if (p > 1)
//             {

//                 //Use backward recursion to compute Q_n^m for n=p-2:m
//                 for (int n=p-2; n>=m; --n)
//                 {
//                     nm_index=geti(n,m);
//                     np1m_index=geti(n+1,m);
//                     np2m_index=geti(n+2,m);
//                     neg_nm_index=geti(n,-m);
//                     neg_nm_coef=(1.*gsl_sf_fact(n-m+0.))/gsl_sf_fact(n+m+0.);

//                     logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;

//                     for (int j=0; j<N; ++j)
//                     {
//                         Q[nm_index][j]=((2*n+3)*u[j]*Q[np1m_index][j]-(n-m+2)*Q[np2m_index][j])/(n+m+1.);
                        
//                         if (m>0)
//                         {
//                             Q[neg_nm_index][j]=neg_nm_coef * Q[nm_index][j];
//                         }
//                     }
//                     logfile << "Q_{" << n << "}^{" <<  m << "} computed." << endl;
//                     logfile << "Q_{" << n << "}^{" << -m << "} computed." << endl;
//                     logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
//                 }
//             }
//         } 
//     }

//     logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
//     logfile << "Q complete." << endl;
//     logfile.close();
// }






// void Dlegendre_otc(int p, vector<double> u, vector<vector<double> > &dP)
// {
//     ofstream logfile("Dlegendre_otc.log");
    
//     // Initialize some stuff
//     int nm_index, np1m_index, np2m_index, n, m;
//     int N=u.size();
//     int sp=(p+1)*(p+1); //Number of P and/or Q functions to compute

//     vector<vector<double> > P;
//     legendre_otc(p+1,u,&P);

//     // Initialize dPas 2D vector
//     dP.resize( sp, vector<double>(N));

//     // Vectors of n and m values corresponding to each index.
//     vector<int> nn(sp);
//     vector<int> mm(sp);
//     for (int i=0; i<sp; ++i)
//     {
//         nn[i]=sqrt(floor(i));
//         mm[i]=i-nn[i]*(nn[i]+1);
//         // logfile << "i=" << i << endl;
//         // logfile << "n=" << nn[i] << ", m=" << mm[i] << endl;
//     }
    
//     // Compute derivatives
//     logfile << "computing derivatives now..." << endl;

//     for (int k=0; k<sp; ++k)
//     {
//         n=nn[k];
//         m=mm[k];
//         nm_index=geti(n,m);
//         np1m_index=geti(n+1,m);

//         logfile << "n=" << n << endl;
//         logfile << "m=" << m << endl;

//         for (int j=0; j<N; ++j)
//         {
//             dP[nm_index][j]=((m-n-1)*P[np1m_index][j]+(n+1)*u[j]*P[nm_index][j])/(1.-u[j]*u[j]); 
//         }

//         logfile << "dP_{" << n << "}^{" << m <<"} computed." << endl; 
//     }
//     logfile << "derivatives complete." << endl;
// }






// void Dlegendre_otc(int p, vector<double> u, vector<vector<double> > &dP, vector<vector<double> > &dQ)
// {
    
//     ofstream logfile("Dlegendre_otc.log");
    
//     // Initialize some stuff
//     int nm_index, np1m_index, np2m_index, n, m;
//     int N=u.size();
//     int sp=(p+1)*(p+1); //Number of P and/or Q functions to compute

//     vector<vector<double> > P;
//     vector<vector<double> > Q;
//     legendre_otc(p+1,u,&P,&Q);

//     // Initialize dP, dQ as 2D vectors
//     dP.resize( sp, vector<double>(N));
//     dQ.resize( sp, vector<double>(N));

//     // Vectors of n and m values corresponding to each index.
//     vector<int> nn(sp);
//     vector<int> mm(sp);
//     for (int i=0; i<sp; ++i)
//     {
//         nn[i]=sqrt(floor(i));
//         mm[i]=i-nn[i]*(nn[i]+1);
//         // logfile << "i=" << i << endl;
//         // logfile << "n=" << nn[i] << ", m=" << mm[i] << endl;
//     }
    
//     // Compute derivatives
//     logfile << "computing derivatives now..." << endl;

//     for (int k=0; k<sp; ++k)
//     {
//         n=nn[k];
//         m=mm[k];
//         nm_index=geti(n,m);
//         np1m_index=geti(n+1,m);

//         logfile << "n=" << n << endl;
//         logfile << "m=" << m << endl;

//         for (int j=0; j<N; ++j)
//         {
//             dP[nm_index][j]=((m-n-1)*P[np1m_index][j]+(n+1)*u[j]*P[nm_index][j])/(1.-u[j]*u[j]); 
//             dQ[nm_index][j]=((m-n-1)*Q[np1m_index][j]+(n+1)*u[j]*Q[nm_index][j])/(1.-u[j]*u[j]);
//         }

//         logfile << "dP_{" << n << "}^{" << m <<"} computed." << endl; 
//         logfile << "dQ_{" << n << "}^{" << m <<"} computed." << endl; 
//     }
//     logfile << "derivatives complete." << endl;
// }
















// //////////////////////////// //////////////////////////// //////////////////////////// //////////////////////////

// void legendre_otc(int p, vector<double> u, vector<vector<double> > *P, vector<vector<double> > *Q, vector<vector<double> > *dP, vector<vector<double> > *dQ)
// {
    
//     ofstream logfile("legendre_otc.log");
    
//     // Initialize some stuff
//     int mm_index, neg_mm_index, mp1m_index, neg_mp1m_index;  
//     int nm_index, nm1m_index, nm2m_index, neg_nm_index;
//     int pp_index, neg_pp_index, pm_index, pm1m_index, neg_pm_index, neg_pm1m_index;
//     int np1m_index, np2m_index, n, m;
//     double mm_coef, Qpm1m_coef, neg_pm_coef, neg_pm1m_coef, neg_mm_coef, neg_nm_coef, neg_pp_coef, Qpp1p_coef; 
//     vector<double> H;
    
//     int N=u.size();
//     int pmax=p;
    

//     //must compute Q in order to find dQ
//     if (dQoption==1 and Qoption==0){Qoption=1;}
    
//     // Check that all u-values are valid: u>1
//     for (int j=0; j<N; ++j)
//     {
//         if (u[j]<=1.) {cerr << "all u-values must be strictly greater than 1."; break; }
//     }
    
//     // To compute derivatives, we need one higher order of P, Q (p+1) for recurrence relation
//     int sp=(p+1)*(p+1); //Number of dP and/or dQ functions to compute
//     if (dPoption == 1 or dQoption == 1){pmax+=1;}
//     int sp=(pmax+1)*(pmax+1); //Number of P and Q functions to compute

//     // Initialize P, Q, dP, dQ as 2D vectors
//     vector<vector<double> >  P( sp, vector<double>(N));
//     vector<vector<double> >  Q( sp, vector<double>(N));
//     vector<vector<double> > dP(dsp, vector<double>(N));
//     vector<vector<double> > dQ(dsp, vector<double>(N));

//     // Vectors of n and m values corresponding to each index.
//     vector<int> nn(sp);
//     vector<int> mm(sp);
//     for (int i=0; i<sp; ++i)
//     {
//         nn[i]=sqrt(floor(i));
//         mm[i]=i-nn[i]*(nn[i]+1);
//         // logfile << "i=" << i << endl;
//         // logfile << "n=" << nn[i] << ", m=" << mm[i] << endl;
//     }

//     // P_0^0(u)=1 
//     for (int j=0; j< N; ++j){P[0][j]=1.;} 
//     logfile << "P_0^0 computed." << endl;

//     //populate with starting values, P_m^m (m=1:p) and P_{m+1}^m (m=1:p-1)
//     if (pmax > 0)
//     {
//         for (int m=0; m<=pmax; ++m)
//         {
//             // logfile << "this m=" << m << endl;

//             if (m>0)
//             {
//                 mm_index = geti(m,m);
//                 // logfile << "found mm index: " << mm_index << endl;

//                 mm_coef = gsl_sf_doublefact(2.*m-1.);
//                 // logfile << "found (2m-1)!!" << endl;  

//                 neg_mm_index=geti(m,-m);
//                 // logfile << "found negative mm index: " << neg_mm_index << endl;

//                 neg_mm_coef=1./gsl_sf_fact(2.*m);
//                 // logfile << "found negative mm coefficient" << endl;

//                 for (int j=0; j< N; ++j)
//                 {
//                     // logfile << "u^2-1=" << u[j]*u[j]-1. << endl;
//                     // logfile << "exponent = " << 0.5*m << endl;
//                     // logfile << "(u^2-1)^(m/2)=" << pow(u[j]*u[j]-1. , 0.5*m) << endl;
//                     // logfile << "C(u^2-1)^(m/2)=" << mm_coef*pow(u[j]*u[j]-1. , 0.5*m) << endl;

//                     P[mm_index][j]=mm_coef*pow(u[j]*u[j]-1. , 0.5*m);     // P_m^m
//                     // logfile << "found P_m^m[" << j <<"]" << endl;

//                     P[neg_mm_index][j]=neg_mm_coef*P[mm_index][j];      //P_m^{-m}
//                     // logfile << "found P_m^{-m}[" << j <<"]" << endl;
//                 }
                
//                 logfile << "P_{" << m << "}^{" <<  m << "} computed." << endl;
//                 logfile << "P_{" << m << "}^{" << -m << "} computed." << endl;
                
                
//             }
            
//             if (m < pmax)
//             {
//                 mm_index = geti(m,m);
//                 mp1m_index=geti(m+1,m);
//                 // logfile << "mp1m_index=" << mp1m_index << endl;
            
//                 neg_mp1m_index=geti(m+1,-m);
//                 // logfile << "neg_mp1m_index=" << neg_mp1m_index << endl;

//                 for (int j=0; j<N; ++j)
//                 {
//                     // logfile << "P[2][" << j << "]=" << P[2][j] << endl;

//                     P[mp1m_index][j] = (2*m+1)* u[j]* P[mm_index][j];       //P_{m+1}^m
//                     // logfile << "P[(m+1,m)[" << j << "]=" << P[mp1m_index][j] << endl;
                    
//                     P[neg_mp1m_index][j] =1./gsl_sf_fact(2.*m+1.) *P[mp1m_index][j];    //P_{m+1}^{-m}
//                     // logfile << "P[(m+1,-m][" << j << "]=" << P[neg_mp1m_index][j] << endl;
//                 }

//                 logfile << "P_{" << m+1 << "}^{" << m << "} computed." << endl;
//                 logfile << "P_{" << m+1 << "}^{" << -m << "} computed." << endl;
//             }
//         }
//     }

//     //Use recursion if higher order is needed.
//     if (pmax > 1)
//     {   
//         for (int m=0; m<=pmax-2; ++m)
//         {
//             //Use forward recursion to compute more P_n^m.
//             for (int n=m+2; n<=pmax; ++n)
//             {
//                 nm_index=geti(n,m);
//                 nm1m_index=geti(n-1,m);
//                 nm2m_index=geti(n-2,m);
//                 neg_nm_index=geti(n,-m);
//                 neg_nm_coef=(gsl_sf_fact(1.*(n-m)))/gsl_sf_fact(1.*(n+m));
                
//                 for (int j=0; j< N; ++j)
//                 {
//                     P[nm_index][j]=( (2*n-1)*u[j]*P[nm1m_index][j]-(n+m-1)*P[nm2m_index][j] )/(n-m);
//                     P[neg_nm_index][j]=neg_nm_coef*P[nm_index][j];
//                 }
//                 logfile << "P_{" << n << "}^{" << m << "} computed." << endl;
//                 logfile << "P_{" << n << "}^{" << -m << "} computed." << endl;
//             }
//         }
//     }

//     logfile << "P complete. " << endl;

//     //--------------------------------------------------------------------------------------
//     // If desired, calculate Q:
//     if (Qoption==1)
//     {
//         // Calculate highest order separately
//         H = cont_frac(pmax+1,pmax,u);
//         pp_index=geti(pmax,pmax);
//         neg_pp_index=geti(pmax,-pmax);
//         neg_pp_coef=1./gsl_sf_fact(2.*pmax);
//         Qpp1p_coef=gsl_sf_fact(2.*pmax)*pow(-1,pmax);
        
//         for (int j=0; j<N; ++j)
//         {
//             Q[pp_index][j]= Qpp1p_coef/ (((2*pmax+1)*u[j]-H[j])*P[pp_index][j]);     //Q_p^p
//             Q[neg_pp_index][j]=neg_pp_coef * Q[pp_index][j];                           //Q_p^{-p}
//         }

//         logfile << "Q_{" << pmax << "}^{" <<  pmax << "} computed." << endl;
//         logfile << "Q_{" << pmax << "}^{" << -pmax << "} computed." << endl;

//         if (pmax > 0)
//         {

//             for (int m=0; m<=pmax-1; ++m)
//             {

//                 logfile << "before cont_frac: m=" << m <<endl;
//                 H=cont_frac(pmax,m,u);
//                 logfile << "after cont_frac: m=" << m <<endl;
//                 pm_index=geti(pmax,m);
//                 logfile << "("<<pmax << "," << m <<") index = " <<pm_index << endl;
//                 pm1m_index=geti(pmax-1,m);
//                 logfile << "("<<pmax-1 << "," << m <<") index = " <<pm1m_index << endl;
//                 neg_pm_index=geti(pmax,-m);
//                 logfile << "("<<pmax << "," << -m <<") index = " << neg_pm_index << endl;
//                 neg_pm1m_index=geti(pmax-1,-m);
//                 logfile << "("<<pmax-1 << "," << -m <<") index = " <<neg_pm1m_index << endl;

//                 Qpm1m_coef=(1.*gsl_sf_fact(pmax+m-1.))/gsl_sf_fact(pmax-m+0.)*pow(-1,m);
//                 neg_pm_coef=(1.*gsl_sf_fact(pmax-m+0.))/gsl_sf_fact(pmax+m+0.);
//                 neg_pm1m_coef=(1.*gsl_sf_fact(pmax-1.-m))/gsl_sf_fact(pmax-1.+m);

//                 logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;

//                 for (int j=0; j<N; ++j)
//                 {
//                     logfile << "j=" << j << endl;

//                     logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
                    
//                     Q[pm1m_index][j]=Qpm1m_coef / (P[pm_index][j]-H[j]*P[pm1m_index][j]);

//                     logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
                    
//                     Q[pm_index][j]=H[j]*Q[pm1m_index][j];

//                     logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;

//                     if (m>0)
//                     {
//                         Q[neg_pm1m_index][j]=neg_pm1m_coef*Q[pm1m_index][j];
//                         Q[neg_pm_index][j]=neg_pm_coef*Q[pm_index][j];
//                     }

//                     logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
//                 }

//                 logfile << "Q_{" << pmax << "}^{" <<  m << "} computed." << endl;
//                 logfile << "Q_{" << pmax-1 << "}^{" <<  m << "} computed." << endl;
//                 logfile << "Q_{" << pmax << "}^{" <<  -m << "} computed." << endl;
//                 logfile << "Q_{" << pmax-1 << "}^{" <<  -m << "} computed." << endl;

//                 logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;

//                 // Use recursion to compute remaining Qnm's
//                 if (pmax > 1)
//                 {

//                     //Use backward recursion to compute Q_n^m for n=p-2:m
//                     for (int n=pmax-2; n>=m; --n)
//                     {
//                         nm_index=geti(n,m);
//                         np1m_index=geti(n+1,m);
//                         np2m_index=geti(n+2,m);
//                         neg_nm_index=geti(n,-m);
//                         neg_nm_coef=(1.*gsl_sf_fact(n-m+0.))/gsl_sf_fact(n+m+0.);

//                         logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;

//                         for (int j=0; j<N; ++j)
//                         {
//                             Q[nm_index][j]=((2*n+3)*u[j]*Q[np1m_index][j]-(n-m+2)*Q[np2m_index][j])/(n+m+1.);
                            
//                             if (m>0)
//                             {
//                                 Q[neg_nm_index][j]=neg_nm_coef * Q[nm_index][j];
//                             }
//                         }
//                         logfile << "Q_{" << n << "}^{" <<  m << "} computed." << endl;
//                         logfile << "Q_{" << n << "}^{" << -m << "} computed." << endl;
//                         logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
//                     }
//                 }
//             }
//         } 
//         logfile << "Q complete." << endl;
//     }

//     logfile << "Q_0^0(" << u[0] << ") = " << Q[0][0] << endl;
    
//     //--------------------------------------------------------------------------------------

//     // Compute derivatives, if desired
//     if (dPoption==1 or dQoption==1)
//     {
//         logfile << "computing derivatives now..." << endl;

//         for (int k=0; k<dsp; ++k)
//         {
//             n=nn[k];
//             m=mm[k];
//             nm_index=geti(n,m);
//             np1m_index=geti(n+1,m);

//             logfile << "n=" << n << endl;
//             logfile << "m=" << m << endl;

//             for (int j=0; j<N; ++j)
//             {
//                 if (dPoption==1) { dP[nm_index][j]=((m-n-1)*P[np1m_index][j]+(n+1)*u[j]*P[nm_index][j])/(1.-u[j]*u[j]); }
//                 if (dQoption==1) { dQ[nm_index][j]=((m-n-1)*Q[np1m_index][j]+(n+1)*u[j]*Q[nm_index][j])/(1.-u[j]*u[j]); }
//             }

//             if (dPoption==1) { logfile << "dP_{" << n << "}^{" << m <<"} computed." << endl; }
//             if (dQoption==1) { logfile << "dQ_{" << n << "}^{" << m <<"} computed." << endl; }
//         }
//     }

//     // Store all arrays into 3D matrix
//     int num_matrices=1+Qoption+dPoption+dQoption;

//     // logfile << "num matrices = " << num_matrices << endl;

//     vector<vector<vector<double> > > PQdPdQ(num_matrices, vector<vector<double> >(dsp, vector<double>(N)));

    
//     // logfile << "PQdPdQ dim 1 size =" << PQdPdQ.size() << endl;
//     // logfile << "PQdPdQ dim 2 size =" << PQdPdQ[0].size() << endl;
//     // logfile << "PQdPdQ dim 3 size =" << PQdPdQ[0][0].size() << endl;

//     for (int i=0; i<dsp; ++i )
//     {
//         for (int j=0; j<N; ++j)
//         {
//             PQdPdQ[0][i][j]=P[i][j];
//             if (Qoption==1) {PQdPdQ[1][i][j]=Q[i][j]; }
//             if (dPoption==1) {PQdPdQ[1+Qoption][i][j]=dP[i][j]; }
//             if (dQoption==1) {PQdPdQ[2+dPoption][i][j]=dQ[i][j]; }
//         }
//     }

    
//     logfile.close();
//     return PQdPdQ;
    
// }