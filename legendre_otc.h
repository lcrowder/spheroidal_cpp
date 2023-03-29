#include <iostream>
#include <vector>
#include <stdio.h>  
#include <math.h> 
using namespace std;


// Lorentz's Algorithm to compute continued fraction
vector<double> cont_frac(int n, int m, vector<double> u)
{
    double tol=1.e-10;
    int max_iters=100000;
    double tiny=1.e-300;
    double a;
    double delta;
    double b;  
    double f0;
    double f1;  
    double c0; 
    double c1;
    double d0;  
    double d1;

    int N=u.size();
    vector<double> H(N);

    cout << "Length of u array is " << N << endl;
    
    // Loop over each u value, since each may need different number of iterations to converge
    for(int i=0 ; i<N; ++i)
    {
        double ui=u[i];

        cout << "u[" << i << "]=" << ui << endl; 

        c0=tiny;
        f0=tiny;
        d0=0.;

        for(int k=1 ; k<=max_iters; ++k)
        {
            a=-(n+k-1+m)/(n+k-m);
            b=-(2*(n+k-1)+1)*ui/(n+k-m);
            d1=b+a*d0;
            if (d1==0.) {d1=tiny;}
            c1=b+a/c0;
            if (c1==0.) {c1=tiny;}
            d1=1./d1;
            delta=c1*d1;
            f1=delta*f0;

            if (fabs(delta-1)<tol)
            {
                cout << "continued fraction algorithm converged in " << k << " iterations" <<endl;
                break;
            }
            
            if (k==max_iters)
            {
                cout << "continued fraction algorithm reached max iterations." <<endl;
            }

            // Update c,d,f values
            c0=c1;
            d0=d1;
            f0=f1;

        }
        H[i]=f1;
        cout << "H["<< i << "]=" << H[i] << endl;
    }
    return H;
}

int geti(int n,int m)
{
    return m+n*(n+1);
}

int double_factorial(int n)
{
    int prod=1;
    for (int i = n; i >=1; i-=2)
    {
        prod = prod*i;
    }
    return prod;
}

int factorial(int n)
{
    int prod=1;
    for (int i = 1; i <=n; ++i)
    {
        prod = i*prod;
    }
    return prod;
}

void legendre_otc(int p, vector<double> u, int Qoption = 1, int dPoption = 1, int dQoption = 1)
{
    // initialize values needed to construct matrices of correct size.
    // To compute derivatives, we need one higher order of P, Q (p+1) for recurrence relation
    int pmax=p;
    int dsp;
    if (dPoption == 1 or dQoption == 1)
    {
        pmax+=1;
        dsp=(p+1)*(p+1);
    }
    int N=u.size();
    int sp=(pmax+1)*(pmax+1);

    cout << "p=" << p <<endl;
    cout << "pmax=" << pmax <<endl;

    // Vectors of n and m values corresponding to each index.
    vector<int> nn(sp);
    vector<int> mm(sp);
    for (int i=0; i<sp; ++i)
    {
        nn[i]=sqrt(floor(i));
        mm[i]=i-nn[i]*(nn[i]+1);
        // cout << "i=" << i << endl;
        // cout << "n=" << nn[i] << ", m=" << mm[i] << endl;
    }

    // Initialize P, Q, dP, dQ as 2D vectors
    vector<vector<double> > P( sp, vector<double>(N,1));
    for (int j=0; j< N; ++j){P[0][j]=1.;} // P_0^0(u)=1 
    cout << "P_0^0 computed." << endl;

    // if ( Qoption == 1 ){vector<vector<double> >  Q( sp, vector<double>(N,1));}
    // if (dPoption == 1 ){vector<vector<double> > dP(dsp, vector<double>(N,1));}
    // if (dQoption == 1 ){vector<vector<double> > dQ(dsp, vector<double>(N,1));}

    //populate with starting values, P_m^m (m=1:p) and P_{m+1}^m (m=1:p-1)
    if (pmax > 0)
    {
        int mm_index;
        double mm_coef;
        int neg_mm_index;
        double neg_mm_coef;
        int mp1m_index;
        int neg_mp1m_index; 

        for (int m=0; m<=pmax; ++m)
        {
            cout << "this m=" << m << endl;

            if (m>0)
            {
                mm_index = geti(m,m);
                cout << "found mm index: " << mm_index << endl;

                mm_coef = double_factorial(2*m-1);
                cout << "found (2m-1)!!" << endl;  

                neg_mm_index=geti(m,-m);
                cout << "found negative mm index: " << neg_mm_index << endl;

                neg_mm_coef=1.0/factorial(2*m);
                cout << "found negative mm coefficient" << endl;

                for (int j=0; j< N; ++j)
                {
                    // cout << "u^2-1=" << u[j]*u[j]-1. << endl;
                    // cout << "exponent = " << 0.5*m << endl;
                    // cout << "(u^2-1)^(m/2)=" << pow(u[j]*u[j]-1. , 0.5*m) << endl;
                    // cout << "C(u^2-1)^(m/2)=" << mm_coef*pow(u[j]*u[j]-1. , 0.5*m) << endl;

                    P[mm_index][j]=mm_coef*pow(u[j]*u[j]-1. , 0.5*m);     // P_m^m
                    // cout << "found P_m^m[" << j <<"]" << endl;

                    P[neg_mm_index][j]=neg_mm_coef*P[mm_index][j];      //P_m^{-m}
                    // cout << "found P_m^{-m}[" << j <<"]" << endl;
                }
                
                cout << "P_{" << m << "}^{" <<  m << "} computed." << endl;
                cout << "P_{" << m << "}^{" << -m << "} computed." << endl;
                
                
            }
            
            if (m < pmax)
            {
                mm_index = geti(m,m);
                mp1m_index=geti(m+1,m);
                cout << "mp1m_index=" << mp1m_index << endl;
            
                neg_mp1m_index=geti(m+1,-m);
                cout << "neg_mp1m_index=" << neg_mp1m_index << endl;

                for (int j=0; j<N; ++j)
                {
                    // cout << "P[2][" << j << "]=" << P[2][j] << endl;

                    P[mp1m_index][j] = (2*m+1)* u[j]* P[mm_index][j];       //P_{m+1}^m
                    cout << "P[(m+1,m)[" << j << "]=" << P[mp1m_index][j] << endl;
                    
                    P[neg_mp1m_index][j] =1./factorial(2*m+1) *P[mp1m_index][j];    //P_{m+1}^{-m}
                    cout << "P[(m+1,-m][" << j << "]=" << P[neg_mp1m_index][j] << endl;
                }

                cout << "P_{" << m+1 << "}^{" << m << "} computed." << endl;
                cout << "P_{" << m+1 << "}^{" << -m << "} computed." << endl;
            }
        }
    }

    //Use recursion if higher order is needed.
    if (pmax > 1)
    {
        int nm_index;
        int nm1m_index;
        int nm2m_index;
        int neg_nm_index;
        double neg_nm_coef;
        
        for (int m=0; m<=pmax-2; ++m)
        {
            //Use forward recursion to compute more P_n^m.
            for (int n=m+2; n<=pmax; ++n)
            {
                nm_index=geti(n,m);
                nm1m_index=geti(n-1,m);
                nm2m_index=geti(n-2,m);
                neg_nm_index=geti(n,-m);
                neg_nm_coef=(1.*factorial(n-m))/factorial(n+m);
                
                for (int j=0; j< N; ++j)
                {
                    P[nm_index][j]=( (2*n-1)*u[j]*P[nm1m_index][j]-(n+m-1)*P[nm2m_index][j] )/(n-m);
                    P[neg_nm_index][j]=neg_nm_coef*P[nm_index][j];
                }
                cout << "P_{" << n << "}^{" << m << "} computed." << endl;
                cout << "P_{" << n << "}^{" << -m << "} computed." << endl;
            }
        }
    }

    cout << "P complete. " << endl;

    //--------------------------------------------------------------------------------------
    // If desired, calculate Q:
    if (Qoption==1)
    {
        vector<vector<double> >  Q( sp, vector<double>(N,1));

        // Calculate highest order separately
        vector<double> Hp = cont_frac(pmax+1,pmax,u);
        int pp_index=geti(pmax,pmax);
        int neg_pp_index=geti(pmax,-pmax);
        double neg_pp_coef=1./factorial(2*pmax);
        double Qpp1p_coef=factorial(2*pmax)*pow(-1,pmax);

        for (int j=1; j<N; ++j)
        {
            Q[pp_index][j]= Qpp1p_coef/ (((2*pmax+1)* u[j]-Hp[j])*P[pp_index][j]);     //Q_p^p
            Q[neg_pp_index][j]=neg_pp_coef * Q[pp_index][j];                           //Q_p^{-p}
        }

        cout << "Q_{" << pmax << "}^{" <<  pmax << "} computed." << endl;
        cout << "Q_{" << pmax << "}^{" << -pmax << "} computed." << endl;

        if (pmax > 0)
        {
            vector<double> H;
            int pm_index;
            int pm1m_index;
            int neg_pm_index;
            int neg_pm1m_index;
            double Qpm1m_coef;
            double neg_pm_coef;
            double neg_pm1m_coef;

            for (int m=0; m<=pmax-1; ++m)
            {

                H=cont_frac(pmax,m,u);
                pm_index=geti(pmax,m);
                pm1m_index=geti(pmax-1,m);
                neg_pm_index=geti(pmax,-m);
                neg_pm1m_index=geti(pmax-1,-m);
                Qpm1m_coef=(1.*factorial(pmax+m-1))/factorial(pmax-m)*pow(-1,m);
                neg_pm_coef=(1.*factorial(pmax-m))/factorial(pmax+m);
                neg_pm1m_coef=(1.*factorial(pmax-1-m))/factorial(pmax-1+m);

                for (int j=0; j<N; ++j)
                {
                    Q[pm1m_index][j]=Qpm1m_coef / (P[pm_index][j]-H[j]*P[pm1m_index][j]);
                    Q[pm_index][j]=H[j]*Q[pm1m_index][j];

                    Q[neg_pm1m_index][j]=neg_pm1m_coef*Q[pm1m_index][j];
                    Q[neg_pm_coef][j]=neg_pm_coef*Q[pm_index][j];
                }

                cout << "Q_{" << p << "}^{" <<  m << "} computed." << endl;
                cout << "Q_{" << p-1 << "}^{" <<  m << "} computed." << endl;
                cout << "Q_{" << p << "}^{" <<  -m << "} computed." << endl;
                cout << "Q_{" << p-1 << "}^{" <<  -m << "} computed." << endl;

                // Use recursion to compute remaining Qnm's
                if (pmax > 1)
                {
                    int nm_index;
                    int np1m_index;
                    int np2m_index;
                    int neg_nm_index;
                    double neg_nm_coef;

                    //Use backward recursion to compute Q_n^m for n=p-2:m
                    for (int n=pmax-2; n<=m; --n)
                    {
                        nm_index=geti(n,m);
                        np1m_index=geti(n+1,m);
                        np2m_index=geti(n+2,m);
                        neg_nm_index=geti(n,-m);
                        neg_nm_coef=(1.*factorial(n-m))/factorial(n+m);

                        for (int j=0; j<N; ++j)
                        {
                            Q[nm_index][j]=((2*n+3)*u[j]*Q[np1m_index][j]-(n-m+2)*Q[np2m_index][j])/(n+m+1);
                            Q[neg_nm_index][j]=neg_nm_coef * Q[nm_index][j];
                        }
                        cout << "Q_{" << n << "}^{" <<  m << "} computed." << endl;
                        cout << "Q_{" << n << "}^{" << -m << "} computed." << endl;
                    }
                }
            }
        } 
    }
    cout << "Q complete." << endl;
    //--------------------------------------------------------------------------------------

}

