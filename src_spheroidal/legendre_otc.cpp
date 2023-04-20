#include <spheroidal/legendre_otc.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <yawg/utils.hpp>
#include <yawg/core.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_sf_gamma.h>
using namespace std;

// Lorentz's Algorithm to compute continued fraction
gsl::vector cont_frac(int n, int m, gsl::vector u)
{
    FILE *logfile = fopen("../../cont_frac.log", "w");

    double tol = 1.e-15;
    int max_iters = 100000;
    double tiny = 1.e-300;
    double a, delta, b, f0, f1, c0, c1, d0, d1;

    int N = u.size();
    gsl::vector H(N);

    fprintf(logfile, "Length of u array is %d\n", N);

    // Loop over each u value, since each may need different number of iterations to converge
    for (int i = 0; i < N; ++i)
    {
        double ui = u(i);

        fprintf(logfile,"u[%d]=%.4f\n", i, ui );

        c0 = tiny;
        f0 = tiny;
        d0 = 0.;

        for (int k = 1; k <= max_iters; ++k)
        {
            a = -(1. * (n + k - 1 + m)) / (n + k - m);
            b = -(2 * (n + k - 1) + 1) * ui / (n + k - m);
            d1 = b + a * d0;
            if (d1 == 0.)
            {
                d1 = tiny;
            }
            c1 = b + a / c0;
            if (c1 == 0.)
            {
                c1 = tiny;
            }
            d1 = 1. / d1;
            delta = c1 * d1;
            f1 = delta * f0;

            if (fabs(delta - 1.) < tol)
            {
                fprintf(logfile , "continued fraction algorithm converged in %d iterations\n", k);
                break;
            }

            if (k == max_iters)
            {
                fprintf(logfile, "continued fraction algorithm reached max iterations.\n");
            }

            // Update c,d,f values
            c0 = c1;
            d0 = d1;
            f0 = f1;
        }
        H(i) = f1;
        fprintf(logfile, "H[%d]=%.5f\n", i, H(i));
    }
    fclose(logfile);
    return H;
}

int geti(int n, int m)
{
    return m + n * (n + 1);
}



/*!
 * \brief Computes associated legendre functions off the cut of the first kind, $$P_n^m(u)$$, $$u>1$$.
 * \param p The order of the spheroidal harmonics
 * \param u The u-values at which to compute the associated legendre functions
 * \param P Matrix of associated legendre functions off the cut of the first kind. Each row corresponds to a different n and m value, and each column corresponds to a different u-value.
 *
 * \note Will cause an error if any u-value is less than or equal to 1.
 */
void legendre_otc(int p, gsl::vector u, gsl::matrix &P)
{

    // Open a log file to print status updates of computations
    FILE *logfile = fopen("../../legendre_otc.log", "w");
    fprintf(logfile, "opened log file.\n");

    // Initialize index and coefficient variables
    int mm_index, neg_mm_index, mp1m_index, neg_mp1m_index;
    int nm_index, nm1m_index, nm2m_index, neg_nm_index;
    int np1m_index, np2m_index, n, m;
    double mm_coef, neg_mm_coef, neg_nm_coef;
    int N = u.size();
    int sp = (p + 1) * (p + 1); // Number of P and/or Q functions to compute

    // Check that all u-values are valid: u>1
    for (int j = 0; j < N; ++j)
    {
        if (u(j) <= 1.)
        {
            cerr << "all u-values must be strictly greater than 1.";
            break;
        }
    }

    P.resize(sp, N);

    // Vectors of n and m values corresponding to each index.
    vector<int> nn(sp);
    vector<int> mm(sp);
    for (int i = 0; i < sp; ++i)
    {
        nn[i] = sqrt(floor(i));
        mm[i] = i - nn[i] * (nn[i] + 1);
    }

    // P_0^0(u)=1
    for (int j = 0; j < N; ++j)
    {
        P.set(0, j, 1.);
    }
    fprintf(logfile, "P_0^0 computed.\n");

    // populate with starting values, P_m^m (m=1:p) and P_{m+1}^m (m=1:p-1)
    if (p > 0)
    {
        for (int m = 0; m <= p; ++m)
        {
            if (m > 0)
            {
                mm_index = geti(m, m);
                mm_coef = gsl_sf_doublefact(2. * m - 1.);
                neg_mm_index = geti(m, -m);
                neg_mm_coef = 1. / gsl_sf_fact(2. * m);

                for (int j = 0; j < N; ++j)
                {
                    P.set(mm_index, j, mm_coef * pow(u(j) * u(j) - 1., 0.5 * m)); // P_m^m
                    P.set(neg_mm_index, j, neg_mm_coef * P.get(mm_index, j)); // P_m^{-m}
                }

                fprintf(logfile, "P_{%d}^{%d} computed.\n",m,m );
                fprintf(logfile, "P_{%d}^{%d} computed.\n",m,-m );
            }

            if (m < p)
            {
                mm_index = geti(m, m);
                mp1m_index = geti(m + 1, m);
                neg_mp1m_index = geti(m + 1, -m);

                for (int j = 0; j < N; ++j)
                {
                    P.set(mp1m_index, j, (2 * m + 1) * u(j) * P.get(mm_index, j)); // P_{m+1}^m
                    P.set(neg_mp1m_index, j, 1. / gsl_sf_fact(2. * m + 1.) * P.get(mp1m_index, j)); // P_{m+1}^{-m}
                }

                fprintf(logfile, "P_{%d}^{%d} computed.\n", m+1, m );
                fprintf(logfile, "P_{%d}^{%d} computed.\n", m+1, -m );
            }
        }
    }

    // Use recursion if higher order is needed.
    if (p > 1)
    {
        for (int m = 0; m <= p - 2; ++m)
        {
            // Use forward recursion to compute more P_n^m.
            for (int n = m + 2; n <= p; ++n)
            {
                nm_index = geti(n, m);
                nm1m_index = geti(n - 1, m);
                nm2m_index = geti(n - 2, m);
                neg_nm_index = geti(n, -m);
                neg_nm_coef = (gsl_sf_fact(1. * (n - m))) / gsl_sf_fact(1. * (n + m));

                for (int j = 0; j < N; ++j)
                {
                    P.set(nm_index, j, ((2 * n - 1) * u(j) * P.get(nm1m_index, j) - (n + m - 1) * P.get(nm2m_index, j)) / (n - m));
                    P.set(neg_nm_index, j, neg_nm_coef * P.get(nm_index, j));
                }
                fprintf(logfile, "P_{%d}^{%d} computed.\n",n,m );
                fprintf(logfile, "P_{%d}^{%d} computed.\n",n,-m );
            }
        }
    }

    fprintf(logfile, "P complete. \n" );
    fclose(logfile);
}



/*!
 * \brief Computes associated legendre functions off the cut of the first and second kind, $$P_n^m(u)$$ and $$Q_n^m(u)$$, for $$u>1$$.
 * \param p The order of the spheroidal harmonics
 * \param u The u-values at which to compute the associated legendre functions
 * \param P Matrix of associated legendre functions off the cut of the first kind. Each row corresponds to a different n and m value, and each column corresponds to a different u-value.
 * \param Q Matrix of associated legendre functions off the cut of the second kind. Each row corresponds to a different n and m value, and each column corresponds to a different u-value.
 * \note Will cause an error if any u-value is less than or equal to 1.
 */
void legendre_otc(int p, gsl::vector u, gsl::matrix &P, gsl::matrix &Q)
{
    FILE *logfile = fopen("../../legendre_otc.log", "w");
    fprintf(logfile, "opened log file.\n");

    // Initialize index and coeficient variables
    int mm_index, neg_mm_index, mp1m_index, neg_mp1m_index;
    int nm_index, nm1m_index, nm2m_index, neg_nm_index;
    int pp_index, neg_pp_index, pm_index, pm1m_index, neg_pm_index, neg_pm1m_index;
    int np1m_index, np2m_index, n, m;
    double mm_coef, Qpm1m_coef, neg_pm_coef, neg_pm1m_coef, neg_mm_coef, neg_nm_coef, neg_pp_coef, Qpp1p_coef;
    gsl::vector H;
    int N = u.size();
    int sp = (p + 1) * (p + 1); // Number of P and Q functions to compute

    // Check that all u-values are valid: u>1
    for (int j = 0; j < N; ++j)
    {
        if (u(j) <= 1.)
        {
            cerr << "all u-values must be strictly greater than 1.";
            break;
        }
    }

    P.resize(sp, N);
    Q.resize(sp, N);

    // Vectors of n and m values corresponding to each index.
    vector<int> nn(sp);
    vector<int> mm(sp);
    for (int i = 0; i < sp; ++i)
    {
        nn[i] = sqrt(floor(i));
        mm[i] = i - nn[i] * (nn[i] + 1);
    }

    // P_0^0(u)=1
    for (int j = 0; j < N; ++j)
    {
        P.set(0, j, 1.);
    }
    fprintf(logfile, "P_0^0 computed.\n");

    // populate with starting values, P_m^m (m=1:p) and P_{m+1}^m (m=1:p-1)
    if (p > 0)
    {
        for (int m = 0; m <= p; ++m)
        {
            if (m > 0)
            {
                mm_index = geti(m, m);
                mm_coef = gsl_sf_doublefact(2. * m - 1.);
                neg_mm_index = geti(m, -m);
                neg_mm_coef = 1. / gsl_sf_fact(2. * m);

                for (int j = 0; j < N; ++j)
                {
                    P.set(mm_index, j, mm_coef * pow(u(j) * u(j) - 1., 0.5 * m)); // P_m^m  
                    P.set(neg_mm_index, j, neg_mm_coef * P.get(mm_index, j)); // P_m^{-m}
                }

                fprintf(logfile, "P_{%d}^{%d} computed.\n",m,m );
                fprintf(logfile, "P_{%d}^{%d} computed.\n",m,-m );
            }

            if (m < p)
            {
                mm_index = geti(m, m);
                mp1m_index = geti(m + 1, m);
                neg_mp1m_index = geti(m + 1, -m);

                for (int j = 0; j < N; ++j)
                {
                    P.set(mp1m_index, j, (2 * m + 1) * u(j) * P.get(mm_index, j)); // P_{m+1}^m 
                    P.set(neg_mp1m_index, j, 1. / gsl_sf_fact(2. * m + 1.) * P.get(mp1m_index, j)); // P_{m+1}^{-m}
                }

                fprintf(logfile, "P_{%d}^{%d} computed.\n", m+1, m );
                fprintf(logfile, "P_{%d}^{%d} computed.\n", m+1, -m );
            }
        }
    }

    // Use recursion if higher order is needed.
    if (p > 1)
    {
        for (int m = 0; m <= p - 2; ++m)
        {
            // Use forward recursion to compute more P_n^m.
            for (int n = m + 2; n <= p; ++n)
            {
                nm_index = geti(n, m);
                nm1m_index = geti(n - 1, m);
                nm2m_index = geti(n - 2, m);
                neg_nm_index = geti(n, -m);
                neg_nm_coef = (gsl_sf_fact(1. * (n - m))) / gsl_sf_fact(1. * (n + m));

                for (int j = 0; j < N; ++j)
                {
                    P.set(nm_index, j, ((2 * n - 1) * u(j) * P.get(nm1m_index, j) - (n + m - 1) * P.get(nm2m_index, j)) / (n - m));
                    P.set(neg_nm_index, j, neg_nm_coef * P.get(nm_index, j));
                }
                fprintf(logfile, "P_{%d}^{%d} computed.\n",n,m );
                fprintf(logfile, "P_{%d}^{%d} computed.\n",n,-m );
            }
        }
    }

    fprintf(logfile, "P complete. \n" );

    //--------------------------------------------------------------------------------------
    // now calculate Qs:
    // Calculate highest order separately, using continued fraction and Wronskian
    H = cont_frac(p+1,p,u);
    pp_index=geti(p,p);
    neg_pp_index=geti(p,-p);
    neg_pp_coef=1./gsl_sf_fact(2.*p);
    Qpp1p_coef=gsl_sf_fact(2.*p)*pow(-1,p);

    for (int j=0; j<N; ++j)
    {
        Q.set(pp_index,j, Qpp1p_coef/ (((2*p+1)*u(j)-H(j))*P.get(pp_index,j)) );     //Q_p^p
        Q.set(neg_pp_index,j, neg_pp_coef * Q.get(pp_index,j));                      //Q_p^{-p}
    }

    fprintf(logfile , "Q_{%d}^{%d} computed.\n",p,p);
    fprintf(logfile , "Q_{%d}^{%d} computed.\n",p,-p);

    if (p > 0)
    {

        for (int m=0; m<=p-1; ++m)
        {
            H=cont_frac(p,m,u);
            pm_index=geti(p,m);
            pm1m_index=geti(p-1,m);
            neg_pm_index=geti(p,-m);
            neg_pm1m_index=geti(p-1,-m);

            Qpm1m_coef=(1.*gsl_sf_fact(p+m-1.))/gsl_sf_fact(p-m+0.)*pow(-1,m);
            neg_pm_coef=(1.*gsl_sf_fact(p-m+0.))/gsl_sf_fact(p+m+0.);
            neg_pm1m_coef=(1.*gsl_sf_fact(p-1.-m))/gsl_sf_fact(p-1.+m);    

            for (int j=0; j<N; ++j)
            {
                Q.set(pm1m_index,j, Qpm1m_coef / (P.get(pm_index,j)-H(j)*P.get(pm1m_index,j)));
                Q.set(pm_index,j, H(j)*Q.get(pm1m_index,j));

                if (m>0)
                {
                    Q.set(neg_pm1m_index,j, neg_pm1m_coef*Q.get(pm1m_index,j));
                    Q.set(neg_pm_index,j, neg_pm_coef*Q.get(pm_index,j));
                }

            }

            fprintf(logfile, "Q_{%d}^{%d} computed.\n", p,m);
            fprintf(logfile, "Q_{%d}^{%d} computed.\n", p-1,m);
            fprintf(logfile, "Q_{%d}^{%d} computed.\n", p,-m);
            fprintf(logfile, "Q_{%d}^{%d} computed.\n", p-1,-m);

            // Use recursion to compute remaining Qnm's
            if (p > 1)
            {
                //Use backward recursion to compute Q_n^m for n=p-2:m
                for (int n=p-2; n>=m; --n)
                {
                    nm_index=geti(n,m);
                    np1m_index=geti(n+1,m);
                    np2m_index=geti(n+2,m);
                    neg_nm_index=geti(n,-m);
                    neg_nm_coef=(1.*gsl_sf_fact(n-m+0.))/gsl_sf_fact(n+m+0.);

                    for (int j=0; j<N; ++j)
                    {
                        Q.set(nm_index, j, ((2*n+3)*u(j)*Q.get(np1m_index,j)-(n-m+2)*Q.get(np2m_index,j))/(n+m+1.));

                        if (m>0)
                        {
                            Q.set(neg_nm_index, j, neg_nm_coef * Q.get(nm_index,j));
                        }
                    }

                    fprintf(logfile, "Q_{%d}^{%d} computed.\n", n,m);
                    fprintf(logfile, "Q_{%d}^{%d} computed.\n", n,-m);
                
                }
            }
        }
    }

    fprintf(logfile, "Q complete.\n");
    fclose(logfile);
}




/*!
 * \brief Computes associated legendre functions off the cut of the first kind, $$P_n^m(u)$$ and their derivatives.
 * \param p The order of the spheroidal harmonics
 * \param u The u-values at which to compute the associated legendre functions
 * \param P Matrix of associated legendre functions off the cut of the first kind. Each row corresponds to a different n and m value, and each column corresponds to a different u-value.
 * \param dP Matrix of derivatives of associated legendre functions off the cut of the first kind. Each row corresponds to a different n and m value, and each column corresponds to a different u-value.
 * \note Will cause an error if any u-value is less than or equal to 1.
 */
void Dlegendre_otc(int p, gsl::vector u, gsl::matrix &P, gsl::matrix &dP)
{
    FILE* logfile= fopen("../../Dlegendre_otc.log","w");

    // Initialize index and coefficient variables
    int nm_index, np1m_index, np2m_index, n, m;
    int N=u.size();
    int sp=(p+1)*(p+1); //Number of dP functions to compute

    gsl::matrix P_pp1;
    legendre_otc(p+1, u, P_pp1); // Need to compute P up to order p+1 in order to get dP to order p

    dP.resize(sp, N);
    P.resize(sp, N);

    // Vectors of n and m values corresponding to each index.
    vector<int> nn(sp);
    vector<int> mm(sp);
    for (int i=0; i<sp; ++i)
    {
        nn[i]=sqrt(floor(i));
        mm[i]=i-nn[i]*(nn[i]+1);
    }

    // Compute derivatives using recurrence relation dPnm= [(m-n-1)*P_{n+1}^m(u) + (n+1)*u*P_n^m(u)]/(1-u^2)
    // --------------------------------------------------------------------------------
    fprintf(logfile, "computing dP:\n");

    for (int k=0; k<sp; ++k)
    {
        n=nn[k];
        m=mm[k];
        nm_index=geti(n,m);
        np1m_index=geti(n+1,m);

        for (int j=0; j<N; ++j)
        {
            dP.set(nm_index, j, ((m-n-1)*P_pp1.get(np1m_index,j)+(n+1)*u(j)*P_pp1.get(nm_index,j))/(1.-u(j)*u(j)) );
            P.set(nm_index, j, P_pp1.get(nm_index, j));
        }

        fprintf(logfile, "dP_{%d}^{%d} computed.\n",n,m);
    }
    fprintf(logfile,  "dP complete.\n");
    fclose(logfile);
}


/*!
 * \brief Computes associated legendre functions off the cut of the first and second kind, $$P_n^m(u)$$ and $$Q_n^m(u)$$, and their derivatives, for $$u>1$$.
 * \param p The order of the spheroidal harmonics
 * \param u The u-values at which to compute the associated legendre functions
 * \param P Matrix of associated legendre functions off the cut of the first kind. Each row corresponds to a different n and m value, and each column corresponds to a different u-value.
 * \param Q Matrix of associated legendre functions off the cut of the second kind. Each row corresponds to a different n and m value, and each column corresponds to a different u-value.
 * \param dP Matrix of derivatives of associated legendre functions off the cut of the first kind. Each row corresponds to a different n and m value, and each column corresponds to a different u-value.
 * \param dQ Matrix of derivatives of associated legendre functions off the cut of the second kind. Each row corresponds to a different n and m value, and each column corresponds to a different u-value.
 * \note Will cause an error if any u-value is less than or equal to 1.
 */
void Dlegendre_otc(int p, gsl::vector u, gsl::matrix &P, gsl::matrix &Q, gsl::matrix &dP, gsl::matrix &dQ)
{
    FILE* logfile= fopen("../../Dlegendre_otc.log","w");

    // Initialize index and coefficient variables
    int nm_index, np1m_index, np2m_index, n, m;
    int N=u.size();
    int sp=(p+1)*(p+1); //Number of dP, dQ functions to compute

    gsl::matrix P_pp1, Q_pp1;
    legendre_otc(p+1, u, P_pp1, Q_pp1); // Need to compute P up to order p+1 in order to get dP to order p

    dP.resize(sp, N);
    P.resize(sp, N);
    dQ.resize(sp, N);
    Q.resize(sp, N);

    // Vectors of n and m values corresponding to each index.
    vector<int> nn(sp);
    vector<int> mm(sp);
    for (int i=0; i<sp; ++i)
    {
        nn[i]=sqrt(floor(i));
        mm[i]=i-nn[i]*(nn[i]+1);
    }

    // Compute derivatives using recurrence relation dPnm= [(m-n-1)*P_{n+1}^m(u) + (n+1)*u*P_n^m(u)]/(1-u^2)
    // --------------------------------------------------------------------------------
    fprintf(logfile, "computing dP, dQ:\n");

    for (int k=0; k<sp; ++k)
    {
        n=nn[k];
        m=mm[k];
        nm_index=geti(n,m);
        np1m_index=geti(n+1,m);

        for (int j=0; j<N; ++j)
        {
            dP.set(nm_index, j, ((m-n-1)*P_pp1.get(np1m_index,j)+(n+1)*u(j)*P_pp1.get(nm_index,j))/(1.-u(j)*u(j)) );
            P.set(nm_index, j, P_pp1.get(nm_index, j));
            dQ.set(nm_index, j, ((m-n-1)*Q_pp1.get(np1m_index,j)+(n+1)*u(j)*Q_pp1.get(nm_index,j))/(1.-u(j)*u(j)) );
            Q.set(nm_index, j, Q_pp1.get(nm_index, j));
        }

        fprintf(logfile, "dP_{%d}^{%d} computed.\n",n,m);
        fprintf(logfile, "dQ_{%d}^{%d} computed.\n",n,m);
    }
    fprintf(logfile,  "dP, dQ complete.\n");
    fclose(logfile);
}
