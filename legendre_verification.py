#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import math

P = np.loadtxt('data/P.txt')
Q = np.loadtxt('data/Q.txt')
dP = np.loadtxt('data/dP.txt')
dQ = np.loadtxt('data/dQ.txt')


p=int(sys.argv[1])
ua=float(sys.argv[2])
ub=float(sys.argv[3])
N=int(sys.argv[4])

[sp,Ndata]=np.shape(P)

if Ndata != N:
    raise Exception("Ndata =/= N. Check that input arguments agree with legendre_otc results.")

if (p+1)**2 != sp:
    raise Exception("(p+1)^2 =/= sp. Check that input arguments agree with legendre_otc results.")

u=np.linspace(ua,ub,N)

print("length(u)=",len(u))
print("N = ",N)

geti= lambda n,m : m+n**2+n  # different from Matlab, 0-based

n=3; m=1; # for this recurrence, n-m must be >= 2
# Recurrence relation: (2n-1)x Q_{n-1}^m = (n-m)Q_n^m +(n+m-1)Q_{n-2}^m
EQ=abs((n-m)*Q[geti(n,m),:]-(2*n-1)* u * Q[geti(n-1,m),:]+(n+m-1)*Q[geti(n-2,m),:])
#  Wronskian Relation: P_n^m Q_{n-1}^m - P_{n-1}^m Q_n^m =(-1)^m (n+m-1)!/(n-m)!
Ew=abs(P[geti(n,m),:]*Q[geti(n-1,m),:]-P[geti(n-1,m),:]*Q[geti(n,m),:]-math.factorial(n+m-1)/math.factorial(n-m)*(-1)**m)



plt.figure()
plt.semilogy(u,EQ)
plt.xlabel('u')
plt.ylabel('log10 abs. Error')
plt.title('Q forward recurrence relation (n='+str(n)+', m='+str(m)+')')
plt.savefig('figures/EQ'+str(n)+str(m)+'.png')

plt.figure()
plt.semilogy(u,Ew)
plt.xlabel('u')
plt.ylabel('log10 abs. Error')
plt.title('Wronskian relation (n='+str(n)+', m='+str(m)+')')
plt.savefig('figures/Ew'+str(n)+str(m)+'.png')


# Test Derivatives with Wronskian Relationship
Wrsnk=np.zeros(np.shape(P))
for k in range(sp):
    n = math.floor(np.sqrt(k))
    m=k-n**2-n
    # fprintf('%d, %d\n',n,m)
    # print(geti(n,m))
    Wrsnk[k,:]=abs(P[k,:]*dQ[k,:]-dP[k,:]*Q[k,:]-math.factorial(n+m)*(-1)**m /(math.factorial(n-m)*(1-u**2)) )/abs(math.factorial(n+m)*(-1)**m /(math.factorial(n-m)*(1-u**2)))


plt.figure()
plt.imshow(np.log10(Wrsnk))
plt.colorbar()
plt.xlabel('u values')
plt.ylabel('(n,m) index')
plt.title('log10 rel. error of Wronskian Relationship with derivatives')
plt.savefig('figures/Ew2.png')