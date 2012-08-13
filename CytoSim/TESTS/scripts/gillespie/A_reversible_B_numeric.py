# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=2>

# Time-Dependent Probability Distribution for the A<->B Reaction

# <markdowncell>

# We will use $p(t)=e^{R*t}p(0)$\, to compute the time evolution of the copy numbers of A

# <codecell>

from scipy import *
from scipy import linalg

# <codecell>

N=10 # The number of states
t=0.48 # time elapsed in seconds
kf=2.5 # forward rate in s^-1
kb=2.5 # backward rate in s^-1

# <codecell>

R=zeros([N+1,N+1])

# <codecell>

for n in r_[1:N]: 
    R[n,n]=-(kf*n+kb*(N-n))
    R[n,n+1]=kf*(n+1)
    R[n,n-1]=kb*(N-(n-1))
R[0,0]=-kb*N
R[10,10]=-kf*N
R[0,1]=kf*1
R[N-1,N]=kf*N
R[N,N-1]=kb*1
print R
print sum(R,axis=0)

# <codecell>

ONES=ones(N+1)
print ONES

# <codecell>

dot(ONES,R)

# <codecell>

T=linalg.expm2(R*t)

# <codecell>

dot(ONES,T)

# <codecell>

p0=zeros(N+1)
p0[0]=1
print p0

# <codecell>

pt=dot(T,p0)
sum(pt)

# <codecell>

print pt

# <markdowncell>

# **Analytical answer from the file: A_reversible_B_explicit_solved.mw**

# <codecell>

p_analyt=array([0.002327151517,0.1940040583e-1, 0.7277956968e-1, .1617947293, .2360416130, .2361326132, .1640442006, 0.781464013e-1, 0.244301648e-1, 0.4525850e-2, 0.377291e-3])
print p_analyt

# <codecell>

allclose(p_analyt,pt)

# <codecell>

for i in r_[0:N+1]: print "P[", i, "] = ", pt[N-i]
