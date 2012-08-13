from __future__ import division
from numpy import *
import scipy
import sys

# This subroutine computes the propagator from a Markovian transition matrix R, for time t
def compute_T(R, t):
	l, u  = linalg.eig(R)
	l, v  = linalg.eig(R.transpose())
	vt=v.transpose()
	lmb = diag(l)
	R_reconstr=(u.dot(diag(l))).dot(vt)
	b1 = allclose(R,R_reconstr)
	b2 = allclose(vt.dot(u),identity(R.shape[0]))
	if(b1==False or b2==False):
		print "An error occurred while trying to exponentiate R, quitting..."
		sys.exit(1)
	T=(u.dot(diag(exp(l*t)))).dot(vt)
	return T.real

def propagate(p0, R, t):
	T = compute_T(R,t)
	pt = (T.dot(p0)).real
	if (allclose(sum(pt),1.0) == False):
		print "The propagated distribution is not normalized, quitting..."
		sys.exit(1)
	return pt

# R = array([[-1,1/2,1/2], [1/2,-1, 1/2], [1/2,1/2, -1]])

R = array([[-1, 0, 1], 
           [1, -1, 0], 
           [0,  1, -1]])

t=1.23 # seconds
T=compute_T(R,t)

alpha=sqrt(3)/2
A=cos(alpha*t)
B=cos(alpha*t-2*pi/3)
C=cos(alpha*t-4*pi/3)

xx1 = 1/3*array([[1, 1, 1], 
                 [1, 1, 1], 
                 [1, 1, 1]])
xx2 = 2/3*exp(-3/2*t)*array([[A, C, B], 
                             [B, A, C], 
                             [C, B, A]])
T_analyt = xx1+xx2
print allclose(T,T_analyt)
