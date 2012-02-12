#!/usr/bin/env python
import sys

import numpy as np

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la
import pickle as pickle
from copy import *
from potential2 import *

import time

from axis import *

n=75
order=4
bctype='x'

lb=0.
ub=10.

y=Axis(bctype,n,lb,ub,'fem',order)#+Axis('x',n,lb,ub,'fem',order)#+Axis('r',2,2.,3.,'fem')

npart=np.floor(n/(order-1))

if bctype=='cos':
        nelem=npart*(order-1)+1
elif bctype=='r' or bctype=='rho':
        nelem=npart*(order-1)
else:
        nelem=npart*(order-1)-1

#Overlap Matrix
overlap=y.overlap()
overlap_inv=y.overlap_inv()

#d|d matrix
B=np.zeros([nelem,nelem])
#Potential Matrix
V1=np.zeros([nelem,nelem])
V2=np.zeros([nelem,nelem])

#Control Parameter
Lambda=1

iter1=0

#Coefficient Matrix
U=np.zeros(nelem)

for e in y.e:
        b=e.matrix('d|d')
        v=e.matrix('fwell')
	v2=e.matrix('pot2')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
                        V1[iter1+k1,iter1+k2]=V1[iter1+k1,iter1+k2]+v[k1][k2]
                        V2[iter1+k1,iter1+k2]=V2[iter1+k1,iter1+k2]+v2[k1][k2]
        iter1=iter1+iter2-1

[evals,evecs]=la.eig(B/2+V1,y.overlap())

[mom_evals, mom_evecs]=la.eig(B/2,overlap)

#Normalization
for l in range(0,int(nelem)):
	norm=np.sqrt(y.FEM_InnerProduct(evecs.T[l],evecs.T[l]))
	evecs.T[l]=evecs.T[l] / norm

#Potential Modification
for k in range(0,int(nelem)):
        tempmat=np.outer(evecs.T[k],evecs.T[k])
        V2=V2+np.dot(tempmat,overlap)

#Scattering Energy
nenergy=20
Etot=np.linspace(4.5,5,nenergy)+2j

#Momentum Eigenstates - Not sure if this works
momentum_eigenstates=np.zeros([nenergy,nelem])+0j
for k in range(0,nenergy):
	momentum_eigenstates[k]=y.FEM_MomentumEigenstate(np.sqrt(Etot[k]))

niter=40
eps=2j

for k in range(0,nenergy):
	print abs(y.FEM_InnerProduct(momentum_eigenstates[0],momentum_eigenstates[k]))


for l in range(0,nenergy):
	E=Etot[l]
	#G0=la.inv(E*y.overlap()-B/2-V2+eps)
	G0_orig=la.inv(E*y.overlap()-B/2)
	#T=V1-V2
	#tempmat=np.dot(V1-V2,G0)

	#Gex=la.inv(E*y.overlap()-B/2-V1+eps)
	#tempmat2=np.dot(Gex,V1)
	#Tex=V1+np.dot(V1,np.dot(Gex,V1))

	#for k in range(0,niter):
	#	T=V1-V2+np.dot(tempmat,T)
	#	t[l][k]=la.norm(T)

	#print A
	#print Aex

	#print a

	#print la.norm(a)
