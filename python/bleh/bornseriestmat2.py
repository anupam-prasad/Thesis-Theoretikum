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

n=300
order=26
bctype='xopen'

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
#	v2=e.matrix('pot2')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
                        V1[iter1+k1,iter1+k2]=V1[iter1+k1,iter1+k2]+v[k1][k2]
#                       V2[iter1+k1,iter1+k2]=V2[iter1+k1,iter1+k2]+v2[k1][k2]
        iter1=iter1+iter2-1

[evals,evecs]=la.eig(B/2+V1,y.overlap())

[mom_evals, mom_evecs]=la.eig(B/2,y.overlap())

#Scattering Energy
nenergy=20
Ptot=np.linspace(0,19,nenergy)*myPi+0j

momentum_eigenstates=np.zeros([nenergy,nelem])+0j
for k in range(0,nenergy):
	momentum_eigenstates[k]=y.FEM_function(np.exp,Ptot[k]*1j)


niter=40
t=np.zeros([nenergy,nenergy])+0j

eps=2j

for k in range(0,nenergy):
	for l in range(0,nenergy):
		tempvec=y.FEM_InnerProduct(momentum_eigenstates[k],V1)
		t[k][l]=y.FEM_InnerProduct(tempvec,momentum_eigenstates[l])
