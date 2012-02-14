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

n=100
order=6
bctype='x'

lb=0.
ub=10.

y=Axis(bctype,n,lb,ub,'fem',order)

ax_cos=Axis('cos',10,0.,1.,'fem',order)
for k in range(1,10):
        ax_cos=ax_cos+Axis('cos',10,0.,1.,'fem',order)

npart=np.floor(n/(order-1))

if bctype=='cos':
        nelem=npart*(order-1)+1
elif bctype=='r' or bctype=='rho':
        nelem=npart*(order-1)
else:
        nelem=npart*(order-1)-1

nelem_cos=nelem+2

#Overlap Matrix
overlap=y.overlap()
overlap_inv=y.overlap_inv()

#d|d matrix
B=np.zeros([nelem,nelem])
B_cos=np.zeros([nelem_cos,nelem_cos])
#Potential Matrix
V1=np.zeros([nelem,nelem])
V2=np.zeros([nelem,nelem])

#Control Parameter
Lambda=1

iter1=0

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

iter1=0
for e in ax_cos.e:
        b=e.matrix('d|d')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B_cos[iter1+k1,iter1+k2]=B_cos[iter1+k1,iter1+k2]+b[k1][k2]
        iter1=iter1+iter2-1

[evals,evecs]=la.eig(B/2+V1,y.overlap())

#[mom_evals, mom_evecs]=la.eig(B/2,overlap)

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
Etot=np.linspace(0,19,nenergy)

#Momentum Eigenstates - Not sure if this works
momentum_eigenstates=np.zeros([nenergy,nelem_cos])+0j
for k in range(0,nenergy):
	momentum_eigenstates[k]=ax_cos.FEM_MomentumEigenstate(np.sqrt(Etot[k]))

print ax_cos.FEM_InnerProduct(momentum_eigenstates[0],momentum_eigenstates[1])

#print Etot
niter=40
eps=2j
