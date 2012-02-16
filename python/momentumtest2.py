#!/usr/bin/env python
import sys

import numpy as np

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la

from axis import *

n=2
order=2
bctype='cos'

lb=0.
ub=1.

ncos=2

ax_cos=Axis(bctype,n,lb,ub,'fem',order)

for k in range(0,ncos-1):
	ax_cos=ax_cos+Axis(bctype,n,lb,ub,'fem',order)

npart=ncos*np.floor(n/(order-1))

if bctype=='cos':
        nelem=npart*(order-1)+1
elif bctype=='r' or bctype=='rho':
        nelem=npart*(order-1)
else:
        nelem=npart*(order-1)-1

#d|d matrix
B_cos=np.zeros([nelem,nelem])
B_test=np.zeros([nelem-2,nelem-2])

iter1=0
for e in ax_cos.e:
        b=e.matrix('d|d')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B_cos[iter1+k1,iter1+k2]=B_cos[iter1+k1,iter1+k2]+b[k1][k2]
        iter1=iter1+iter2-1

#Scattering Energy
Etot=np.linspace(0,nelem-1,nelem)*2*myPi+0j

#Momentum Eigenstates - Not sure if this works
momentum_eigenstates=np.zeros([nelem,nelem])+0j
for k in range(0,int(nelem)):
	momentum_eigenstates[k]=ax_cos.FEM_MomentumEigenstate(Etot[k])

[mom_evals,mom_evecs]=la.eig(B_cos/2,ax_cos.overlap())

print mom_evals
print mom_evecs.T

#print momentum_eigenstates[0]
