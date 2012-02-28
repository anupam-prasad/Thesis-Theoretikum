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

n=400
order=26

lb=0.
ub=10.

y=Axis('xopen',n,lb,ub,'fem',order)

#Overlap Matrix
overlap=y.overlap()
overlap_inv=y.overlap_inv()

#d|d matrix
B=np.zeros([y.len(),y.len()])
#Potential Matrix
V1=np.zeros([y.len(),y.len()])
V2=np.zeros([y.len(),y.len()])

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

[evals,evecs]=la.eig(B/2+V1,y.overlap())

#[mom_evals, mom_evecs]=la.eig(B/2,overlap)

#Normalization
for l in range(0,int(y.len())):
	norm=np.sqrt(y.FEM_InnerProduct(evecs.T[l],evecs.T[l]))
	evecs.T[l]=evecs.T[l] / norm

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

#Potential Modification
Lambda=1
for k in range(0,int(y.len())):
        V2=V2+Lambda*y.FEM_outer(evecs.T[k],evecs.T[k])

#Scattering Energy
nenergy=40
Ptot=np.linspace(0,19,nenergy)*myPi+0j

#Momentum Eigenstates - Not sure if this works - It works.. kind of.
momentum_eigenstates=np.zeros([nenergy,y.len()])+0j
for k in range(0,nenergy):
	momentum_eigenstates[k]=y.FEM_function(np.exp,Ptot[k]*1j)

#print Etot
niter=40
eps=2j
