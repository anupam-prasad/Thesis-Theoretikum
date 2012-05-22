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

n=200
order=41

lb=0.
ub=10.

y=Axis('xopen',n,lb,ub,'fem',order)
#d|d matrix
B=np.zeros([y.len(),y.len()])

#Potential Matrix
V1=np.zeros([y.len(),y.len()])
V2=np.zeros([y.len(),y.len()])

V0=.001

for n in range(len(y.e)):
	i0=y.e[n].i0
	i1=i0+y.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y.e[n].matrix('d|d')
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y.e[n].matrix('fwell', np.array([4.,6.,V0]))
	V2[i0:i1,i0:i1]=V2[i0:i1,i0:i1]+y.e[n].matrix('fwell', np.array([6.,7.,V0]))


[evals,evecs]=la.eig(B/2. + V1 + V2,y.overlap())


#Plane Wave States
n0=0
nenergy=30
mom_evecs=np.zeros([y.len(),nenergy])+0j
mom_evals=np.zeros(nenergy)+0j
for k in range(nenergy):
	mom_evecs[:,k]=y.FEM_function(np.exp,(n0+k)*myPi*1j/10)
	mom_evals[k]=.5*(k*myPi/10)*(k*myPi/10)
	y.FEM_plot(mom_evecs[:,k].real,mom_evecs[:,k].imag)

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

#Normalization and Potential Modification
for k in range(0,y.len()):
	evecsnorm=np.sqrt(y.FEM_InnerProduct(evecs[:,k],evecs[:,k]))
	evecs[:,k]=evecs[:,k] / evecsnorm



