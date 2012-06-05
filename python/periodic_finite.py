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

n=2
order=2

lb=0.
ub=1.

y=Axis('xopen',n,lb,ub,'fem',order)
#d|d matrix
B=np.zeros([y.len(),y.len()])

#Potential Matrix
V1=np.zeros([y.len(),y.len()])
V2=np.zeros([y.len(),y.len()])

V0=10.

for n in range(len(y.e)):
	i0=y.e[n].i0
	i1=i0+y.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y.e[n].matrix('d|d')
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y.e[n].matrix('fwell', np.array([0.,1.,V0]))
#	V2[i0:i1,i0:i1]=V2[i0:i1,i0:i1]+y.e[n].matrix('gaussian', np.array([6.,7.,V0]))


#print B
#print y.overlap()

[evals,evecs]=la.eig(B/2. + V1 + V2,y.overlap())
[cos_evals,cos_evecs]=la.eig(B/2.,y.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

perm=np.argsort(cos_evals)
cos_evals=cos_evals[perm]
cos_evecs=cos_evecs[:,perm]

print cos_evecs
print cos_evals
exit('here')
#Normalization and Potential Modification
for k in range(0,y.len()):
	evecsnorm=np.sqrt(y.FEM_InnerProduct(evecs[:,k],evecs[:,k]))
	cos_evecsnorm=np.sqrt(y.FEM_InnerProduct(cos_evecs[:,k],cos_evecs[:,k]))

	cos_evecs[:,k]=np.sqrt(ub-lb) * cos_evecs[:,k] / cos_evecsnorm	
	if evals[k]<0:
		evecs[:,k]=evecs[:,k] / evecsnorm
	else:
		evecs[:,k]=np.sqrt(ub-lb) * evecs[:,k] / evecsnorm
	y.FEM_plot(cos_evecs[:,k])
