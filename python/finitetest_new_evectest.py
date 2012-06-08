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
order=41

lb=-5.
ub=5.

y=Axis('xopen',n,lb,ub,'fem',order)

#d|d matrix
B=np.zeros([y.len(),y.len()])

#Potential Matrix
V1=np.zeros([y.len(),y.len()])
V2=np.zeros([y.len(),y.len()])

V0=1.

for n in range(len(y.e)):
	i0=y.e[n].i0
	i1=i0+y.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y.e[n].matrix('d|d')
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y.e[n].matrix('fwell', np.array([-1.,1.,V0]))
#	V2[i0:i1,i0:i1]=V2[i0:i1,i0:i1]+y.e[n].matrix('gaussian', np.array([6.,7.,V0]))


#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

#Normalization and Potential Modification
for k in range(0,y.len()):
	evecsnorm=np.sqrt(y.FEM_InnerProduct(evecs[:,k],evecs[:,k]))
	cos_evecsnorm=np.sqrt(y.FEM_InnerProduct(cos_evecs[:,k],cos_evecs[:,k]))

	cos_evecs[:,k]=np.sqrt(ub-lb) * cos_evecs[:,k] / cos_evecsnorm	
	if evals[k]<0:
		evecs[:,k]=evecs[:,k] / evecsnorm
	else:
		evecs[:,k]=np.sqrt(ub-lb) * evecs[:,k] / evecsnorm
	

eps=1e-7j
n_bound=sum(evals<0)+150

Vmod=np.dot(y.overlap_inv(),V1+V2)
for l in range(11):
	mom_evec=y.FEM_function(np.exp,np.sqrt(2*cos_evals[l])*1j)
#	G0_test=np.zeros([y.len(),y.len()])
#	for k in range(21):
#		G0_test=G0_test+y.FEM_Outer(cos_evecs[:,k],cos_evecs[:,k]) / (cos_evals[l]+eps-cos_evals[k])

	print abs(np.dot(mom_evec.conjugate(),np.dot(B/2.,mom_evec)) / (ub-lb))
