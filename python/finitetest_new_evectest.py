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

V0=1

for n in range(len(y.e)):
	i0=y.e[n].i0
	i1=i0+y.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y.e[n].matrix('d|d')
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y.e[n].matrix('gaussian', np.array([1.9,2.1,V0]))
	V2[i0:i1,i0:i1]=V2[i0:i1,i0:i1]+y.e[n].matrix('gaussian', np.array([7.9,8.1,V0]))


[evals,evecs]=la.eig(B/2. + V1 + V2,y.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

#Normalization and Potential Modification
for k in range(0,y.len()):
	evecsnorm=np.sqrt(y.FEM_InnerProduct(evecs[:,k],evecs[:,k]))
	evecs[:,k]=evecs[:,k] / evecsnorm

n_bound=sum(evals<0)+3
eps=1e-10j

#Plane Wave States
B0_mod=np.dot(B/2.,y.overlap_inv())
B_orig=np.dot(B/2. + V1 + V2,y.overlap_inv())

Vmod_left=np.dot(V1+V2,y.overlap_inv())
Vmod_right=np.dot(y.overlap_inv(),V1+V2)

#mom_evals=.5*myPi*myPi/100
mom_evals=evals[n_bound]
mom_evecs=y.FEM_function(np.exp,np.sqrt(2*mom_evals)*1j)

#print abs(np.dot(mom_evecs,np.dot(B/2.,mom_evecs)))
Emat=y.overlap() * (mom_evals+eps)

#Free Green's Function
G0=la.inv(-B/2. + Emat)
#G0=np.dot(G0temp,y.overlap())

Gtemp_orig=la.inv(-B_orig + Emat)
G_orig=np.dot(Gtemp_orig,y.overlap())

#T_orig=V1+V2+np.dot(Vmod_left,np.dot(G_orig,Vmod_right))
tempvec=mom_evecs
for k in range(20):
	tempvec=mom_evecs+np.dot(np.dot(G0,V1+V2),tempvec)

tempvec=tempvec / np.sqrt(y.FEM_InnerProduct(tempvec,tempvec))

print np.dot(evecs[:,n_bound],np.dot(B/2.+V1+V2,evecs[:,n_bound]))
print np.dot(tempvec.conjugate(),np.dot(B/2.+V1+V2,tempvec))
#print y.FEM_InnerProduct(tempvec,tempvec)

y.FEM_plot(abs(mom_evecs))
y.FEM_plot(abs(evecs[:,n_bound]),abs(tempvec))
y.FEM_plot(tempvec.real,tempvec.imag)
