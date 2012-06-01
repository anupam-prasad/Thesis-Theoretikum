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

lb=-25.
ub=25.

y=Axis('xopen',n,lb,ub,'fem',order)
#d|d matrix
B=np.zeros([y.len(),y.len()])

#Potential Matrix
V1=np.zeros([y.len(),y.len()])
V2=np.zeros([y.len(),y.len()])

V0=-1.

for n in range(len(y.e)):
	i0=y.e[n].i0
	i1=i0+y.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y.e[n].matrix('d|d')
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y.e[n].matrix('fwell', np.array([-10.,10.,V0]))
#	V2[i0:i1,i0:i1]=V2[i0:i1,i0:i1]+y.e[n].matrix('gaussian', np.array([6.,7.,V0]))


[evals,evecs]=la.eig(B/2. + V1 + V2,y.overlap())
#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

#Normalization and Potential Modification
for k in range(0,y.len()):
	evecsnorm=np.sqrt(y.FEM_InnerProduct(evecs[:,k],evecs[:,k]))
	if evals[k]<0:
		evecs[:,k]=evecs[:,k] / evecsnorm
	else:
		evecs[:,k]=np.sqrt(ub-lb) * evecs[:,k] / evecsnorm	
#	if evals[k]>=-V0:
	print 2*myPi / np.sqrt(2*evals[k])
	print 2*myPi / np.sqrt(2*(evals[k]+V0))
	y.FEM_plot(evecs[:,k])

exit('here')
n_bound=sum(evals<0)+150
eps=1e-7j

B0_mod=np.dot(y.overlap_inv(),B/2.)
B_orig=np.dot(y.overlap_inv(),B/2. + V1 + V2)
Vmod=np.dot(y.overlap_inv(),V1+V2)

#Plane Wave States
#mom_evals=myPi*myPi*64/200
mom_evals=evals[n_bound]
mom_evecs=y.FEM_function(np.exp,np.sqrt(2*mom_evals)*1j)

Emat=np.eye(y.len()) * (mom_evals+eps)

#Free Green's Function
G0=la.inv(-B0_mod + Emat)
G_orig=la.inv(-B_orig + Emat)

tempvec=mom_evecs
for k in range(20):
	tempvec=mom_evecs+np.dot(np.dot(G0,Vmod),tempvec)

tempvec = tempvec / np.sqrt(y.FEM_InnerProduct(tempvec,tempvec) / (ub-lb))

tempvec2=mom_evecs+np.dot(np.dot(G_orig,Vmod),mom_evecs)
tempvec2 = tempvec2 / np.sqrt(y.FEM_InnerProduct(tempvec2,tempvec2) / (ub-lb))

print abs(np.dot(tempvec.conjugate(),np.dot(B/2.+V1+V2,tempvec))) / 50.
print abs(np.dot(tempvec2.conjugate(),np.dot(B/2.+V1+V2,tempvec2))) / 50.
y.FEM_absplot(tempvec2,tempvec)
