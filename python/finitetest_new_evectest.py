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

lb=0.
ub=100.

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
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y.e[n].matrix('gaussian', np.array([49.,50.,V0]))
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

n_bound=sum(evals<0)
eps=1e-10j

B0_mod=np.dot(y.overlap_inv(),B/2.)
B_orig=np.dot(y.overlap_inv(),B/2. + V1 + V2)

#Plane Wave States
mom_evals=myPi*myPi*64/200
#mom_evals=evals[n_bound]
mom_evecs=y.FEM_function(np.exp,np.sqrt(2*mom_evals)*1j)

Emat=np.eye(y.len()) * (mom_evals+eps)

#Free Green's Function
G0=la.inv(-B0_mod + Emat)

G_orig=la.inv(-B_orig + Emat)
#G_orig=np.dot(Gtemp_orig,y.overlap())

tempvec=mom_evecs
Gam=np.zeros([y.len(),y.len()])
Gam0=np.zeros([y.len(),y.len()])
#for k in range(n_bound):
#	Gam=Gam+y.FEM_Outer(evecs[:,k],evecs[:,k]) / (mom_evals+eps-evals[k])
#	Gam0=Gam0+y.FEM_Outer(evecs[:,k],evecs[:,k]) / y.FEM_InnerProduct(evecs[:,k],np.dot(G0,evecs[:,k]))
#	tempvec=tempvec-y.FEM_InnerProduct(evecs[:,k],mom_evecs) * evecs[:,k]

Gam=np.dot(y.overlap_inv(),np.dot(Gam,y.overlap_inv()))
Gam0=np.dot(y.overlap_inv(),np.dot(Gam0,y.overlap_inv()))

G_orig_mod=G_orig-Gam	
G0_mod=G0-np.dot(G0,np.dot(Gam0,G0))

for k in range(10):
	tempvec=mom_evecs+np.dot(np.dot(G0_mod,V1+V2),tempvec)

#tempvec=mom_evecs+np.dot(np.dot(G_orig,V1+V2),mom_evecs)
#tempvec=np.sqrt(ub-lb) * tempvec / np.sqrt(y.FEM_InnerProduct(tempvec,tempvec))

tempvec2=tempvec+np.dot(np.dot(G_orig_mod,V1+V2),tempvec)
tempvec2=tempvec2 / np.sqrt(y.FEM_InnerProduct(tempvec2,tempvec2)/(ub-lb))

#tempvec=tempvec+np.dot(np.dot(Gam,V1+V2),mom_evecs)
tempvec=tempvec / np.sqrt(y.FEM_InnerProduct(tempvec,tempvec)/(ub-lb))

y.FEM_absplot(tempvec,tempvec2)
exit('here')

tempvec2=np.dot((np.eye(y.len())+np.dot(G_orig_mod,V1+V2)),mom_evecs)
tempvec2=tempvec2 / np.sqrt(y.FEM_InnerProduct(tempvec2,tempvec2))
#print np.dot(evecs[:,n_bound],np.dot(B/2.+V1+V2,evecs[:,n_bound]))
print evals[:n_bound+1]
print np.dot(tempvec.conjugate(),np.dot(B/2.+V1+V2,tempvec))
print np.dot(mom_evecs.conjugate(),np.dot(B/2.,mom_evecs)) / 10
#print y.FEM_InnerProduct(tempvec,tempvec)

y.FEM_absplot(tempvec,evecs[:,n_bound])
