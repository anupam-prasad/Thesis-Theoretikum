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

n=800
order=5

lb=0.
ub=50.

y=Axis('r',n,lb,ub,'fem',order)

B=np.zeros([y.len(),y.len()])

V=np.zeros([y.len(),y.len()])

V0=1

for n in range(len(y.e)):
	i0=y.e[n].i0
	i1=i0+y.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y.e[n].matrix('d|d')
	V[i0:i1,i0:i1]=V[i0:i1,i0:i1]+y.e[n].matrix('yukawa', np.array([0.,1.,V0]))

[mom_evals,mom_evecs]=la.eig(B/2.,y.overlap())

[evals,evecs]=la.eig(B/2.+V,y.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(mom_evals)
mom_evals=np.append([0],mom_evals[perm])
mom_evecs=mom_evecs[:,perm]

temp_evecs=np.zeros([y.len(),y.len()+1])
unit_vec=y.FEM_ones()
temp_evecs[:,0]=unit_vec
temp_evecs[:,1:]=mom_evecs
mom_evecs=temp_evecs

perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

nbound=sum(evals<0)

eps = 0j
niter=2
nenergy=3

Tmatrix_elements=np.zeros([niter+1,nenergy])+0j

V_right=np.dot(y.overlap_inv(),V)
V_left=np.dot(V,y.overlap_inv())

Vmod_right=np.dot(y.overlap_inv(),V)
Vmod_left=np.dot(V,y.overlap_inv())

Bmod = np.dot(B/2.,y.overlap_inv())
Borig_mod = np.dot(B/2.+V,y.overlap_inv())

for k in range(nenergy):		
	Emat = np.eye(y.len()) * (mom_evals[k] + eps)
	G0=np.dot(la.inv(Emat - Bmod),y.overlap())
	Gorig=np.dot(la.inv(Emat - Borig_mod),y.overlap())

	Torig=V+np.dot(V_left,np.dot(Gorig,V_right))
	Tmatrix_elements[0,k] = np.dot(unit_vec,np.dot(Torig,mom_evecs[:,k]))

	Tmatrix_elements[1,k] = np.dot(unit_vec,np.dot(V,mom_evecs[:,k]))
	Tmat_right=Vmod_right
	for l in range(2,niter+1):
		Tmat = V + np.dot(Vmod_left,np.dot(G0,Tmat_right))
		Tmatrix_elements[l,k] = np.dot(unit_vec,np.dot(Tmat,mom_evecs[:,k]))
		Tmat_right=np.dot(y.overlap_inv(),Tmat)

print abs(Tmatrix_elements)
