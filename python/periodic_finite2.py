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

n=600
order=5

lb=-50.
ub=50.

ycos=Axis('xopen',n,lb,ub,'fem',order)
ysin=Axis('x',n,lb,ub,'fem',order)

Bcos=np.zeros([ycos.len(),ycos.len()])
Bsin=np.zeros([ysin.len(),ysin.len()])

Vcos=np.zeros([ycos.len(),ycos.len()])
Vsin=np.zeros([ysin.len(),ysin.len()])

V0=1

for n in range(len(ycos.e)):
	i0=ycos.e[n].i0
	i1=i0+ycos.e[n].n
	Bcos[i0:i1,i0:i1]=Bcos[i0:i1,i0:i1]+ycos.e[n].matrix('d|d')
	Vcos[i0:i1,i0:i1]=Vcos[i0:i1,i0:i1]+ycos.e[n].matrix('fwell', np.array([-1.,1.,V0]))

[cos_evals,cos_evecs]=la.eig(Bcos/2.,ycos.overlap())

for n in range(len(ysin.e)):
	i0=ysin.e[n].i0
	i1=i0+ysin.e[n].n
	Bsin[i0:i1,i0:i1]=Bsin[i0:i1,i0:i1]+ysin.e[n].matrix('d|d')
	Vsin[i0:i1,i0:i1]=Vsin[i0:i1,i0:i1]+ysin.e[n].matrix('fwell', np.array([-1.,1.,V0]))	

[sin_evals,sin_evecs]=la.eig(Bsin/2.,ysin.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(cos_evals)
cos_evals=cos_evals[perm]
cos_evecs=cos_evecs[:,perm]

perm=np.argsort(sin_evals)
sin_evals=np.append([0],sin_evals[perm])
sin_evecs=sin_evecs[:,perm]

temp_evecs=np.zeros([ysin.len(),ysin.len()+1])
temp_evecs[:,1:]=sin_evecs
sin_evecs=temp_evecs

eps = 1j
niter=2
nenergy=5

Tmatrix_elements=np.zeros([niter,nenergy])

for k in range(5):
	Emat = np.eye(ycos.len()) * (cos_evals[k] + eps)
	Bmod = np.dot(Bcos/2.,ycos.overlap_inv())
	Tmat=Vcos
	Tmatrix_elements[0,k] = np.dot(cos_evecs[:,k],np.dot(Tmat,cos_evecs[:,k]))
	G0=np.dot(la.inv(Emat - Bmod),ycos.overlap())
	for l in range(niter-1,niter):
		Tmat = Vcos + np.dot(Vcos,np.dot(G0,Tmat))
		Tmatrix_elements[l,k] = np.dot(cos_evecs[:,k],np.dot(Tmat,cos_evecs[:,k]))

for k in range(5):
	Emat = np.eye(ysin.len()) * (sin_evals[k] + eps)
	Bmod = np.dot(Bsin/2.,ysin.overlap_inv())
	Tmat=Vsin
	Tmatrix_elements[0,k] = Tmatrix_elements[0,k] + sin_evecs[:,k],np.dot(Tmat,sin_evecs[:,k])
	G0=np.dot(la.inv(Emat - Bmod),ysin.overlap())
	for l in range(niter-1,niter):
		Tmat = Vsin + np.dot(Vsin,np.dot(G0,Tmat))
		Tmatrix_elements[l,k] = Tmatrix_elements[l,k] + np.dot(sin_evecs[:,k],np.dot(Tmat,sin_evecs[:,k]))

