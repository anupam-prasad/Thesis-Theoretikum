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

[prob_cos_evals, prob_cos_evecs]=la.eig(Bcos/2.+Vcos,ycos.overlap())
[prob_sin_evals, prob_sin_evecs]=la.eig(Bsin/2.+Vsin,ysin.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(cos_evals)
cos_evals=cos_evals[perm]
cos_evecs=cos_evecs[:,perm]

perm=np.argsort(sin_evals)
sin_evals=np.append([0],sin_evals[perm])
sin_evecs=sin_evecs[:,perm]

perm=np.argsort(prob_cos_evals)
prob_cos_evals=prob_cos_evals[perm]
prob_cos_evecs=prob_cos_evecs[:,perm]

perm=np.argsort(prob_sin_evals)
prob_sin_evals=prob_sin_evals[perm]
prob_sin_evecs=prob_sin_evecs[:,perm]

nbound_cos=sum(prob_cos_evals<0)
nbound_sin=sum(prob_sin_evals<0)

temp_evecs=np.zeros([ysin.len(),ysin.len()+1])
temp_evecs[:,1:]=sin_evecs
sin_evecs=temp_evecs

Lambda=0
Gam_cos=np.zeros([ycos.len(),ycos.len()])
Gam_sin=np.zeros([ysin.len(),ysin.len()])

for k in range(nbound_cos):
	Gam_cos=Gam_cos+ycos.FEM_Outer(prob_cos_evecs[:,k],prob_cos_evecs[:,k])

for k in range(nbound_sin):
	Gam_sin=Gam_sin+ysin.FEM_Outer(prob_sin_evecs[:,k],prob_sin_evecs[:,k])

eps = 0j
niter=2
nenergy=3

Tmatrix_elements=np.zeros([niter+1,nenergy])+0j

V_right=np.dot(ycos.overlap_inv(),Vcos)
V_left=np.dot(Vcos,ycos.overlap_inv())

Vmod_right=np.dot(ycos.overlap_inv(),Vcos-Lambda*Gam_cos)
Vmod_left=np.dot(Vcos-Lambda*Gam_cos,ycos.overlap_inv())

Bmod = np.dot(Bcos/2.+Lambda*Gam_cos,ycos.overlap_inv())
Borig_mod = np.dot(Bcos/2.+Vcos,ycos.overlap_inv())
for k in range(nenergy):
	Emat = np.eye(ycos.len()) * (cos_evals[k] + eps)
	G0=np.dot(la.inv(Emat - Bmod),ycos.overlap())
	Gorig=np.dot(la.inv(Emat - Borig_mod),ycos.overlap())

	Torig=Vcos+np.dot(V_left,np.dot(Gorig,V_right))
	Tmatrix_elements[0,k] = np.dot(cos_evecs[:,k],np.dot(Torig,cos_evecs[:,k]))

	Tmatrix_elements[1,k] = np.dot(cos_evecs[:,k],np.dot(Vcos,cos_evecs[:,k]))
	Tmat_right=Vmod_right
	for l in range(2,niter+1):
		Tmat = Vcos + np.dot(Vmod_left,np.dot(G0,Tmat_right))
		Tmatrix_elements[l,k] = np.dot(cos_evecs[:,k],np.dot(Tmat,cos_evecs[:,k]))
		Tmat_right=np.dot(ycos.overlap_inv(),Tmat)

V_right=np.dot(ysin.overlap_inv(),Vsin)
V_left=np.dot(Vsin,ysin.overlap_inv())

Vmod_left=np.dot(ysin.overlap_inv(),Vsin-Lambda*Gam_sin)
Vmod_right=np.dot(Vsin-Lambda*Gam_sin,ysin.overlap_inv())

Bmod = np.dot(Bsin/2.+Lambda*Gam_sin,ysin.overlap_inv())
Borig_mod = np.dot(Bsin/2.+Vsin,ysin.overlap_inv())
for k in range(nenergy):
	Emat = np.eye(ysin.len()) * (sin_evals[k] + eps)
	G0=np.dot(la.inv(Emat - Bmod),ysin.overlap())
	Gorig=np.dot(la.inv(Emat - Borig_mod),ysin.overlap())

	Torig=Vsin+np.dot(V_left,np.dot(Gorig,V_right))
	Tmatrix_elements[0,k] += np.dot(sin_evecs[:,k],np.dot(Torig,sin_evecs[:,k]))

	Tmatrix_elements[1,k] += np.dot(sin_evecs[:,k],np.dot(Vsin,sin_evecs[:,k]))
	Tmat_right=Vmod_right
	for l in range(2,niter+1):
		Tmat = Vsin + np.dot(Vmod_left,np.dot(G0,Tmat_right))
		Tmatrix_elements[l,k] += np.dot(sin_evecs[:,k],np.dot(Tmat,sin_evecs[:,k]))
		Tmat_right=np.dot(ysin.overlap_inv(),Tmat)

print abs(Tmatrix_elements)
