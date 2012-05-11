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

Gam=np.zeros([y.len(),y.len()])

V0=20

iter1=0

t=np.linspace(0,10,1001)
tpot=potential(t,'fwell',np.array([4.,6.,V0]))

plt.plot(t,tpot)
plt.show()


for e in y.e:
        b=e.matrix('d|d')
        v=e.matrix('fwell', np.array([4.,6.,V0]))
#        v2=e.matrix('gaussian', np.array([3.,4,V0]))
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
                        V1[iter1+k1,iter1+k2]=V1[iter1+k1,iter1+k2]+v[k1][k2]
#                        V2[iter1+k1,iter1+k2]=V2[iter1+k1,iter1+k2]+v2[k1][k2]
        iter1=iter1+iter2-1

[evals,evecs]=la.eig(B/2. + V1 + V2,y.overlap())

#Momentum Eigenstates
[cos_evals,cos_evecs]=la.eig(B/2.,y.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(cos_evals)
cos_evals=cos_evals[perm]
cos_evecs=cos_evecs[:,perm]

perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

#Normalization and Potential Modification
#Lambda=100.
for k in range(0,y.len()):
	cosnorm=np.sqrt(2 * y.FEM_InnerProduct(cos_evecs[:,k],cos_evecs[:,k]) / (ub-lb))
	cos_evecs[:,k]=cos_evecs[:,k] / cosnorm
	
	evecsnorm=np.sqrt(y.FEM_InnerProduct(evecs[:,k],evecs[:,k]))
	evecs[:,k]=evecs[:,k] / evecsnorm

	#Potential Modification
#	if evals[k]<0:
#	        Gam=Gam + y.FEM_Outer(evecs[:,k],evecs[:,k])

niter=15
eps=1e-1j

#Free Green's Operator
n0=0
nenergy=20
store1=np.zeros([niter+2,nenergy])+0j

Bmod=np.dot(B/2.,y.overlap_inv())
B_orig=np.dot(B/2. + V1 + V2,y.overlap_inv())

Vtemp=np.dot(V1+V2,y.overlap_inv())
Vmod=np.dot(y.overlap_inv(),Vtemp)

Tmat=Vmod


for k in range(nenergy):
	for l in range(0,y.len()):
		if evals[l]<0:	Gam=Gam+y.FEM_Outer(evecs[:,l],evecs[:,l])/(cos_evals[k+n0]-evals[l])
		else:	break

	Emat=np.eye(y.len()) * (cos_evals[k+n0]+eps)

	Gtemp=la.inv(-Bmod + Emat)
	Gtemp2=la.inv(-B_orig + Emat)

	G0=np.dot(Gtemp,y.overlap())
	G_orig=np.dot(Gtemp2,y.overlap())
	Gexact=G_orig-Gam

	vec1=np.dot(Gexact,cos_evecs[:,k+n0])
	store1[0,k]=np.dot(cos_evecs[:,k+n0],vec1)

	vec1=np.dot(G_orig,cos_evecs[:,k+n0])
	store1[1,k]=np.dot(cos_evecs[:,k+n0],vec1)

	Gmat=G0-Gam
	VG=np.dot(Gmat,Vmod)
	
	for l in range(niter):
		Gmat=G0+np.dot(VG,Gmat)
		vec1=np.dot(Gmat,cos_evecs[:,k+n0])
		store1[l+2,k]=np.dot(cos_evecs[:,k+n0],vec1)

#print abs(store1[0,:] * Lambda)
#

for k in range(nenergy):
	print abs(store1[:,k])

#print Lambda, niter

#f=open('test5plot/test5results2','w')
#pickle.dump(store1,f)

