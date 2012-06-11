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

V0=100

iter1=0

#t=np.linspace(0,10,1001)
#tpot=potential(t,'fwell',np.array([4.5,5.5,.001]))
#
#plt.plot(t,tpot)
#plt.show()
#
#exit('here')

for e in y.e:
        b=e.matrix('d|d')
        v=e.matrix('gaussian', np.array([6.,7.,V0]))
        v2=e.matrix('gaussian', np.array([3.,4.,V0]))
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
                        V1[iter1+k1,iter1+k2]=V1[iter1+k1,iter1+k2]+v[k1][k2]
                        V2[iter1+k1,iter1+k2]=V2[iter1+k1,iter1+k2]+v2[k1][k2]
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
Lambda=1000
#Potential Modification
for k in range(0,y.len()):
	cosnorm=np.sqrt(2 * y.FEM_InnerProduct(cos_evecs[:,k],cos_evecs[:,k]) / (ub-lb))
	cos_evecs[:,k]=cos_evecs[:,k] / cosnorm
	
	evecsnorm=np.sqrt(y.FEM_InnerProduct(evecs[:,k],evecs[:,k]))
	evecs[:,k]=evecs[:,k] / evecsnorm
	if evals[k]<0:
		Gam=Gam + y.FEM_Outer(evecs[:,k],evecs[:,k])
		print evals[k]


niter=30
eps=1e-8j

n0=60
nenergy=30
store1=np.zeros([niter+3,nenergy])+0j

#Free Green's Operator
Bmod=np.dot(B/2. + Lambda*Gam,y.overlap_inv())

Bmod1=np.dot(B/2. + V1 + V2 + Lambda*Gam,y.overlap_inv())
B_orig=np.dot(B/2. + V1 + V2,y.overlap_inv())

Vtemp=np.dot(V1+V2,y.overlap_inv())
Vmod=np.dot(y.overlap_inv(),Vtemp)

Gamtemp=np.dot(Lambda*Gam,y.overlap_inv())
Gam_ver=np.dot(y.overlap_inv(),Gamtemp)

for k in range(nenergy):
	Emat=np.eye(y.len()) * (cos_evals[k+n0]+eps)

	Gtemp=la.inv(-Bmod + Emat)
	Gtemp1=la.inv(-Bmod1 + Emat)
	Gtemp2=la.inv(-B_orig + Emat)

	G0=np.dot(Gtemp,y.overlap())
	Gexact=np.dot(Gtemp1,y.overlap())
	G_orig=np.dot(Gtemp2,y.overlap())

	vec1=np.dot(Gexact,cos_evecs[:,k+n0])
	store1[0,k]=np.dot(cos_evecs[:,k+n0],vec1)

	vec1=np.dot(G_orig,cos_evecs[:,k+n0])
	store1[1,k]=np.dot(cos_evecs[:,k+n0],vec1)

	reverse_trans=np.eye(y.len(),y.len())
	for l in range(y.len()):
		if evals[l] < 0:
			A=np.dot(y.FEM_Outer(evecs[:,l],evecs[:,l]),y.overlap_inv())
			reverse_trans=reverse_trans-Lambda*A / (cos_evals[k+n0]-evals[l]+eps)
		else: break
	
	Gmat=G0
	vec1=np.dot(Gmat,cos_evecs[:,k+n0])
	store1[3,k]=np.dot(cos_evecs[:,k+n0],vec1)
	VG=np.dot(G0,Vmod)

	for l in range(1,niter):
		Gstore=Gmat
		Gmat=G0+np.dot(VG,Gmat)
		Gver=np.dot(reverse_trans,Gmat)
		vec1=np.dot(Gver,cos_evecs[:,k+n0])
		store1[l+3,k]=np.dot(cos_evecs[:,k+n0],vec1)		

	Gver=np.dot(reverse_trans,Gmat)
	vec1=np.dot(Gver,cos_evecs[:,k+n0])
	store1[2,k]=np.dot(cos_evecs[:,k+n0],vec1)

for k in range(nenergy):
	print abs(store1[:,k]), cos_evals[k+n0]
	raw_input()

#if Lambda==0:
#	fname='test5plot/test5results_gaussian_unmod'
#else:
#	fname='test5plot/test5results_gaussian_mod'
#
#f=open(fname,'w')
#pickle.dump(store1,f)

