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

V0=100

for n in range(len(y.e)):
	i0=y.e[n].i0
	i1=i0+y.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y.e[n].matrix('d|d')
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y.e[n].matrix('gaussian', np.array([3.,4.,V0]))
	V2[i0:i1,i0:i1]=V2[i0:i1,i0:i1]+y.e[n].matrix('gaussian', np.array([6.,7.,V0]))


[evals,evecs]=la.eig(B/2. + V1 + V2,y.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

#Normalization and Potential Modification
Lambda=1e8
Gam_init=np.zeros([y.len(),y.len()])
for k in range(0,y.len()):
	evecsnorm=np.sqrt(y.FEM_InnerProduct(evecs[:,k],evecs[:,k]))
	evecs[:,k]=evecs[:,k] / evecsnorm
	if evals[k]<0:
		Gam_init=Gam_init+y.FEM_Outer(evecs[:,k],evecs[:,k])


n0=0
n_bound=sum(evals<0)
nenergy=10
niter=10
eps=1ej

store1=np.zeros([niter+3,nenergy])+0j

#Plane Wave States
mom_evecs=np.zeros([y.len(),nenergy])+0j
mom_evals=np.zeros(nenergy)+0j

B0_mod=np.dot(B/2.,y.overlap_inv())

B_orig=np.dot(B/2. + V1 + V2,y.overlap_inv())

Vmod_left=np.dot(V1+V2,y.overlap_inv())
Vmod_right=np.dot(y.overlap_inv(),V1+V2)

for k in range(nenergy):
	mom_evecs[:,k]=y.FEM_function(np.exp,(n0+k)*myPi*1j/10)
	mom_evals[k]=.5*((n0+k)*myPi/10)*((n0+k)*myPi/10)

	Emat=np.eye(y.len()) * (mom_evals[k]+eps)

        #Free Green's Function
	G0temp=la.inv(-B0_mod + Emat)
        G0=np.dot(G0temp,y.overlap())

	Gam=np.zeros([y.len(),y.len()])
	mat1=np.zeros([y.len(),y.len()])
        for l in range(n_bound):
		Gam=Gam+y.FEM_Outer(evecs[:,l],evecs[:,l])/(mom_evals[k]+eps-evals[l])
		mat1=mat1+y.FEM_Outer(evecs[:,l],evecs[:,l]) / np.dot(evecs[:,l],np.dot(G0,evecs[:,l]))
	 
	mat1=np.dot(y.overlap_inv(),np.dot(mat1,y.overlap_inv()))
 	G0_mod=G0-np.dot(G0,np.dot(mat1,G0))

        Gtemp_orig=la.inv(-B_orig + Emat)
        G_orig=np.dot(Gtemp_orig,y.overlap())

        Gexact=G_orig-Gam

        Texact=V1+V2+np.dot(Vmod_left,np.dot(Gexact,Vmod_right))
        T_orig=V1+V2+np.dot(Vmod_left,np.dot(G_orig,Vmod_right))

        vec1=np.dot(Texact,mom_evecs[:,k])
        store1[0,k]=np.dot(mom_evecs[:,k].conjugate(),vec1)

        vec1=np.dot(T_orig,mom_evecs[:,k])
        store1[1,k]=np.dot(mom_evecs[:,k].conjugate(),vec1)

        Tmat=Vmod_right
        vec1=np.dot(V1+V2,mom_evecs[:,k])
        store1[3,k]=np.dot(mom_evecs[:,k].conjugate(),vec1)
        VG=np.dot(Vmod_left,G0_mod)

        for l in range(1,niter):
                Tmat=V1+V2+np.dot(VG,Tmat)
                vec1=np.dot(Tmat,mom_evecs[:,k])
                store1[l+3,k]=np.dot(mom_evecs[:,k].conjugate(),vec1)
		Tmat=np.dot(y.overlap_inv(),Tmat)

for k in range(nenergy):
        print abs(store1[:,k]), mom_evals[k]
        raw_input()
