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

V0=0

for n in range(len(y.e)):
	i0=y.e[n].i0
	i1=i0+y.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y.e[n].matrix('d|d')
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y.e[n].matrix('fwell', np.array([3.,4.,V0]))
	V2[i0:i1,i0:i1]=V2[i0:i1,i0:i1]+y.e[n].matrix('fwell', np.array([6.,7.,V0]))


[evals,evecs]=la.eig(B/2. + V1 + V2,y.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

#Normalization and Potential Modification
Lambda=0
Gam=np.zeros([y.len(),y.len()])+0j
for k in range(0,y.len()):
	evecsnorm=np.sqrt(y.FEM_InnerProduct(evecs[:,k],evecs[:,k]))
	evecs[:,k]=evecs[:,k] / evecsnorm
	if evals[k]<0:
		Gam=Gam+y.FEM_Outer(evecs[:,k],evecs[:,k])

n0=0
n_bound=sum(evals<0)
nenergy=10
niter=20
eps=-1e-6j

store1=np.zeros([niter+3,nenergy])+0j

#Plane Wave States
mom_evecs=np.zeros([y.len(),nenergy])+0j
mom_evals=np.zeros(nenergy)+0j

B0=np.dot(B/2.+Lambda*Gam,y.overlap_inv())+0j
B_orig=np.dot(B/2. + V1 + V2,y.overlap_inv())+0j
B_exact=np.dot(B/2. + V1 + V2 + Lambda*Gam,y.overlap_inv())+0j

Vmod_left=np.dot(V1+V2,y.overlap_inv())+0j
Vmod_right=np.dot(y.overlap_inv(),V1+V2)+0j

for k in range(nenergy):
	mom_evecs[:,k]=y.FEM_function(np.exp,(n0+k)*myPi*1j/10)
	mom_evals[k]=np.dot(mom_evecs[:,k].conjugate(),np.dot(B/2.,mom_evecs[:,k])) / y.FEM_InnerProduct(mom_evecs[:,k],mom_evecs[:,k])

	Emat=np.eye(y.len()) * (mom_evals[k]+eps)

        #Free Green's Function
	G0temp=la.inv(-B0 + Emat)
        G0=np.dot(G0temp,y.overlap())

        Gtemp_orig=la.inv(-B_orig + Emat)
        G_orig=np.dot(Gtemp_orig,y.overlap())

        Gtemp_exact=la.inv(-B_exact + Emat)
        Gexact=np.dot(Gtemp_exact,y.overlap())

        vec1=np.dot(Gexact,mom_evecs[:,k])
        store1[0,k]=np.dot(mom_evecs[:,k].conjugate(),vec1)

        vec1=np.dot(G_orig,mom_evecs[:,k])
        store1[1,k]=np.dot(mom_evecs[:,k].conjugate(),vec1)

 	Gmat=G0
        vec1=np.dot(Gmat,mom_evecs[:,k])
        store1[3,k]=np.dot(mom_evecs[:,k].conjugate(),vec1)
        VG=np.dot(G0,Vmod_right)

        for l in range(1,niter):
                Gmat=G0+np.dot(VG,Gmat)
                vec1=np.dot(Gmat,mom_evecs[:,k])
                store1[l+3,k]=np.dot(mom_evecs[:,k].conjugate(),vec1)
		Gmat=np.dot(y.overlap_inv(),Gmat)

for k in range(nenergy):
        print abs(store1[:,k]), mom_evals[k]
        raw_input()
