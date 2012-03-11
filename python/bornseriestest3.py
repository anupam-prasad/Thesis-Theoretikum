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
ub=10.

y=Axis('xopen',n,lb,ub,'fem',order)
ysin=Axis('x',n,lb,ub,'fem',order)
#d|d matrix
B=np.zeros([y.len(),y.len()])
Bsin=np.zeros([ysin.len(),ysin.len()])
#Potential Matrix
V1=np.zeros([y.len(),y.len()])
V2=np.zeros([y.len(),y.len()])

iter1=0

for e in y.e:
        b=e.matrix('d|d')
        v=e.matrix('fwell', np.array([4.,6.,0.]))
	v2=e.matrix('pot2')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
                        V1[iter1+k1,iter1+k2]=V1[iter1+k1,iter1+k2]+v[k1][k2]
                        V2[iter1+k1,iter1+k2]=V2[iter1+k1,iter1+k2]+v2[k1][k2]
        iter1=iter1+iter2-1

iter1=0

for e in ysin.e:
        b=e.matrix('d|d')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        Bsin[iter1+k1,iter1+k2]=Bsin[iter1+k1,iter1+k2]+b[k1][k2]
        iter1=iter1+iter2-1

[evals,evecs]=la.eig(B/2+V1,y.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

#Momentum Eigenstates
[cos_evals,cos_evecs]=la.eig(B/2.,y.overlap())
[sin_evals_temp,sin_evecs_temp]=la.eig(Bsin/2.,ysin.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(cos_evals)
cos_evals=cos_evals[perm]
cos_evecs=cos_evecs[:,perm]

cos_evals=cos_evals[0:ysin.len()+1]
cos_evecs=cos_evecs[:,0:ysin.len()+1]

perm=np.argsort(sin_evals_temp)
sin_evals_temp=sin_evals_temp[perm]
sin_evecs_temp=sin_evecs_temp[:,perm]

sin_evals=np.zeros([ysin.len()+1])+0j
sin_evecs=np.zeros([y.len(),ysin.len()+1])+0j

for k in range(0,ysin.len()):
	sin_evals[k+1]=sin_evals_temp[k]
	sin_evecs[1:ysin.len()+1,k+1]=sin_evecs_temp[:,k]

#Normalization
mom_evals=np.zeros([ysin.len()+1])+0j
for k in range(0,ysin.len()+1):
	cosnorm=np.sqrt(2. * y.FEM_InnerProduct(cos_evecs[:,k],cos_evecs[:,k]))
	sinnorm=np.sqrt(2. * y.FEM_InnerProduct(sin_evecs[:,k],sin_evecs[:,k]))

	cos_evecs[:,k]=cos_evecs[:,k] / cosnorm
	if k>0:
		sin_evecs[:,k]=sin_evecs[:,k] / sinnorm
	

mom_evecs=cos_evecs-sin_evecs*1j

for k in range(0,ysin.len()+1):
	vec1=np.dot(B/2.,mom_evecs[:,k])
	mom_evals[k]=abs(np.dot(mom_evecs[:,k].conjugate(),vec1))

#Potential Modification
Lambda=1
for k in range(0,int(y.len())):
        V2=V2+Lambda*y.FEM_Outer(evecs.T[k],evecs.T[k])

niter=5
eps=1j

#Free Green's Operator
#Bmod=np.dot(B/2.,y.overlap_inv())
#Vmod=np.dot(V1,y.overlap_inv())
#epsmat=np.eye(y.len())*eps
store1=np.zeros([2,40])+0j
Tmat=V1

for k in range(40):
#	Emat=np.eye(y.len())*(mom_evals[k]+eps)
	Emat=y.overlap()*(mom_evals[k]+eps)
	Gtemp=la.inv(-B/2. + Emat)
	G0=np.dot(Gtemp,y.overlap())
#	tempmat=np.dot(V1,G0)

#	tempmat2 = la.inv(np.eye(y.len())-tempmat)
#	Texact = np.dot(tempmat2,V1)
#	Gexact = la.inv(-Bmod-Vmod+Emat)

#	vec1=np.dot(G0,evecs.T[k])
#	store1[0,k]=np.dot(evecs.T[k].conjugate(),vec1) / 10.
#	store1[1,k]=np.dot(evecs.T[k],vec1) / 10.

	vec1=np.dot(G0,mom_evecs[:,k].real)
#	store1[0,k]=np.dot(mom_evecs[:,k].conjugate(),vec1)
	store1[0,k]=y.FEM_InnerProduct(mom_evecs[:,k].real,vec1)
	store1[1,k]=np.dot(mom_evecs[:,k],vec1)

#	for l  in range(0,niter):
#		Tmat=Tmat+np.dot(tempmat,Tmat)	

#	vec1=np.dot(Tmat,momentum_eigenstates[k])
#	store1[0,k]=np.dot(momentum_eigenstates[0].conjugate(),vec1)
#	store1[1,k]=np.dot(momentum_eigenstates[0],vec1)

print abs(store1[0,:])

#print abs(store1[1,:])

#print abs(store1[0,:]) / abs(store1[1,:])
