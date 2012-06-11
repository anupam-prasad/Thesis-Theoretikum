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

n=75
order=4
bctype='x'

lb=0.
ub=10.

y=Axis(bctype,n,lb,ub,'fem',order)#+Axis('x',n,lb,ub,'fem',order)#+Axis('r',2,2.,3.,'fem')

free_ax=Axis('cos',n,-1.,1.,'fem',order)

npart=np.floor(n/(order-1))

if bctype=='cos':
        nelem=npart*(order-1)+1
elif bctype=='r' or bctype=='rho':
        nelem=npart*(order-1)
else:
        nelem=npart*(order-1)-1

#Overlap Matrix
overlap=y.overlap()
overlap_inv=y.overlap_inv()

#d|d matrix
B=np.zeros([nelem,nelem])

iter1=0

for e in y.e:
        b=e.matrix('d|d')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
        iter1=iter1+iter2-1

#Scattering Energy
nenergy=20
Ptot=np.linspace(0,19,nenergy)*myPi

def FEM_MomentumEigenstate(momentum,n,order,lb,ub):

	ax=Axis('cos',n,0.,1.,'fem',order)
	nelem=ax.n

	#Convert momentum eigenstate into finite element representation
	fem_vec=np.zeros(int(nelem))+0j
	invert_vec=np.zeros(int(nelem))+0j
	iter1=0

	for e in ax.e:
		(x,w) = e.quadrature(add=10)
		x=(ub-lb)*x+lb
		w=w*2/(ub-lb)
		(temp_a,temp_b)=e.val(x[0])
		iter2=len(temp_a)
		quad_iter=len(x)
		v=np.zeros([quad_iter,iter2])+0j

		for k in range(0,quad_iter):
			(a,b) = e.val(x[k])+0j
			c=a * np.exp(-momentum*x[k]*1j)
			for l in range(0,iter2):
				v[k][l]=c[l]*w[k]

		integral=v.sum(axis=0)
		for k in range(0,iter2):
			invert_vec[iter1+k]=invert_vec[iter1+k]+integral[k]
		iter1=iter1+iter2-1

	fem_vec=np.dot(invert_vec,ax.overlap_inv())
	return fem_vec


#Momentum Eigenstates - Not sure if this works
#momentum_eigenstates=np.zeros([nenergy,nelem])+0j
momentum_eigenstates_corrected=np.zeros([nenergy,int(nelem)])+0j
for k in range(0,nenergy):
	momentum_eigenstates=FEM_MomentumEigenstate(Ptot[k],n,order,lb,ub)
#	for l in range(0,int(nelem)-2):
	momentum_eigenstates_corrected[k]=momentum_eigenstates[1:nelem+1]

print y.FEM_InnerProduct(momentum_eigenstates_corrected[0], momentum_eigenstates_corrected[1])
#
#for k in range(0,nenergy):
#	vec=np.dot(B/2,momentum_eigenstates_corrected[k]) / np.dot(overlap,momentum_eigenstates_corrected[k])
#	print abs(vec[1:int(nelem)-1])
#	raw_input()
