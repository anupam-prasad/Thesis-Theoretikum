#!/usr/bin/env python
import sys

import numpy as np

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la
from copy import *

from axis import *

n=200
order=6
bctype='x'

lb=-5.
ub=5.
	
y=Axis(bctype,n,lb,ub,'fem',order)

npart=np.floor(n/(order-1))

if bctype=='cos':
	nelem=npart*(order-1)+1
elif bctype=='r' or bctype=='rho':
	nelem=npart*(order-1)
else:
	nelem=npart*(order-1)-1

overlap=y.overlap()
overlap_inv=y.overlap_inv()

#d|d matrix
B=np.zeros([nelem,nelem])
#Potential Matrix
V=np.zeros([nelem,nelem])

iter1=0

for e in y.e:
	b=e.matrix('d|d')
	v=e.matrix('qho')
	iter2=int(np.sqrt(np.size(b)))
	for k1 in range(0,iter2):
		for k2 in range(0,iter2):
			B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
			V[iter1+k1,iter1+k2]=V[iter1+k1,iter1+k2]+v[k1][k2]
	iter1=iter1+iter2-1

[evals,evecs]=la.eig(B/2.+V,overlap)

print evals[181]
y.FEM_plot(evecs.T[181])
