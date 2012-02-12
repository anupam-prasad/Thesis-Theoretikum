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

n=4
order=2
bctype='cos'

lb=-1.
ub=1.

y=Axis(bctype,n,lb,ub,'fem',order)

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

momentum=0

def free_wavefun(Momentum,x):
	a=np.exp(-Momentum*x*1j)
	return a

#Convert momentum eigenstate into finite element representation
V=np.zeros(int(nelem))+0j
iter1=0
for e in y.e:
	(x,w) = e.quadrature(add=10)
	(temp_a,temp_b)=e.val(x[0])
	iter2=len(temp_a)
	quad_iter=len(x)
	v=np.zeros([quad_iter,iter2])+0j
	
	for k in range(0,quad_iter):
		(a,b) = e.val(x[k])+0j
		c=a * np.exp(-momentum*x[k]*1j)
		for l in range(0,iter2):
			v[k][l]=c[l]*w[k]
#		raw_input()
#		v[k]=c[0]*w[k]

	integral=v.sum(axis=0)
	for k in range(0,iter2):
		V[iter1+k]=V[iter1+k]+integral[k]
	iter1=iter1+iter2-1


#raw_input()
y.plot()
plt.show()	

print V
