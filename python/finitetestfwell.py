#!/usr/bin/env python
import sys

import numpy as np

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la
from copy import *

from axis import *

n=300
order=26
bctype='xopen'

lb=0.
ub=10.
	
y=Axis(bctype,n,lb,ub,'fem',order)

overlap=y.overlap()
overlap_inv=y.overlap_inv()

#d|d matrix
B=np.zeros([y.len(),y.len()])
#Potential Matrix
V=np.zeros([y.len(),y.len()])

iter1=0

for e in y.e:
	b=e.matrix('d|d')
	v=e.matrix('fwell', np.array([4.,6.,100]))
	iter2=int(np.sqrt(np.size(b)))
	for k1 in range(0,iter2):
		for k2 in range(0,iter2):
			B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
			V[iter1+k1,iter1+k2]=V[iter1+k1,iter1+k2]+v[k1][k2]
	iter1=iter1+iter2-1

[evals,evecs]=la.eig(B/2.,overlap)

sorted_indices=abs(evals).argsort()
evals_sorted=evals[sorted_indices]
evecs_sorted=evecs.T[sorted_indices]
print evals_sorted
#print sorted1/sorted1[1]
y.FEM_plot(evecs_sorted[3])

#for k in range(0,y.len()):
#print y.FEM_InnerProduct(evecs.T[0],evecs.T[k])
