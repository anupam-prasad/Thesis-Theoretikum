#!/usr/bin/env python

#Trying to construct momentum eigenstate by appending many cosine elements 
#together
import sys

import numpy as np
#import my_pyplot
#import matplotlib as plt

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la
from copy import *
from potential2 import *

#from math import *
#import cmath

#from array import array

from axis import *

n=10
order=6
bctype='cos'

lb=0.
ub=1.

y=Axis(bctype,n,lb,ub,'fem',order)
for k in range(1,10):	
	y=y+Axis(bctype,n,lb,ub,'fem',order)

ax2=Axis('x',100,0.,10.,'fem',order)

npart=10*np.floor(n/(order-1))

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

#Coefficient Matrix
U=np.zeros(nelem)

for e in y.e:
#	for b in e.b:
	b=e.matrix('d|d')
	v=e.matrix('pot')
	iter2=int(np.sqrt(np.size(b)))
	for k1 in range(0,iter2):
		for k2 in range(0,iter2):
			B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
	iter1=iter1+iter2-1

#y.plot()
#plt.show()

#raw_input()

#ax2.plot()
#plt.show()
