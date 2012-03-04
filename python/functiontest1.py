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

#d|d matrix
B=np.zeros([y.len(),y.len()])

iter1=0

for e in y.e:
        b=e.matrix('d|d')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
        iter1=iter1+iter2-1

Ptot=myPi/2+0j
momentum_eigenstates=np.zeros([2,y.len()])+0j
for k in range(4,6):
	momentum_eigenstates[k-4]=y.FEM_function(np.cos,Ptot*k)+y.FEM_function(np.sin,Ptot*k)*1j

#y.FEM_plot(momentum_eigenstates[1,:].real,momentum_eigenstates[1,:].imag)

print y.FEM_InnerProduct(momentum_eigenstates[0],momentum_eigenstates[0])

val1=np.dot(B/2,momentum_eigenstates[0])

n1=abs(np.dot(momentum_eigenstates[0].conjugate(),val1))


val1=np.dot(B/2,momentum_eigenstates[1])

n2=abs(np.dot(momentum_eigenstates[1].conjugate(),val1))

print n2 / n1
