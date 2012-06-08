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
order=11

lb=0.
ub=1.

y=Axis('x',n,lb,ub,'fem',order)+Axis('x',n,ub,2*ub,'fem',order)

B=np.zeros([y.len(),y.len()])

for n in range(len(y.e)):
	i0=y.e[n].i0
	i1=i0+y.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y.e[n].matrix('d|d')

[mom_evals,mom_evecs]=la.eig(B/2.,y.overlap())

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(mom_evals)
mom_evals=mom_evals[perm]
mom_evecs=mom_evecs[:,perm]

print mom_evals[:10]

y.FEM_plot(mom_evecs[:,2],mom_evecs[:,3])
