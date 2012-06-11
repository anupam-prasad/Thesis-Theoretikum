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
order=6
bctype='x'

lb=-5.
ub=5.

y=Axis(bctype,n,lb,ub,'fem',order)#+Axis('x',n,lb,ub,'fem',order)#+Axis('r',2,2.,3.,'fem')

#y.plot()
#plt.show()
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
#Potential Matrix
V=np.zeros([nelem,nelem])

iter1=0

#Coefficient Matrix
U=np.zeros(nelem)

Pin=y.FEM_MomentumEigenstate(5)
Pin=y.FEM_Normalize(Pin)

Pout=y.FEM_MomentumEigenstate(100)
Pout=y.FEM_Normalize(Pout)

for e in y.e:
#       for b in e.b:
        b=e.matrix('d|d')
        v=e.matrix('qho')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
                        V[iter1+k1,iter1+k2]=V[iter1+k1,iter1+k2]+v[k1][k2]
        iter1=iter1+iter2-1

[evals,evecs]=la.eig(B/2+V,y.overlap())

npoints=int(nelem*100)
t=np.linspace(lb,ub,npoints)
wavefun=np.zeros([nelem,npoints])
iter3=0
lenstore=0

