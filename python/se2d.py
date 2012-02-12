#! /usr/bin/env python
import sys

from numpy import *
import matplotlib.pyplot as plt

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la

from axis import *
from tensor import *
from math import *

from matplotlib import colors, ticker

pot='q^2/2'
neval=40

ax=[Axis('x',26,-5.,5.,kind='fem',order=8),Axis('y',26,-5.,5.,kind='fem',order=8)]

ham=[]
ovr=[]
for i in range(2):
    ham.append(np.zeros((ax[i].n,ax[i].n)))
    ovr.append(np.zeros((ax[i].n,ax[i].n)))
    for e in ax[i].e:
        ham[i][e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('d|d')+0.5*e.matrix('|q^2|')
        ovr[i][e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('|')

Ham=Tensor(TensorProd(ham[0],ovr[1]))+Tensor(TensorProd(ovr[0],ham[1]))
Ovr=Tensor(TensorProd(ovr[0],ovr[1]))
print 'eigenvalues'
print np.sort(la.eigvals(Ham.array(),Ovr.array()).real)[:neval]


arho=Axis('rho',50,0.,10.,kind='fem',order=5)
aphi=Axis('phi',8,0.,2.*myPi,kind='trigon')

hrho=np.zeros((arho.n,arho.n))
qrho=np.zeros((arho.n,arho.n))
srho=np.zeros((arho.n,arho.n))
for e in arho.e:
    hrho[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('d|d')+0.5*e.matrix('|q^2|')
    qrho[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('|1/q^2|')    
    srho[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('|')

dphi=np.zeros((aphi.n,aphi.n))
sphi=np.zeros((aphi.n,aphi.n))
for e in aphi.e:
    dphi[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('d|d')
    sphi[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('|')

Ham=Tensor(TensorProd(sphi,hrho))+Tensor(TensorProd(dphi,qrho))
Ovr=Tensor(TensorProd(sphi,srho))

plt.matshow(np.log10(Ham.array()),vmin=-4,vmax=1)
plt.show()

print 'eigenvalues'
print np.sort(la.eigvals(Ham.array(),Ovr.array()).real)[:neval]

