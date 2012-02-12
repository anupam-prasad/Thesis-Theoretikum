#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from axis import *

# input
n=200
L=100.
kind='fem'
base='legendre'
order=20

# harmonic osillcator
print 'harmonic oscillator'
ax=Axis('r',n,0,L,kind,order=order)
ham=np.zeros((ax.n,ax.n))
ovr=np.zeros((ax.n,ax.n))
for e in ax.e:
    ham[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('1/qdq|1/qdq')+e.matrix('coulomb(1.)')
    ovr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('|')

(val,vec)=la.eig(ham,ovr)      # solve eigenproblem
print np.sort(val.real)[:10]   # show results (real part of eigenvalues, sorted !)

ax.plot()
plt.show()

exit()

# hydrogen atom (L=0)
print 'hydrogen atom'
ax=Axis('r',n,0.,L,kind,order=order)
ham=np.zeros((ax.n,ax.n))
ovr=np.zeros((ax.n,ax.n))
for e in ax.e:
    ham[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('1/qdq|1/qdq')+e.matrix('coulomb(1.)')
    ovr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('|')


(val,vec)=la.eig(ham,ovr)      # solve eigenproblem
print np.sort(val.real)[:10]   # show results (real part of eigenvalues, sorted !)
