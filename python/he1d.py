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

read_input_open(sys.argv[1])

neval=read_input(3,'.number of eigenvalues',doc='how many eigenvalues to print')
nadd=read_input(5,'.extra quadrature points',doc='how many points to add for quadrature')


ham=[]
ovr=[]
ax=[]
for i in range(2):
    ax.append(axis_read(i))
    ham.append(np.zeros((ax[i].n,ax[i].n)))
    ovr.append(np.zeros((ax[i].n,ax[i].n)))
    for e in ax[i].e:
        ham[i][e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('d|d')+e.matrix('cou1d(2.,0.5)')
        ovr[i][e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('|')

read_input_finish()

print 'single particle',np.sort(la.eigvals(ham[0].real,ovr[0].real)).real[:neval]


Ham=Tensor(TensorProd(ham[0],ovr[1]))+Tensor(TensorProd(ovr[0],ham[1]))
Ovr=Tensor(TensorProd(ovr[0],ovr[0]))

Vee=np.zeros((ax[0].n,ax[1].n,ax[0].n,ax[1].n))
for e in ax[0].e:
    print 'interaction',e.i0
    for f in ax[1].e:
        Vee[e.i0:e.i0+e.n,f.i0:f.i0+f.n,e.i0:e.i0+e.n,f.i0:f.i0+f.n]+=interaction(e,f,'cou1d(1.,0.5)',nadd)

Hmat=Ham.array().real-np.reshape(Vee,(ax[0].n*ax[1].n,ax[0].n*ax[1].n)).real

ev,vec=la.eig(Hmat,Ovr.array().real)
ls=np.argsort(ev.real)
print 'two particles',ev[ls[:neval]].real

#def plotwf(ax,vec,box,pts=100):
#    # convert to plot grid
#    x0=np.linspace(box[0,0],box[1,0],pts)
#    x1=np.linspace(box[0,1],box[1,1],pts)
#    ij=-1
#    for x in x0:
#        for y in x1:
            



#plt.matshow(np.log10(Ham.array()),vmin=-4,vmax=1)
#plt.show()
