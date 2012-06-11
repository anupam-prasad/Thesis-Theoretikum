#!/usr/bin/env python
import sys

import numpy as np

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la

from axis import *

n=18
order=7
bctype='cos'

lb=0.
ub=1.

ncos=10

ax_cos=Axis(bctype,n,lb,ub,'fem',order)
test_axis=Axis('x',ncos*n,0.,ncos*ub,'fem',order)

for k in range(0,ncos-1):
	ax_cos=ax_cos+Axis(bctype,n,lb,ub,'fem',order)

npart=ncos*np.floor(n/(order-1))

if bctype=='cos':
        nelem=npart*(order-1)+1
elif bctype=='r' or bctype=='rho':
        nelem=npart*(order-1)
else:
        nelem=npart*(order-1)-1

#d|d matrix
B_cos=np.zeros([nelem,nelem])
B_test=np.zeros([nelem-2,nelem-2])

iter1=0
for e in ax_cos.e:
        b=e.matrix('d|d')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B_cos[iter1+k1,iter1+k2]=B_cos[iter1+k1,iter1+k2]+b[k1][k2]
        iter1=iter1+iter2-1

iter1=0
for e in test_axis.e:
        b=e.matrix('d|d')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B_test[iter1+k1,iter1+k2]=B_test[iter1+k1,iter1+k2]+b[k1][k2]
        iter1=iter1+iter2-1

#Scattering Energy
nenergy=20
Etot=np.linspace(0,nenergy-1,nenergy)*2*myPi+0j

#Momentum Eigenstates - Not sure if this works
momentum_eigenstates=np.zeros([nenergy,nelem])+0j
for k in range(0,nenergy):
	momentum_eigenstates[k]=ax_cos.FEM_MomentumEigenstate(Etot[k])

mom_eigenstates_test=momentum_eigenstates

#Might Need To Rewrite Inner Product Function
#for k in range(0,nenergy):
#	print np.abs(ax_cos.FEM_InnerProduct(momentum_eigenstates[k],momentum_eigenstates[k]))

[mom_evals,mom_evecs]=la.eig(B_cos/2,ax_cos.overlap())

[mom_evals_sin,mom_evecs_sin]=la.eig(B_test/2,test_axis.overlap())

#Might Need To Rewrite Inner Product Function
#for k in range(0,int(nelem)):

#print np.dot(B_cos/2,mom_evecs.T[175]) / np.dot(ax_cos.overlap(),mom_evecs.T[175])

append1=np.append(np.array([0]),mom_evecs_sin.T[169])
append1=np.append(append1,np.array([0]))

#print np.dot(B_cos/2,mom_evecs.T[175]) / np.dot(ax_cos.overlap(),mom_evecs.T[175])

for k in range(0,int(nelem)):
#	mom_eigenstates_test[k,0]=mom_eigenstates_test[k,0].real+0j
#	mom_eigenstates_test[k,180]=mom_eigenstates_test[k,180].real+0j
	mom_evecs.T[k,0]=mom_evecs.T[k,0].real+0j
	mom_evecs.T[k,180]=mom_evecs.T[k,180].real+0j

testval=np.dot(B_cos/2,mom_evecs.T[2]) / np.dot(ax_cos.overlap(),mom_evecs.T[2])
print testval
