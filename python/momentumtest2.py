#!/usr/bin/env python
import sys

import numpy as np

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la

from axis import *

n=100
order=20
bctype='xopen'

lb=-1.
ub= 1.

ncos=1

ax_cos=Axis(bctype,n,lb,ub,'fem',order)

#d|d matrix
B_cos=np.zeros([ax_cos.len(),ax_cos.len()])
B_test=np.zeros([ax_cos.len()-2,ax_cos.len()-2])

iter1=0
for e in ax_cos.e:
        b=e.matrix('d|d')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B_cos[iter1+k1,iter1+k2]=B_cos[iter1+k1,iter1+k2]+b[k1][k2]
        iter1=iter1+iter2-1

#Scattering Energy
Etot=np.linspace(0,ax_cos.len()-1,ax_cos.len())*myPi+0j

#Momentum Eigenstates - Not sure if this works
momentum_eigenstates=np.zeros([ax_cos.len(),ax_cos.len()])+0j
print "sizes",ax_cos.len(),momentum_eigenstates.size,Etot.size
for k in range(0,int(ax_cos.len())):
	momentum_eigenstates[k,:]=ax_cos.FEM_function(np.exp,Etot[k]*1j)

#[mom_evals,mom_evecs]=la.eig(B_cos/2,ax_cos.overlap())

#print mom_evals
#print mom_evecs

#ax_cos.plot()
#plt.show()
#raw_input()
#for k in range(0,int(ax_cos.len())):
#	ax_cos.FEM_functionplot(mom_evecs.T[k])
#	ax_cos.FEM_functionplot(momentum_eigenstates[k].imag)
count=0
for k in range(0,10):
	for l in range(0,10):
		inner=ax_cos.FEM_InnerProduct(momentum_eigenstates[k],momentum_eigenstates[l])
		if k!=l and (abs(inner.imag)>1e-8 or abs(inner.real)>1e-8):
			print inner, k,l
			count=count+1
print count
#for k in range(0,int(ax_cos.len())):
#	ax_cos.FEM_plot(momentum_eigenstates[k].real)
