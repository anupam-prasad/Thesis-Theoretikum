#!/usr/bin/env python
import sys

import numpy as np

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la

from axis import *

n=300
order=26
bctype='xopen'

lb=0.
ub=10.

ax_cos=Axis(bctype,n,lb,ub,'fem',order)

#Scattering Energy
Ptot=np.linspace(0,ax_cos.len()-1,ax_cos.len())*myPi+0j

#Momentum Eigenstates
momentum_eigenstates=np.zeros([ax_cos.len(),ax_cos.len()])+0j

for k in range(0,int(ax_cos.len())):
	momentum_eigenstates[k]=ax_cos.FEM_function(np.exp,Ptot[k]*1j)

#ax_cos.FEM_plot(momentum_eigenstates[1].real)
count=0
inner=np.zeros(190)

for k in range(0,20):
	for l in range(k+1,20):
		inner[count]=abs(ax_cos.FEM_InnerProduct(momentum_eigenstates[k],momentum_eigenstates[l]))
		count=count+1

print np.mean(inner), np.std(inner), max(inner), min(inner)
