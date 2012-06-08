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
ub=20.

y_cos=Axis('xopen',n,lb,ub,'fem',order)
y_sin=Axis('x',n,lb,ub,'fem',order)

#d|d matrix
B_cos=np.zeros([y_cos.len(),y_cos.len()])
B_sin=np.zeros([y_sin.len(),y_sin.len()])

B=np.zeros([y_sin.len()+y_cos.len(),y_sin.len()+y_cos.len()])

overlap=np.zeros([y_sin.len()+y_cos.len(),y_sin.len()+y_cos.len()])
overlap_inv=np.zeros([y_sin.len()+y_cos.len(),y_sin.len()+y_cos.len()])

#Potential Matrix
V1_cos=np.zeros([y_cos.len(),y_cos.len()])
V2_cos=np.zeros([y_cos.len(),y_cos.len()])

V1_sin=np.zeros([y_sin.len(),y_sin.len()])
V2_sin=np.zeros([y_sin.len(),y_sin.len()])

V=np.zeros([y_cos.len()+y_sin.len(),y_cos.len()+y_sin.len()])

V0=1

for n in range(len(y_cos.e)):
	i0=y_cos.e[n].i0
	i1=i0+y_cos.e[n].n
	B_cos[i0:i1,i0:i1]=B_cos[i0:i1,i0:i1]+y_cos.e[n].matrix('d|d')
	V1_cos[i0:i1,i0:i1]=V1_cos[i0:i1,i0:i1]+y_cos.e[n].matrix('fwell', np.array([9.,11.,V0]))
#	V2_cos[i0:i1,i0:i1]=V2_cos[i0:i1,i0:i1]+y_cos.e[n].matrix('gaussian', np.array([6.,7.,V0]))


for n in range(len(y_sin.e)):
	i0=y_sin.e[n].i0
	i1=i0+y_sin.e[n].n
	B_sin[i0:i1,i0:i1]=B_sin[i0:i1,i0:i1]+y_sin.e[n].matrix('d|d')
	V1_sin[i0:i1,i0:i1]=V1_sin[i0:i1,i0:i1]+y_sin.e[n].matrix('fwell', np.array([9.,11.,V0]))
#	V2_sin[i0:i1,i0:i1]=V2_sin[i0:i1,i0:i1]+y_sin.e[n].matrix('gaussian', np.array([6.,7.,V0]))

B[:y_cos.len(),:y_cos.len()]=B_cos
B[y_cos.len():,y_cos.len():]=B_sin

overlap[:y_cos.len(),:y_cos.len()]=y_cos.overlap()
overlap[y_cos.len():,y_cos.len():]=y_sin.overlap()

overlap_inv[:y_cos.len(),:y_cos.len()]=y_cos.overlap()
overlap_inv[y_cos.len():,y_cos.len():]=y_sin.overlap()

V[:y_cos.len(),:y_cos.len()]=V1_cos + V2_cos
V[y_cos.len():,y_cos.len():]=V1_sin + V2_sin

[mom_evals_temp,mom_evecs_temp]=la.eig(B/2.,overlap)
[evals,evecs]=la.eig(B/2. + V,overlap)

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(mom_evals_temp)
mom_evals_temp=mom_evals_temp[perm]
mom_evecs_temp=mom_evecs_temp[:,perm]

perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

mom_evecs=np.zeros([y_sin.len()+y_cos.len(),(y_sin.len()+y_cos.len())/2])
mom_evals=np.zeros([(y_sin.len()+y_cos.len())/2])

mom_evecs[:,0]=mom_evecs_temp[:,0]
for l in range(1,(y_sin.len()+y_cos.len())/2):
	mom_evecs[:,l]=mom_evecs_temp[:,2*l-1]+mom_evecs_temp[:,2*l]

for k in range(10):
	print np.dot(mom_evecs[:,k],np.dot(overlap,mom_evecs[:,k]))

y_cos.FEM_plot(mom_evecs[:y_cos.len(),1])
y_sin.FEM_plot(mom_evecs[y_cos.len():,1])
exit('here')
