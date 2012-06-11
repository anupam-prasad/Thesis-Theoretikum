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

n=500
order=6

lb=-50.
ub=50.

y_cos=Axis('xopen',n,lb,ub,'fem',order)
y_sin=Axis('x',n,lb,ub,'fem',order)

#d|d matrix
B=np.zeros([y_sin.len()+y_cos.len(),y_sin.len()+y_cos.len()])

overlap=np.zeros([y_sin.len()+y_cos.len(),y_sin.len()+y_cos.len()])
overlap_inv=np.zeros([y_sin.len()+y_cos.len(),y_sin.len()+y_cos.len()])

#Potential Matrix
V1=np.zeros([y_cos.len()+y_sin.len(),y_cos.len()+y_sin.len()])
V2=np.zeros([y_cos.len()+y_sin.len(),y_cos.len()+y_sin.len()])

V0=1

for n in range(len(y_cos.e)):
	i0=y_cos.e[n].i0
	i1=i0+y_cos.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y_cos.e[n].matrix('d|d')
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y_cos.e[n].matrix('fwell', np.array([-5.,5.,V0]))
#	V2[i0:i1,i0:i1]=V2[i0:i1,i0:i1]+y_cos.e[n].matrix('gaussian', np.array([6.,7.,V0]))

for n in range(len(y_sin.e)):
	i0=y_sin.e[n].i0+y_cos.len()
	i1=i0+y_sin.e[n].n
	B[i0:i1,i0:i1]=B[i0:i1,i0:i1]+y_sin.e[n].matrix('d|d')
	V1[i0:i1,i0:i1]=V1[i0:i1,i0:i1]+y_sin.e[n].matrix('fwell', np.array([-5.,-5.,V0]))
#	V2[i0:i1,i0:i1]=V2[i0:i1,i0:i1]+y_sin.e[n].matrix('gaussian', np.array([6.,7.,V0]))

overlap[:y_cos.len(),:y_cos.len()]=y_cos.overlap()
overlap[y_cos.len():,y_cos.len():]=y_sin.overlap()

overlap_inv[:y_cos.len(),:y_cos.len()]=y_cos.overlap()
overlap_inv[y_cos.len():,y_cos.len():]=y_sin.overlap()

[mom_evals_temp,mom_evecs_temp]=la.eig(B/2.,overlap)
[evals,evecs]=la.eig(B/2. + V1 + V2,overlap)

#Sorting eigenvalues/eigenvectors in ascending order
perm=np.argsort(mom_evals_temp)
mom_evals_temp=mom_evals_temp[perm]
mom_evecs_temp=mom_evecs_temp[:,perm]

perm=np.argsort(evals)
evals=evals[perm]
evecs=evecs[:,perm]

mom_evecs=np.zeros([y_sin.len()+y_cos.len(),(y_sin.len()+y_cos.len())/2])+0j
mom_evals=np.zeros([(y_sin.len()+y_cos.len())/2])+0j

mom_evecs[:,0]=mom_evecs_temp[:,0]
for l in range(1,(y_sin.len()+y_cos.len())/2):
	mom_evecs[:,l]=mom_evecs_temp[:,2*l-1]+mom_evecs_temp[:,2*l]
	mom_evecs[:,l]=np.sqrt(ub-lb) * mom_evecs[:,l]/np.sqrt(np.dot(mom_evecs[:,l],np.dot(overlap,mom_evecs[:,l])))
	mom_evals[l]=mom_evals_temp[2*l]

eps=1j
Bmod=np.dot(B/2.,overlap_inv)

for l in range(2):
	Emat=np.eye(y_cos.len()+y_sin.len()) * (mom_evals[l]+eps)
	G0temp = la.inv(-Bmod+Emat)
	G0 = np.dot(G0temp,overlap)
#	print G0
	raw_input()
	tempvec1=np.dot(G0,mom_evecs[:,l])
	print abs(np.dot(mom_evecs[:,l],tempvec1)) / (ub-lb)
