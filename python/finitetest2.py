#!/usr/bin/env python
import sys

import numpy as np
#import my_pyplot
#import matplotlib as plt

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la
from copy import *

#from math import *
#import cmath

#from array import array

from axis import *

def pot(x):
	V=np.zeros(len(x))
	for k in range(len(x)):
		if x[k] < 1. and x[k] > 0.:
#                       V[k]=-10
#			V[k]=-np.exp(-x[k]*x[k]/.01)
			V[k]=-1/(x[k])
		else:
			V[k]=0

	return V

n=4
order=2
bctype='rho'

lb=0.
ub=1.
	
y=Axis(bctype,n,lb,ub,'fem',order)#+Axis('x',9,0.,3.,'fem',4)#+Axis('r',2,2.,3.,'fem')

npart=np.floor(n/(order-1))

if bctype=='cos':
	nelem=npart*(order-1)+1
elif bctype=='r' or bctype=='rho':
	nelem=npart*(order-1)
else:
	nelem=npart*(order-1)-1

#print npart,nelem
plt.clf()

overlap=y.overlap()
overlap_inv=y.overlap_inv()

#d|d matrix
B=np.zeros([nelem,nelem])
#Potential Matrix
V=np.zeros([nelem,nelem])

iter1=0

#Coefficient Matrix
U=np.zeros(nelem)

for e in y.e:
#	for b in e.b:
	b=e.matrix('d|d')
	v=e.matrix('pot')
	iter2=int(np.sqrt(np.size(b)))
	for k1 in range(0,iter2):
		for k2 in range(0,iter2):
			B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
			V[iter1+k1,iter1+k2]=V[iter1+k1,iter1+k2]+v[k1][k2]
	iter1=iter1+iter2-1

#print V
#[P,L,U]=la.lu(np.dot(overlap_inv,B))
[evals,evecs]=la.eig(B+V,overlap)

evals=np.abs(evals)

#for k in range(0,int(nelem)):
#	for l in range(k+1,int(nelem)):
#		if evals[k]>evals[l]:
#			temp1=evals[k]
#			evals[k]=evals[l]
#			evals[l]=temp1
#			temp2=evecs[k]
#			evecs[k]=evecs[l]
#			evecs[l]=temp2
#axis vector
npoints=100*npart
t=np.linspace(lb,ub,npoints)

#print t[100]
#raw_input()
#Wavefunction vector
wavefun=np.zeros([nelem,npoints])

iter3=0
iter5=0
lenstore1=0
lenstore2=0
k=0

for e in y.e:
	lb=e.x0
	ub=e.x1
#	print e.matrix('|V|',pot)
	for k in range(iter3*100,(iter3+1)*100):
		v,d=e.val(t[k])
#		if iter3==0:
#			iter5=0
#		else:
#			iter5=iter3*(order-1)-1
#		if lenstore != len(v):
			
		lenstore=len(v)

		for iter4 in range(0,len(v)):
			for l in range(0,int(nelem)):
	 			wavefun[l][k]=wavefun[l][k]+evecs[iter5][l]*v[iter4]
			iter5=iter3*(order-1)+iter4+1
		iter5=iter5-1
	iter3=iter3+1

#y.plot()
#plt.show()
#print evecs
#for l in range(0,int(nelem)):
#	print evecs[0][l]
#raw_input()

for k in range(0,int(nelem)):
	plt.plot(t,wavefun[k])
#	plt.plot(t,pot(t))
	plt.show()
	raw_input()
	plt.clf()
