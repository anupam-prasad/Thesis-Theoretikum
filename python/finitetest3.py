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
from potential2 import *

#from math import *
#import cmath

#from array import array

from axis import *

n=4
order=2
bctype='x'

lb=-5.
ub=5.
	
y=Axis(bctype,n,lb,ub,'fem',order)#+Axis(bctype,n,lb,ub,'fem',order)#+Axis('r',2,2.,3.,'fem')

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
[evals,evecs]=la.eig(B/2+V,overlap)

#print V
#raw_input()

#evals=np.abs(evals)
evalstest=np.sort(evals)
#evalstest=evalstest
#print evalstest
n1=np.linspace(1,len(evalstest),len(evalstest))

#print evalstest

print evecs

V = np.zeros([nelem,nelem])
for column1 in evecs.T:
	for column2 in evecs.T:
		temp1=np.outer(column1,column2)
		temp1=temp1*overlap
#		print temp1.sum()

#axis vector
npoints=100*npart
t=np.linspace(lb,ub,npoints)

#print t[100]
#raw_input()
#Wavefunction vector
wavefun=np.zeros([nelem,npoints])

iter3=0
iter5=0
k=0

#for e in y.e:
#	lb=e.x0
#	ub=e.x1
#	for k in range(iter3*100,(iter3+1)*100):
#		v,d=e.val(t[k])
#		if iter3==0:
#			iter5=0
#		else:
#			iter5=iter3*(order-1)-1
#		if lenstore != len(v):
#      			lenstore=len(v)

#		for iter4 in range(0,len(v)):
#			for l in range(0,int(nelem)):
#	 			wavefun[l][k]=wavefun[l][k]+evecs.T[iter5+iter4][l]*v[iter4]
#	iter5=iter5+len(v)-1
#	iter3=iter3+1

#y.plot()
#plt.show()
#print evecs
#for l in range(0,int(nelem)):
#	print evecs[0][l]
#raw_input()

#t1=np.linspaci
