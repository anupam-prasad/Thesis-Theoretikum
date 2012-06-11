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

n=100
order=5
bctype='x'

lb=-4.
ub=4.

y=Axis(bctype,n,lb,ub,'fem',order)#+Axis('x',n,lb,ub,'fem',order)#+Axis('r',2,2.,3.,'fem')

#y.plot()
#plt.show()
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
V1=np.zeros([nelem,nelem])
V2=np.zeros([nelem,nelem])

#Control Parameter
Lambda=1

iter1=0

#Coefficient Matrix
U=np.zeros(nelem)

for e in y.e:
#       for b in e.b:
        b=e.matrix('d|d')
        v=e.matrix('pot')
	v2=e.matrix('pot2')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
                        V1[iter1+k1,iter1+k2]=V1[iter1+k1,iter1+k2]+v[k1][k2]
                        V2[iter1+k1,iter1+k2]=V2[iter1+k1,iter1+k2]+v2[k1][k2]
        iter1=iter1+iter2-1

#Scattering Energy
E=2+1j

niter=40
a=np.zeros(niter)

[evals,evecs]=la.eig(B/2+V1,y.overlap())

G0=la.inv(E*y.overlap()-B/2-V2)
A=G0
tempmat=np.dot(G0,V1-V2)

Aex=la.inv(E*y.overlap()-B/2-V1)
aex=la.norm(Aex)

for k in range(0,niter):
	Aprev=A
	A=G0+np.dot(tempmat,A)
	a[k]=la.norm(A-Aprev)

	#print A
	#print Aex

	#print a

	#print la.norm(a)


#f=open('bornseries/values3','w')
##pickle.dump(a,f)
#print evals

#t=range(0,niter)
#print len(a[0])
#for l in range(0,nenergy):
#	plt.plot(t,a[l])
#	plt.show()
#	time.sleep(.5)
#	plt.clf()

evals_sorted=np.sort(abs(evals))

x=np.linspace(0,nelem-1,int(nelem))

#print abs(x-evals_sorted)

#tempvec1=np.dot(y.overlap(),evecs[:,87])
#tempvec2=np.dot(Aex,tempvec1)

#print tempvec2 / evecs[:,87]

#print abs(evals[87])
