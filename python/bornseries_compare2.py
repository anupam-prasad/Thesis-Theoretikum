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
bctype='r'

lb=0.
ub=8.

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
nenergy=21
Etot=np.linspace(6,8,nenergy)+2j

niter=20
a=np.zeros([nenergy,niter])
a2=np.zeros([nenergy,niter])
a_exact1=np.zeros(nenergy)
a_exact2=np.zeros(nenergy)

[evals,evecs]=la.eig(B/2+V1,y.overlap())

Gamma=10.
#Vmod=np.zeros([nelem,nelem])
#for l in range(0,int(nelem)):
#	test1=np.outer(evecs[l],evecs[l])
#	V2=V2+np.dot(test1,y.overlap())*Gamma/(evals[l]*evals[l])

#Momentum Eigenstates
k1=10
k2=5
Vk1=y.FEM_MomentumEigenstate(k1)
Ek1=k1*k1/2

Vk2=y.FEM_MomentumEigenstate(k2)
Ek1=k2*k2/2

for l in range(0,nenergy):
	E=Etot[l]
	G0_orig=la.inv(E*y.overlap()-B/2)

	Gexact1=la.inv(E*y.overlap()-B/2-V1)
	Gexact2=la.inv(E*y.overlap()-B/2-V2)

	a_exact1[l]=la.norm(Gexact1)
	a_exact2[l]=la.norm(Gexact2)

	A=G0_orig
	A2=G0_orig

	tempmat1=np.dot(G0_orig,V1)
	tempmat2=np.dot(G0_orig,V2)

#	Aex=la.inv(E*y.overlap()-B/2-V)
#	aex=la.norm(Aex)
	
	for k in range(0,niter):
		Aprev=A
		A=G0_orig+np.dot(tempmat1,A)
		a[l][k]=la.norm(A-Aprev)
#		a[l][k]=la.norm(A)
		
		Aprev2=A2
		A2=G0_orig+np.dot(tempmat2,A2)
		a2[l][k]=la.norm(A2-G0_orig)
#		a2[l][k]=la.norm(A2)

	#print A
	#print Aex

	#print a

	#print la.norm(a)


f=open('bornseries/comparison_values3','w')
pickle.dump([a,a2],f)
#print evals

t=range(0,niter)
#print len(a[0])

#print a2
#raw_input()

for l in range(0,int(nelem)):
	mat1=np.outer(evecs[l],evecs[l])
	mat2=mat1*y.overlap()
	print mat2.sum()

for l in range(0,nenergy):
	plt.plot(t,a[l], 'bs')
	plt.plot(t,a[l], 'blue')

	plt.plot(t,a2[l],'rs')
	plt.plot(t,a2[l],'red')

	plt.xlabel('Number of Iterations')
	plt.ylabel('norm(G_n-G_(n-1))')
	ptitle='Born Series Convergence (Green\'s Function)'
	plt.title(ptitle)
	ypos=np.mean([a[l][0], a2[l][0]])
	plt.text(15, ypos, 'Energy = '+str(Etot[l]))
	plt.grid(True)
#	plt.plot(t,np.ones(niter)*a_exact1[l])
#	plt.plot(t,np.ones(niter)*a_exact2[l])
	plt.show()
	time.sleep(.5)
	plt.clf()
