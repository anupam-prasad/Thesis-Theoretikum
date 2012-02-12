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
        v=e.matrix('fwell')
#	v2=e.matrix('pot2')
        iter2=int(np.sqrt(np.size(b)))
        for k1 in range(0,iter2):
                for k2 in range(0,iter2):
                        B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
                        V1[iter1+k1,iter1+k2]=V1[iter1+k1,iter1+k2]+v[k1][k2]
#                       V2[iter1+k1,iter1+k2]=V2[iter1+k1,iter1+k2]+v2[k1][k2]
        iter1=iter1+iter2-1

#Scattering Energy
nenergy=21
momentum=np.linspace(0,10,nenergy)
Etot=momentum*momentum/2+.1j

niter=15
a=np.zeros([nenergy,niter])
a2=np.zeros([nenergy,niter])
a_exact=np.zeros(nenergy)

b=np.zeros([nenergy,niter])
b2=np.zeros([nenergy,niter])
b_exact=np.zeros(nenergy)
[evals,evecs]=la.eig(B/2+V1,overlap)

alpha=.99

#V2=alpha*V1+(1-alpha)*V1
#Vmod=np.zeros([nelem,nelem])
#I=np.identity(nelem,float)
#for k in range(0,int(nelem)):
#	for l in range(0,int(nelem)):
#	test1=np.outer(evecs.T[k],evecs.T[k])
#	test1=test1.T
#	test2=test1*overlap
#	norm_factor=test2.sum()
#	V2=V2+evals[k]*np.dot(test1,overlap) / norm_factor

#raw_input()
for l in range(0,nenergy):
	E=Etot[l]
	G0_orig=la.inv(E*overlap-B/2)
	G0_mod=la.inv(E*overlap-B/2-alpha*V1)
	Gexact=la.inv(E*overlap-B/2-V1)

	a_exact[l]=la.norm(Gexact)

	A=G0_orig
	A2=G0_mod

#	tempmat1=np.dot(G0_orig,V1)
#	tempmat2=np.dot(G0_mod,V1-V2)

	tempmat1=np.dot(G0_orig,V1)
	tempmat2=np.dot(G0_mod,(1-alpha)*V1)

	storemat1=G0_orig
	storemat2=G0_mod
#	Aex=la.inv(E*overlap-B/2-V)
#	aex=la.norm(Aex)
	
	for k in range(0,niter):
		b[l][k]=la.norm(A)
		Aprev=A
		A=G0_orig+np.dot(tempmat1,A)
#		storemat1=np.dot(storemat1,tempmat1)
#		A=A+storemat1
		a[l][k]=la.norm(A-Aprev)
		
		b2[l][k]=la.norm(A2)
		Aprev2=A2
		A2=G0_mod+np.dot(tempmat2,A2)
#		storemat2=np.dot(storemat2,tempmat2)
#		A=A+storemat2
		a2[l][k]=la.norm(A2-Aprev2)

	#print A
	#print Aex

	#print a

	#print la.norm(a)


f=open('bornseries/comparison_values','w')
pickle.dump([a,a2,b,b2,a_exact],f)
#print evals

t=range(1,niter+1)
#print len(a[0])

#print a2
#raw_input()

for l in range(0,nenergy):
	plt.figure(1)

	plt.subplot(211)
#	plt.axis([0,niter+1,-.5, (np.mean(b2[l])+1.)])
#	plt.plot(t,b[l],'bo',t,b2[l],'rs',t,np.ones(niter)*a_exact[l],'green')
#	plt.plot(t,b[l],'b--',t,b2[l],'r--')
	plt.plot(t,b2[l],'r^',t,np.ones(niter)*a_exact[l],'green')
        plt.ylabel('norm(G_n)')
        ptitle='Born Series Convergence (Green\'s Function Norm)'
        plt.title(ptitle, fontsize=13)
        ypos=np.mean([b[l][0], b2[l][0]])
#        plt.text(8, ypos, 'Energy = '+str(Etot[l]))
        plt.grid(True)
	
	plt.subplot(212)
	plt.axis([0,niter+1,-.5, max(a2[l])+1.])
#	plt.plot(t,a[l], 'bo',t,a2[l],'rs',t,np.zeros(niter),'green')
#	plt.plot(t,a[l],'b--',t,a2[l],'r--')
	plt.plot(t,a2[l],'r^',t,np.zeros(niter),'green')
	plt.ylabel('norm(G_n-G_(n-1))')
        plt.xlabel('Number of Iterations. E = '+str(Etot[l]))
        ptitle='Born Series Convergence (Green\'s Function Difference)'
        plt.title(ptitle, fontsize=13)
        ypos=np.mean([a[l][0], a2[l][0]])
#       plt.text(8, ypos, 'Energy = '+str(Etot[l]))
        plt.grid(True)
	

	plt.show()
	time.sleep(.5)
	plt.clf()
