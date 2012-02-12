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
		if x[k] < 1.:
			V[k]=0
		elif x[k] > 0.:
			V[k]=0
		else:
			V[k]=myInf

#	V=-np.exp(-x*x/(2*myPi*5))
	return V
	
n=4
#x=[]

#x.append(FiniteElement(0.,1.,0,n))
#x.append(FiniteElement(1.,2.,0,n))

#v1=np.zeros((n,100))
#d1=np.zeros((n,100))
#
#v2=np.zeros((n,100))
#d2=np.zeros((n,100))
#
#t1=np.zeros((100))
#t2=np.zeros((100))
#
#for k in range(0,100):
#        t1[k]=k/100.
#        t2[k]=(k/100.)+1
#        [v1[:,k],d1[:,k]]=x[0].val(k/100.)
#        [v2[:,k],d2[:,k]]=x[1].val((k/100.)+1)

#for k in range(0,n):
#        plt.plot(t1,v1[k,:])
#	plt.plot(t2,v2[k,:])
#plt.plot(t,v[0,:])
#plt.show()
#[v,d]=x.val(.01)

#print t2
#for k in range(0,10):
#       start=k/10.0
#       end=(k+1)/10.0
#       print k,start, end
#       x.append(FiniteElement(start,end,0,3))

#print x[0].tr

y=Axis('rho',3,0.,1.,'fem',2)#+Axis('x',9,0.,3.,'fem',4)#+Axis('r',2,2.,3.,'fem')

#y=Axis('r',4,0.,1.,'fem')#+Axis('r',2,1.,2.,'fem')

for e in y.e:
#       for b in e.b:
        b=e.matrix('d|d')
        v=e.matrix('pot')
        iter2=int(np.sqrt(np.size(b)))
	print iter2
#	print b
	for a in b:
		print a
#        for k1 in range(0,iter2):
#               for k2 in range(0,iter2):
#                      B[iter1+k1,iter1+k2]=B[iter1+k1,iter1+k2]+b[k1][k2]
#                        V[iter1+k1,iter1+k2]=V[iter1+k1,iter1+k2]+v[k1][k2]
#        iter1=iter1+iter2-1

raw_input()
#print 
t=np.linspace(-5,5,1000)
V=pot(t)

#v,d=y.val(.01)
#plt.clf()
#plt.plot(t,V)
#plt.show()

#print y.overlap()

y.plot()
#plt.show()

#a=np.random.rand(4,4)

#[P,L,U]=la.lu(a)

#A=np.dot(L,U)

#print y.show('blah')
#print y.range()
#for e in y.e:
#	print e.val(.01)

print y.overlap()
plt.show()
