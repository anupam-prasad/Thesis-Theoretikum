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
	
y=Axis('cos',1,0.,1.,'fem',1)#+Axis('x',9,0.,3.,'fem',4)#+Axis('r',2,2.,3.,'fem')

#y=Axis('r',4,0.,1.,'fem')#+Axis('r',2,1.,2.,'fem')

#print 
t=np.linspace(-5,5,1000)
V=pot(t)

#v,d=y.val(.01)

#plt.clf()
#plt.plot(t,V)
#plt.show()

#print y.overlap()

y.plot()
plt.show()

#a=np.random.rand(4,4)
a=y.overlap()

[P,L,U]=la.lu(a)

#print y.overlap_inv()

for e in y.e:
#	for b in e.b:
	print e.matrix('d|d')

print y.overlap()
