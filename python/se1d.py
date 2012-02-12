#! /usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt

from qmoperator import *
from axis import *

# debug
import basisfunction as bf
np.set_printoptions(precision=7,suppress=True,linewidth=132)

def ho_eigen(ax,nval=5):
    hamiltonian=Operator('kinetic',ax)+Operator('+1 |q^2|',ax)
    print 'n   eigenvalue'
    for n in range(nval):
        wf=hamiltonian.eigenstate(n)
        print n,(wf | (hamiltonian*wf)) / (wf|wf)

np=5
L=8.


print "x,b"
x=bf.LegendreScaled(5,-1.,1.).quadrature(0)[0]

print "values"
for y in x: print LegendreScaled(5,-1.,1.).val(y)[0]


print 'legendre'
ho_eigen(Axis('x',np,-L/2,L/2,kind='legendre'))

print 'grid (2nd order FD)'
ho_eigen(Axis('x',np,-L/2,L/2,kind='fd'))

print 'monomial'
ho_eigen(Axis('x',np,-L/2,L/2,kind='monomial'))


print 'fem'
ho_eigen(Axis('x',40,-10.,10.,kind='fem',order=10))


print 'fd grid'
ho_eigen(Axis('x',30,-5.,5.,kind='fd',order=2))

print 'fem orders'
for n in range(2,15):
    ax=Axis('x',30,-5.,5.,kind='fem',order=n)
    print 'order,points:', n,ax.n
    ho_eigen(ax)
    

