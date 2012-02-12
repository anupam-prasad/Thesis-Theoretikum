#!/usr/bin/env python
import sys

import numpy as np
import my_pyplot # pyplot with include guard
import matplotlib.pyplot as plt

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la
from copy import *

from math import *
import cmath
import myfunction as mf
import orthopol as op
from basisfunction import *

class FiniteElement(BasisFunction):
    """complete definition of a finite element
       x0 lower boundary of element
       x1 upper boundary of element
       i0 lowest global index of element functions
       jc Jacobian for integration
       tr transformation from Legendre polynomials to FE functions
       lb0  all left boundary values = 0
       rb0 all right boundary values = 0"""

    # setup of the finite element
    def __init__(self,x0,x1,i0,nf,jc='1',lb0=False,rb0=False,fempar=None):
        """
        x0,x1  ...interval
        i0     ...global index of first function in element
        nf     ...number of functions
        jc     ...jacobian string
        lb0    ...True for lower boundary =0
        ub0    ...True for upper boundary =0
        kind   ...'legendre','laguerre'
        """
        self.par=fempar
        if   self.par is None: self.par=['legendre'];self.basis=LegendreScaled(nf,x0,x1,jc)
        elif self.par[0]=='legendre': self.basis=LegendreScaled(nf,x0,x1,jc)
        elif self.par[0]=='laguerre': self.basis=LaguerreExpon(nf,self.par[1],jc,x0)
        else: exit('no FiniteElement for function kind= '+self.par[0])

        self.l=0
        self.x0=x0
        self.x1=x1
        self.i0=i0
        self.jc=jc
        self.transformation(self.basis.n,lb0,rb0)
        self.n=np.shape(self.tr)[1]

    def __str__(self): return (self.par[0]+' '+str(self.n)+' '+str(self.x0)+' '
                               +str(self.x1)+' '+str(self.i0)+' '+self.jc)

    # quadrature rule for FE
    def quadrature(self,add=10): return self.basis.quadrature(self.n+add)

    # value of all FEM functions at x
    def val(self,x,n=None): 
        # no centrifugal term
        if self.l==0: return np.dot(self.basis.val(x),self.tr)

        # multiply by a centrifugal term
        v,d=np.dot(self.basis.val(x),self.tr)
        v*=(x/self.x1)**self.l
        d=self.l*(x/self.x1)**(self.l-1)*v+(x/self.x1)**self.l*d
        return v,d

    # True if lower/upper boundary condition = 0
    def lb(self): return self.val(self.x0)[0][0]<1.e-10
    def ub(self): return self.val(self.x1)[0][-1]<1.e-10

    def transformation(self,order,lb0,ub0):
        """set up a fem transformation matrix for general functions"""

        # 2 x n matrix of (untransformed) boundary values
        self.tr=np.eye(order)
        b=np.matrix([self.val(self.range()[0])[0],self.val(self.range()[1])[0]])

        # transformation matrix
        binv=b[:,:2].I
        self.tr=np.zeros((order,order))
        self.tr[:2,[0,-1]]=binv
        self.tr[2:,1:-1]=np.eye(order-2)
        self.tr[:2, 1:-1]=-binv*b[:,2:]
 
        # for zero boundary condition, remove first or last function
        if lb0 and not self.basis.lb(): self.tr=self.tr[:,1:]
        if ub0 and not self.basis.ub(): self.tr=self.tr[:,:-1]

        np.set_printoptions(precision=5,linewidth=120,suppress=True)

    def range(self): return self.basis.range()
    
    def overlap(self,e):
        """True if elements overlap"""
        return self.x0<e.x1 and self.x1>e.x0


    def centrifugal(self,l):
        """centrifugal factor x^l if x0=0"""
        return self
        if self.x0!=0 or l==0: return self
        b=deepcopy(self)
        b.l=l
        return b

#n=4
#x=FiniteElement(0.,1.,0,n)
#v=np.zeros((n,100))
#d=np.zeros((n,100))
#t=np.zeros((100))

#for k in range(0,100):
#	t[k]=k/100.
#	[v[:,k],d[:,k]]=x.val(k/100.)

#for k in range(0,n):
#	plt.plot(t,v[k,:])
#plt.plot(t,v[0,:])
#plt.show()
#[v,d]=x.val(.01)

#for k in range(0,10):
#	start=k/10.0
#	end=(k+1)/10.0
#	print k,start, end
#	x.append(FiniteElement(start,end,0,3))

#print x[0].tr
