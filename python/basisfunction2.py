#!/usr/bin/env python
import sys

import numpy as np
import my_pyplot # pyplot with include guard

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la

from math import *
import cmath
import myfunction as mf
import orthopol as op
from potential import *
from my_constants import myPi

class BasisFunction(mf.MyFunction):
    """basis functions on a fixed (possibly infinite) interval
    """

    # Jacobian function
    def jacobian(self,q):
        if   self.jc == '1': return 1.
        elif self.jc == 'q': return q
        elif self.jc == 'q^2': return q*q
        else: sys.exit('jacobian not defined "'+self.jc+'"')

    def matrix(self,name,br=None,pot='0'):
        """matrix elements between all basis functions
        |        ...int j(q)dq v^* v 
        d|d      ...int j(q)dq dv^* dv
        |1/q^2|  ...int j(q)dq v^* 1/q^2 v
        |1/q|    ...int j(q)dq v^* 1/q v
        |V|      ...int j(q)dq v^* potential v
        """
        # get quadrature grid and values and derivatives on it
        (x,w)=self.quadrature(add=10)
        v=np.zeros((len(x),self.n))
        d=np.zeros((len(x),self.n))
        for i in range(len(x)): 
            a,b=self.val(x[i],self.n)
            v[i,:]=a
            d[i,:]=b

        if br==None:
            vr=v
            dr=d
        else:
            vr=np.zeros((len(x),br.n))
            dr=np.zeros((len(x),br.n))
            for i in range(len(x)): 
                a,b=br.val(x[i],self.n)
                vr[i,:]=a
                dr[i,:]=b
            

        # local auxiliary routines
        def bracket(cv,v,f):
            m=np.zeros((np.shape(cv)[1],np.shape(v)[1]))
            for i in range(np.shape(v)[1]):
                for j in range(np.shape(cv)[1]):
                    m[i,j]=m[i,j]+np.dot(cv[:,i]*f,v[:,j])
            return m
 
        # different operators
        if   name == '|':       return  bracket(v,vr,w*self.jacobian(x))
        elif name == 'd|d':     return  bracket(d,dr,w*self.jacobian(x))
        elif name == '|dd':     return -bracket(d,dr,w*self.jacobian(x))
        elif name == '|1/q^2|': return  bracket(v,vr,w*self.jacobian(x)/x**2)
        elif name == '|q|':     return  bracket(v,vr,w*self.jacobian(x)*x)
        elif name == '|q^2|':   return  bracket(v,vr,w*self.jacobian(x)*x**2)
        elif name == '|1/q|':   return  bracket(v,vr,w*self.jacobian(x)/x)
        elif name == '1/qdq|1/qdq': 
            return  (bracket(d,dr,w*self.jacobian(x))+
                     bracket(d,vr,w*self.jacobian(x)/x)+
                     bracket(v,dr,w*self.jacobian(x)/x)+
                     bracket(v,vr,w*self.jacobian(x)/x**2))
        
        else:  return bracket(v,vr,w*self.jacobian(x)*Potential(name).v(x))

    
    # most functions have open boundary conditions
    def lb(self): return False
    def ub(self): return False

    # range: default is interval limits
    def range(self): return self.x0,self.x1

class LaguerreExpon(BasisFunction):
    """n functions L_n(2(x-x0)k) cmath.exp(-(x-x0)k) 
    """
    def __init__(self,n,k,jc='q^2',x0=0.):
        """
        n ...number of functions
        k ...exp(-x*k)
        j ...jacobian string
        """
        self.name='laguerreexpon'
        self.o=x0
        self.s=k
        self.n=n
        self.jc=jc

    def quadrature(self,add=10): 
        """base points and weights
        number of points = size of basis + add
        """
        b,w =so.la_roots(self.n+add,0.)
        for i in range(len(b)): w[i]=w[i]*cmath.exp(b[i])
        return b/(2.*self.s)+self.o, w

    def val(self,x,n=None):
        xo=x-self.o
        if n is None: n=self.n
        v,d=op.LaguerrePolynomial().val(2.*self.s*xo,n)
        return v*cmath.exp(-xo*self.s),(2.*self.s*d-v*self.s)*cmath.exp(-xo*self.s)
    
    # upper boundary is always  =0
    def ub(self): return True
    
    # range is infinite, return characteristic scale
    def range(self): return self.o,self.o+8./self.s

class FDGrid(BasisFunction):
    """n functions L_n(2x/k) cmath.exp(-x/k) 
    """
    def __init__(self,n,x0,x1,jc='1',order=1):
        """
        n.....number of points
        x0,x1.grid boundaries
        j....jacobian string
        """
        self.name='fdgrid'
        self.o=0.
        self.s=1.
        self.n=n
        self.jc=jc
        self.b=np.linspace(x0,x1,n)
        self.order=order

    def quadrature(self,add=0): return self.b,np.ones(self.n)

    def val(self,x,n=None):
        v=np.zeros(self.n)
        d=np.zeros(self.n)
        for i in range(self.n):
            if x == self.b[i]: 
                v[i]=1.
                if self.order == 1:
                    d[i]=-1./(self.b[0]-self.b[1])
                    if i > 0: d[i-1]=1./(self.b[0]-self.b[1])
                elif self.order == 2:
                    d[i]=-1.5/(self.b[1]-self.b[0])
                    if i > 0:        d[i-1]=2./(self.b[1]-self.b[0])
                    if i > 1:        d[i-2]=-0.5/(self.b[1]-self.b[0])
                elif self.order!=self.n: sys.exit('only up to order 2 or order=n (FFT)')
        if self.order==self.n:
            exit('not implemented yet')
        
        return v,d

    # Dirichlet boundary condtions
    def lb(self): return True
    def ub(self): return True

class Monomial(BasisFunction):
    """(x-(x0+x1)/2)^n basis on fixed interval [x0,x1]
    """
    # setup of the finite element
    def __init__(self,n,x0,x1,jc='1'):
        """
        n    ...highest order
        x0,x1...interval
        """
        self.name='monomial'
        self.s=(x1-x0)/2.
        self.o=(x1+x0)/2.
        self.n=n
        self.jc=jc

    # quadrature rule for monomial
    def quadrature(self,add=5): 
        b,w = op.LegendrePolynomial().quadrature(self.n+add)
        b=b*self.s+self.o
        w=w*self.s
        return b,w

    def val(self,x,n=None):
        if n is None: n=self.n
        v=np.zeros(n)
        d=np.zeros(n)
        v[0]=1.
        for i in range(1,n):
            v[i]=v[i-1]*x
            d[i]=v[i-1]*i
        return v,d

class CosSin(BasisFunction):
    """cos 2nx, sin(2n+1)x basis on fixed interval [0,2pi]
    """
    # setup of the finite element
    def __init__(self,n,jc='1'):
        """
        n    ...highest order
        x0,x1...interval
        """
        self.name='cossin'
        self.s=1.
        self.o=0.
        self.n=n
        self.jc=jc

    # quadrature rule for monomial
    def quadrature(self,add=5): 
        b = np.arange(self.n+add)*2.*pi/float(self.n+add)
        w = np.ones(len(b))*2.*pi/float(len(b))
        return b,w

    def val(self,x,n=None):
        if n is None: n=self.n
        v,d=np.zeros(n),np.zeros(n)
        for i in range(len(v)):
            if i%2 == 0:
                v[i]=( cos(x*float(i/2)))
                d[i]=(-sin(x*float(i/2))*float(i/2))
            else:
                v[i]=( sin(x*float(i/2+1)))
                d[i]=( cos(x*float(i/2+1))*float(i/2+1))
        return v,d
                

class LegendreScaled(BasisFunction):
    """Legendre polynomial P_n shifted and scaled to interval [x0,x1]
    """
    # setup of the finite element
    def __init__(self,n,x0,x1,jc='1'):
        """
        n    ...highest order
        x0,x1...interval
        """
        self.name='legendrescaled'
        self.s=(x1-x0)/2.
        self.o=(x1+x0)/2.
        self.n=n
        self.jc=jc
        self.x0=x0
        self.x1=x1

    # quadrature rule for monomial
    def quadrature(self,add=5): 
        b,w = op.LegendrePolynomial().quadrature(self.n+add)
        b=b*self.s+self.o
        w=w*self.s
        return b,w

    def val(self,x,n=None):
        if n is None: n=self.n
        v,d=op.LegendrePolynomial().val((x-self.o)/self.s,n)
        return v,d/self.s

class LegendreAssociated:
    """associated legendre functions"""
    def __init__(self,lmax,m):
        """P^(m)_l, l=m,lmax
        lmax.....number of functions
        """
        self.name='legendre associated'
        self.m=m
        self.lmax=lmax
        self.n=lmax-m+1
        self.jc='1'

    def val(self,q):
        """values and derivatives up to ORDER n (degree n-1)
        """
        def b(i,m): return (2.*float(i)-1.)/float(i-m)
        def c(i,m): return (-float(i+m)+1.)/float(i-m)
        def double_factorial(n):
            m=n
            f=1
            while (m>1): 
                f*=m
                m-=2
            return f

        # values and derivatives
        v=np.zeros((self.lmax+1))
        d=np.zeros((self.lmax+1))
        m=abs(self.m)
        v[m]=double_factorial(2*m-1)*(1-q*q)**(0.5*m)
        d[m]=double_factorial(2*m-1)*(1-q*q)**(0.5*m-1.)*m*(-q)
        if m%2==1: v[m],d[m]=-v[m],-d[m]
        for i in range(self.m+1,self.lmax+1):
            v[i]=b(i,m)*q*v[i-1]              +c(i,m)*v[i-2]
            d[i]=b(i,m)*q*d[i-1]+b(i,m)*v[i-1]+c(i,m)*d[i-2]
        if self.m<0 and m%2==1: v,d=-v,-d
        return v[m:],d[m:]

    def normY(self,l,m):
        """norm for spherical harmonics"""
        def prod(m,n): 
            p=1
            for k in range(m,n+1): p*=k
            return p

        if m>0:
            return sqrt((2.*l+1)/(prod(l-m+1,l+m))/(4*myPi))
        else:
            return sqrt((2.*l+1)*(prod(l+m+1,l-m))/(4*myPi))

def interaction(e,f,pot,add=5):
    """returns a 4-index array of interaction potential"""
    x,v=e.quadrature(add)
    y,w=f.quadrature(add)
    ve=np.zeros((len(x),e.n))
    vf=np.zeros((len(y),f.n))
    for i in range(len(x)): ve[i,:]=e.val(x[i],e.n)[0]
    for i in range(len(y)): vf[i,:]=f.val(y[i],e.n)[0]
    # Note: this is a lot of operations
    # the following arrangement of loops restricst the number of operations 
    inter=np.zeros((e.n,f.n,e.n,f.n))
    p=Potential(pot)
    for i in range(len(x)):
        je=e.jacobian(x[i])
        for j in range(len(y)):
            potij=p.pot(x[i]-y[j])*v[i]*w[j]*je*f.jacobian(y[j])
            for k in range(e.n):
                vk=ve[i,k]*potij
                for l in range(f.n):
                    vkl=vk*vf[j,l]
                    for m in range(e.n):
                        vklm=vkl*ve[i,m]
                        for n in range(f.n):
                            inter[k,l,m,n]+=vklm*vf[j,n]
    return inter

#x=LaguerreExpon(1000,0)
#print x.quadrature()
