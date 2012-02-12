import sys

import numpy
from numpy import *
import matplotlib.pyplot as plt

import scipy.linalg as la

import __main__ 

class TensorProd:

    def __init__(self,a,b):
        self.a=a
        self.b=b

    def dim(self): 
        return (shape(self.a)[0]*shape(self.b)[0],shape(self.a)[1]*shape(self.b)[1])

    def array(self):
        """
        convert tensor product to array
        """
        a=zeros(self.dim())
        ij0=-1
        for i0 in range(shape(self.a)[0]):
            for j0 in range(shape(self.b)[0]):
                ij0=ij0+1
                ij1=-1
                for i1 in range(shape(self.a)[1]):
                    for j1 in range(shape(self.b)[1]):
                        ij1=ij1+1
                        a[ij0,ij1]=self.a[i0,i1]*self.b[j0,j1]
        return a

class Tensor:
    """
    linear combination of tensor products and multi-index arrays
    """
    def __init__(self,t,c=1.):
        self.t=[]
        self.t.append(t)
        self.t.append(c)

    def copy(self):
        s=Tensor(self.t[0],self.t[1])
        for i in range(2,len(self.t)):
            s.t.append(self.t[i])
        return s

    def __add__(self,a):

        # check dimensions
        if not self.dim() == a.dim():
            sys.exit('dimensions do not match')

        s=self.copy()
        for i in range(0,len(a.t)):
            ati=a.t[i]
            s.t.append(ati)
        return s

    def dim(self):
        t=self.t[0]
        if type(t) == numpy.ndarray: 
            return shape(t)
        elif t.__class__ == __main__.TensorProd:
            return t.dim()
        else: sys.exit('no dimensions defined for component')

    def array(self):
        """
        convert tensor into single array
        """
        a=zeros(self.dim())
        for i in range(0,len(self.t),2):
            t=self.t[i]
            if t.__class__ == __main__.TensorProd: a=a+t.array()*self.t[i+1]
            elif type(t) == numpy.ndarray:         a=a+t*self.t[i+1]
            else:
                print type(t),t.__class__
                sys.exit('unidentified tensor component')
        return a

def test():
    s=Tensor(identity(4),3.)
    t=s+Tensor(identity(4),1.5)
    
    print 's:',s.array()
    print 't:',t.array()
    
    p=t+Tensor(TensorProd(identity(2),identity(2)),3.)
    print 'p:',p.array()

    x=s+Tensor(TensorProd(identity(2),identity(3)))

