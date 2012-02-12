#!/usr/bin/env python
"""general 1d MC integrator"""
import sys 
import numpy as np
import scipy.linalg as linalg
import scipy.linalg.flapack as lapack
import mytimer
import copy
import my_pyplot
from math import *

def mcintegral(a,b,func,eps,maxpoint=1e9):
    """MC integration over [a,b]"""

    np.random.seed(1)
    prev=np.zeros((3))
    Mmod=1
    icur=0
    mcsum=0.
    N=0
    while N<int(maxpoint):
        if N%Mmod==0:
            Mmod*=2
            icur=(icur+1)%len(prev)

            # primitive check for convergence
            res=mcsum*(b-a)/(N+1)
            if all(abs(prev-res)<eps*abs(res)): return res

            # store result
            prev[icur]=res

            # get new random numbers
            x=a+(b-a)*np.random.random((Mmod))

        mcsum=mcsum+func(x[N%Mmod])
        N+=1

    else:
        print 'maximum number of '+int(maxpoint)+' points exceeded'

if __name__  == "__main__":
    """standalone: tests"""
    a=0.
    b=1.
    def func(x): return x*x

    print mcintegral(a,b,func,1.e-3)#/(1.-exp(-1.))
