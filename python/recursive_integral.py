#! /usr/bin/env python
import sys
import numpy as np
import scipy.special.orthogonal as so
from copy import *
from integral_nd import *

max_level=10      # allow recursions only to this level

def vol_sub(vol,n):
    """return n'th subvolume"""
    sub=deepcopy(vol)
    m=n
    for p in range(np.shape(vol)[0]):
        sub[p,m%2]=(vol[p,0]+vol[p,1])/2
        m=m/2
    return  sub

def recursive_integral(vol,integrand,params=[],acc_rel=1.e-12,acc_abs=1.e-12,
                       nquad=8,value_previous=None,level=0):
    """"
    vol        ...integration volume (lower,upper bounds)
    integrand  ...integrand function
    acc_rel    ...relative accuracy
    acc_abs    ...abssolute accuracy
    level      ...never specify, internal use only
    """
#    print "level",level,vol

    # check recursive level
    if level > max_level: exit('recursive integration level exceeded '+str(max_level))

    # first entry, quadrature over full interval for comparison
    if value_previous==None: value_previous=integral_nd(vol,integrand,params,nquad,[])

    # integrate over split volumes 
    val_sub=[]
    for i in range(2**np.shape(vol)[0]): 
        val_sub.append(integral_nd(vol_sub(vol,i),integrand,params,nquad,[]))
    value=np.zeros(np.shape(val_sub[0]))
    for i in range(len(val_sub)): value+= val_sub[i]

    # convergence criterion: relative or absolute error
    if max(abs((value-value_previous)/value))<acc_rel or  max(abs((value-value_previous)))<acc_abs: 
        return value

    # get full precision integrals on each sub-interval
    value=0
    for i in range(2**np.shape(vol)[0]): 
        value+=recursive_integral(vol_sub(vol,i),integrand,params,acc_rel,acc_abs/2**np.shape(vol)[0],\
                                      nquad,val_sub[i],level=level+1)

    level=0 # before returning reset level
    return value

def test():

    def func0(x,par=[]): return np.array([1.,1/(x*x+1)])
    def func1(x,par=[]): return np.array([x*x*x])
    def func01(x,par=[]): return np.array([1.,x[0]*x[1],x[0]/(x[1]*x[1]+1)])

    print 'd=1',recursive_integral(np.array([[0.,1.]]),func0,nquad=8)
    print 'd=2',recursive_integral(np.array([[0.,1.],[0.,1.]]),func01,nquad=8)

test()


