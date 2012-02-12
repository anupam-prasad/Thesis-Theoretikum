#! /usr/bin/env python
import sys
import numpy as np
import scipy.special.orthogonal as so
from copy import *

quad01_table={}   # keep quadrature rules in this dictionary

def quad01(n,kind='gauss-legendre'):
    """get n-point quadrature rules on [0,1] and keep them in table
    
    """

    try:
        x,w=quad01_table[kind+'_'+str(n)]
    except:
        # if table entry does not exist, generate it
        if kind=='gauss-legendre': 
            x,w=so.p_roots(n)
            x=(x+1)/2
            w=w/2
        else:
            sys.exit('undefined quadrature kind="'+kind+'"')
        quad01_table[kind+'_'+str(n)]=(x,w)

    return x,w

def integral_nd(vol,integrand,params=[],nquad=4,kind="gauss-legendre",xvals=[]):
    """volume integral by recurrence
    vol    ...integration volume
    params ...extra parameters can be transmitted to the function
    nquad  ...number of quadrature points
    kind   ...type of quadrature to use
    xvals  ...for recursive integration: those values that are fixed, not integrated over
    """

   # print 'shape',np.shape(vol),'xvals',xvals
    if np.shape(vol)[0]==0: 
        return integrand(np.array(xvals),params)
    else:
        xq,wq=quad01(nquad)

        xvals.append(0.)
        for i in range(len(xq)):
            xvals[-1]=xq[i]*(vol[0,1]-vol[0,0])+vol[0,0]
            ww       =wq[i]*(vol[0,1]-vol[0,0])
            if(i==0): integrals  = ww*integral_nd(vol[1:,:],integrand,params,nquad,np.array(xvals))
            else:     integrals += ww*integral_nd(vol[1:,:],integrand,params,nquad,np.array(xvals))
        xvals.pop()
#        print vol[0,0],vol[0,1],integrals
        return integrals

def test():

    def func0(x,par=[]):  return 1/(x*x+1)
    def func1(x,par=[]):  return x*x*x
    def func01(x,par=[]): return func0(x[0])*func1(x[1]) 

    print 'd=1',integral_nd(np.array([[0.,1.]]),func0,nquad=12)
    print 'd=2',integral_nd(np.array([[0.,1.],[0.,1.]]),func01,nquad=12)

#test()
