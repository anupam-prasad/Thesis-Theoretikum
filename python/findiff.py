# collection of finite difference methods
import sys
import numpy as np
import scipy.linalg as la

methods={                           # dx^2 schemes with increasing order
    'd2_o2': (-2.,1.),              # define method as "tuple" of diagonal and off-diagonal values
    'd2_o4': (-5./2.,4./3.,-1./12.),
    'd2_o6': (-49./18.,3./2.,-3./20.,1./90.),
    'd2_o8': (-205./72.,8./5.,-1./5.,8./315.,-1./560.)
   }

def method(n,order=2,deriv=2):                               # return a matrix for the FD scheme
        m=np.zeros((n,n))                                    # create 0-matrix
        noff=len(methods['d'+str(deriv)+'_o'+str(order)])-1  # number of off-diagonals in method = length of definig "tuple" -1
        d=np.zeros((2*noff+1))                               # create auxiliary vector    
        d[noff:]=methods['d'+str(deriv)+'_o'+str(order)]     # upper half = method tuple   
        d[noff::-1]=d[noff:]                                 # lower half = reversed method tuple
        
        for i in range(n):                                   # put into matrix (watch out for index limits)
            m[i,max(0,i-noff):min(n,i+noff+1)]=d[max(noff-i,0):min(noff+n-i,2*noff+1)]
        return m



def general(x,n):
    """finite-difference derivative scheme for general grid
    x..... grid points
    n..... number of points in scheme
    tries to center the scheme, switches to backwar/forward at the upper lower edge
    """
    
    # need n+1 points for n-th order scheme
    if len(x)<n+1: sys.exit('cannot do '+str(n)+' order scheme with only '+str(len(x))+' points')
    
    # for too closely lying values we will get numerical problems
    for i in range(len(x)):
        for j in range(i):
            if abs(x[j]-x[i])<1.e-10*(max(x)-min(x)): sys.exit('x-values too close')


    fd=np.zeros((len(x),len(x))) 
    for i in range(len(x)):

        # center interval (as well as possible)
        i0=max(i-n/2,0)        
        i1=min(i0+n+1,len(x))
        i0=min(i0,i1-n-1)

        # get powers of grid points
        xp=np.zeros((n+1,n+1)) 
        xp[0,:]=np.ones((n+1))
        xp[1,:]=x[i0:i1]-sum(x[i0:i1])/(n+1) # evaluate at shifted points (reduce floating errors)
        for j in range(2,n+1): xp[j,:]=xp[j-1,:]*xp[1,:] # higher powers of x_i
        dp=np.zeros((n+1))                    
        for j in range(1,n+1): dp[j]=float(j)*xp[j-1,i-i0] # d/dx x_i^j = j x_i^j

        # solve for coefficients
        fd[i,i0:i1]=la.solve(xp,dp)

    return fd
