#! /usr/bin/env python
import numpy as np
import scipy.linalg as la
import findiff as fd  # simple home-made module for FD schemes

# "input"
ng=64        # number of grid points
L=10.        # size of simulation box
o=8    # order of the FD scheme

def potential(x): return x*x/2.  # potential 

h=L/float(ng) # grid spacing
qh=0.5/(h*h)

kin=np.zeros((ng,ng)) # kinetic energy matrix
for i in range(ng):
    kin[i,i]=2.*qh
    if i>0: 
        kin[i,i-1]=-qh
        kin[i-1,i]=-qh
    if i<ng-1:
        kin[i+1,i]=-qh
        kin[i,i+1]=-qh

kin=-qh*fd.method(ng,deriv=2) # overwrite higher order schemes

pot=np.zeros((ng,ng)) # potential energy matrix
for i in range(ng): pot[i,i]=potential(-L/2+float(i)*h)

(val,vec)=la.eig(kin+pot)      # solve eigenproblem
print np.sort(val.real)[:10]   # show results (real part of eigenvalues, sorted !)
