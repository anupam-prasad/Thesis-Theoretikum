#! /usr/bin/env python

# explore the Krylov subspace

import sys
import numpy as np
from math import * 
import scipy.linalg as la
from matplotlib.pyplot import *

eps_machine=1.e-14 # just guessing... 

# numpy print style
np.set_printoptions(precision=12,suppress=True,linewidth=132)

# parameters
M=3

#----------------------------------------------------
# input: if no command line argument is given use default
try: N=int(sys.argv[1])
except: N=100
try: shift=float(sys.argv[2])
except: shift=0.
try: Kmax=int(sys.argv[3])
except: Kmax=100
#-----------------------------------------------------

# check input
if M>N: exit('M>N')
exact=np.zeros((M))
evals=np.zeros((Kmax,M))

# reproducible "random" series
np.random.seed(1)

# generate random symmetric matrix
A=np.matrix(np.random.random((N,N)))
for n in range(N): A[n,n:]=A[n:,n].T # symmetrize

# get the exact eigenvalues with highest modulus
ev0=la.eig(A-shift*np.eye(N))[0].real
ls=np.argsort(abs(ev0))
exact[max(-M,-N):]=ev0[ls[max(-M,-N):]]+shift


# generate random initial vectors
dplt=np.zeros((Kmax,M))
rplt=np.zeros((Kmax-1,M))
vk=np.matrix(np.random.random((N,M)))
for K in range(Kmax):

    # Schmidt orthonormalize
    for m in range(M): 
        vk[:,m]=vk[:,m]-vk[:,:m]*(vk[:,:m].T*vk[:,m])
        vk[:,m]=vk[:,m]/np.sqrt(vk[:,m].T*vk[:,m])        

    vk0=vk
    vk=A*vk
    ek=(vk.T*vk0).diagonal()
    vk=vk-shift*vk0
    dplt[K,:]=abs(np.sort(ek.real)/np.sort(exact)-1.)

np.set_printoptions(precision=12,suppress=True,linewidth=132)
print 'appr.',np.sort(ek.real)
print 'exact ',np.sort(exact)
print dplt[Kmax-1,:]
print dplt[Kmax-2,:]
print dplt[Kmax-2,:-1]/dplt[Kmax-1,:-1]
xlabel('Number of iterations')
ylabel('|DE/E|')
title('Convergence of the power method')
axes().set_aspect('equal', 'datalim')
axes().set_yscale('log')
for m in range(M-1): plot(range(Kmax), dplt[:,m], '-', lw=1)
#for m in range(M): plot(range(Kmax-100), rplt[99:,m], '-', lw=1)

show()

var = raw_input("Hit any key to quit... ")

