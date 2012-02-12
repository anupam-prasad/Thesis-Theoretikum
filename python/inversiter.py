#! /usr/bin/env python

# inverse iteration 

import sys
import numpy as np
from math import * 
import scipy.linalg as la
import mytimer as mt

# get timer
t=mt.create(10)

eps_machine=1.e-14 # just guessing... 

# numpy print style
np.set_printoptions(precision=12,suppress=True,linewidth=132)

# parameters
M=1
Klu=20
Kpower=1

# input 
try: N=int(sys.argv[1])
except: N=100
try: eguess=float(sys.argv[2])
except: eguess=5.

if M>N: exit('M>N')

# "reproducible" random series
np.random.seed(1)

# generate random symmetric matrix
A=np.matrix(np.random.random((N,N)))
for n in range(N): A[n,n:]=A[n:,n].T # symmetrize

# get the exact eigenvalues with highest modulus
t[0].start('la.eig')
ev0=la.eig(A)[0].real
t[0].stop()

# find nearest value
def nearest(array,value): return array[(np.abs(array-value)).argmin()]

np.set_printoptions(precision=12,suppress=True,linewidth=132)

# generate random initial vectors
kplt=0 # plot counter
eplt=np.zeros((M,Klu*Kpower))
vk=np.matrix(np.random.random((N,M)))
for L in range(Klu):
    # get the LU decompostion of A-shift
    t[5].start('inverse iteration')
    t[2].start('LU factor')
    lu,piv=la.lu_factor(A-eguess*np.eye(N))
    t[2].stop()

    for K in range(Kpower):
        t[3].start('LU solve')
        vk=np.matrix(la.lu_solve((lu,piv),vk,overwrite_b=True))
        # Schmidt orthonormalize
        for m in range(M): 
            vk[:,m]=vk[:,m]-vk[:,:m]*(vk[:,:m].T*vk[:,m])
            vk[:,m]=vk[:,m]/np.sqrt(vk[:,m].T*vk[:,m]) 
        t[3].stop()        

        t[4].start('A*vk')
        ek=(vk.T*A*vk).diagonal()
        t[4].stop()
        
    eguess=ek[0,0]
    t[5].stop()

    print ' appr.,exact',L,eguess,ev0[(np.abs(ev0-eguess)).argmin()]
print ' ---------------'
mt.table()

xlabel('Number of iterations')
ylabel('|DE/E|')
title('About as simple as it gets, folks')
axes().set_aspect('equal', 'datalim')
for m in range(M):
    plot(range(kplt), eplt[:,m], '-', lw=2)

show()

def inversiter(eguess,wfguess,A,S):

    vk=np.matrix(wfguess)
    olu,opi=la.lu_factor(Ovr)
    
    for L in range(Klu):
        # get the LU decompostion of A-e*S
        lu,piv=la.lu_factor(A-eguess*S)
        
        for K in range(Kpower):
            vk=np.matrix(la.lu_solve((lu,piv),vk,overwrite_b=True))
            ek=(vk.T*A*vk).diagonal()/(vk.T*S*vk)
        
        if (eguess-ek[0,0])/eguess<1.e-9: return ek[0,0],wk[:,0]
        eguess=ek[0,0]
        vk=S*vk
    
