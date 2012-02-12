#!/usr/bin/env python
import sys 
import random 
import copy
import numpy as np
from math import *

def energy(S,J,B): return J*(sum(S[:N-1]*S[1:])+S[N-1]*S[0])-B*sum(S)

def IsingMetropolis(N,J,B,kT,M,printinterval=1):
    """Metropolis algorithm for the Ising model"""
    random.seed(1)

    # start from "hot" state
    S=np.ones((N),'int')
    E=energy(S,J,B)
    sumS=0.
    sumE=0
    expEkT=[exp((-4*J-2*B)/kT),
            1.,
            exp((-4*J+2*B)/kT),
            exp((    -2*B)/kT),
            1.,
            exp((    +2*B)/kT),
            exp((+4*J-2*B)/kT),
            1.,
            exp((+4*J+2*B)/kT)
            ]
    for n in range(M):
        i=int(N*random.random())

        if  (S[i]==S[(i+1)%N]) != (S[i]==S[(i-1)%N]): 
            Ex=E
            dE=4
        elif S[i]==S[(i+1)%N]: 
            Ex=E-4*J
            dE=7
        else: 
            Ex=E+4*J
            dE=1

        if S[i]==1: 
            Ex+=2*B
            dE-=1
        else: 
            Ex-=2*B
            dE+=1

        if Ex<=E:
            S[i]=-S[i]
            E=Ex
        elif random.random()<expEkT[dE]: 
            S[i]=-S[i]
            E=Ex

        sumS+=S[0]
        sumE+=E
        if n%printinterval==1: print sumS/n,sumE/n

    return sumS/M, sumE/M
        
def exact(N,J,B,kT,M):

    # loop through all states S
    sumS=0.
    sumE=0.
    sumN=0.
    S=np.zeros((N))
    for n in range(2**N):
        nn=n
        for i in range(N):
            S[i]=(2*(nn%2)-1)
            nn/=2
        E=energy(S,J,B)
        expEkT=exp(-E/kT)   
        sumS+=S[0]*expEkT
        sumE+=E*expEkT
        sumN+=expEkT

    return sumS/sumN, sumE/sumN

if __name__ == "__main__": 

    J=1.
    B=1.
    kT=1.
    try:
        N=int(sys.argv[1])
        M=int(sys.argv[2])
    except:
        N=10
        M=100000
        
    print '\n using N= '+str(N)+' spins and M= '+str(M)+' metropolis steps'
#    print 'lim N->infty',sinh(B/kT)/sqrt(cosh(B/kT)**2-(1-exp(-4./kT*J)))
    print '               <s>                      <E>'
    print '       exact',exact(N,J,B,kT,M)
    print '  metropolis',IsingMetropolis(N,J,B,kT,M)

    try:
        N=int(sys.argv[3])
        M=1000000
        print '\n using N= '+str(N)+' spins and M= '+str(M)+' metropolis steps'
        print IsingMetropolis(N,J,B,kT,M,M/10)
    except:
        pass
