#! /usr/bin/env python
import sys
import numpy as np
from basisfunction import *
from recursive_integral import *

def coulomb_radial_integrand(rr,par):
    """ evaluate 
    bas0[ki](r1)*bas0[kj](r1)[min(r1,r2)^l/max(r1,r2)^(l+1)] bas1[mi](r2)*bas1[mj](r2) 
    for ki<=kj<=K, mi<=mj<=M, l=0,...,lmax
    """
    # get the basis function values at r1 and r2
    vals=[par[0].val(rr[0])[0],par[1].val(rr[1])[0]]

    # get all products of function values
    vij=[[],[]]
    for n in range(2):
        w=vals[n]
        for i in range(len(w)):
            for j in range(i+1):
                vij[n].append(w[i]*w[j])

    integs=np.zeros((len(vij[0])*len(vij[1])*(par[2]+1)))
    if min(rr)==0: return integs

    q1=min(rr)/max(rr)
    rrl=rr[0]**2*rr[1]**2/(q1*max(rr))

    ijl=-1
    for l in range(par[2]+1):
        rrl*=q1 # get next [r_<]^l/[r_>]^(l+1)

#......................................
# debug stuff...
#        # harmonic oscillator
#        if   l==0: rrl= 3*sum(rr**2)
#        elif l==1: rrl=-100*rr[0]*rr[1]
#        else: rrl=0.

#        if l==0: rrl=1.
#        else:    rrl=0.

        #zero
#        rrl=0
#        rrl*=(rr[0]*rr[1])**2
#......................................

        for ij0 in range(len(vij[0])):
            vij0_rrl=vij[0][ij0]*rrl
            for ij1 in range(len(vij[1])):
                ijl+=1
                integs[ijl]=vij0_rrl*vij[1][ij1]
    return integs

def coulomb_integrals_radial(bas0,bas1,lmax,acc=1.e-6,nquad=6):
    """coulomb integrals for the radial basis functions
    """
    vol=np.array([[bas0.x0,bas0.x1],[bas1.x0,bas1.x1]])

    if bas0.overlap(bas1): 
        # limited accuracies on the diagonal
        acc0=1.e-3
        nqu0=3
    else:
        acc0=acc
        nqu0=nquad
        
    val=recursive_integral(vol,coulomb_radial_integrand,params=[bas0,bas1,lmax],
                           acc_rel=acc0,acc_abs=acc0,nquad=nqu0)

    # return reshaped array
    return val.reshape((lmax+1,len(val)/(lmax+1)))

def test():
    bas0=LegendreScaled(3,4.,6.,jc='q^2')
    bas1=LegendreScaled(3,4.,6.,jc='q^2')
    ints=coulomb_integrals_radial(bas0,bas1,6,acc=1.e-3,nquad=3)
    print ints
    print len(ints),'integrals'

test()
