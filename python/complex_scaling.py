#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from read_input import *
from axis import *
from tensor import *
from myeigen import inversiter

# open the first string after the program name as input file 
read_input_open(sys.argv[1])

# read system parameters
field=read_input(0.,'.system parameters',1,doc='field strength (a.u.)')
lmin =read_input(0, '.system parameters',2,doc='minimal angular momentum')
lmax =read_input(0, '.system parameters',3,doc='maximal angular momentum')
mqn  =read_input(0, '.system parameters',4,doc='m quantum number')
theta=read_input(0.,'.complex scaling angle',1,doc='how to find roots: full,invit')
steps=read_input(3, '.complex scaling angle',2,doc='repeat calculation steps times with different angles')

lmin=max(mqn,lmin) # minimal angular momentum cannot be below m-quantum number

ax=axis_read(0) # read the axis

# check input, write docu-file
read_input_finish()

# hydrogen atom (L=0)
print 'hydrogen atom'
hamr=np.zeros((ax.n,ax.n))
kinr=np.zeros((ax.n,ax.n))
potr=np.zeros((ax.n,ax.n))
angr=np.zeros((ax.n,ax.n))
dipr=np.zeros((ax.n,ax.n))
ovrr=np.zeros((ax.n,ax.n))
for e in ax.e:
    kinr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('1/qdq|1/qdq')
    potr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('coulomb(1.)')
    angr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('|1/q^2|')
    dipr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('|q|')
    ovrr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('|')
hamr=kinr+potr

ovro=np.eye((lmax+1-lmin))
ango=np.zeros((lmax+1-lmin,lmax+1-lmin))
dipo=np.zeros((lmax+1-lmin,lmax+1-lmin))
ango[0,0]=lmin*(lmin+1)
for l in range(lmin,lmax+1):
    ango[l-lmin,l-lmin]=l*(l+1)
    dipo[l-lmin-1,l-lmin  ]=sqrt((l*l-mqn*mqn)/float(4*l*l-1))
    dipo[l-lmin,  l-lmin-1]=dipo[l-lmin-1,l-lmin]

Ovr=Tensor(TensorProd(ovro,ovrr))

#np.set_printoptions(precision=9,suppress=True,linewidth=132)
#print np.sort(val.real)[:10]   # show results (real part of eigenvalues, sorted !)

plotsym=('<','>','x')
for s in range(steps):
    th=theta*(s+1)/float(steps)
    Tmat=exp(-2j*th)*(Tensor(TensorProd(ovro,kinr))+Tensor(TensorProd(ango,angr))).array()
    Vmat=exp(-1j*th)* Tensor(TensorProd(ovro,potr)).array()
    Hmat=Tmat+Vmat
    Dmat=exp( 1j*th)*Tensor(TensorProd(dipo,dipr)).array()
    Omat=Ovr.array()

    val=la.eigvals(Hmat+field*Dmat,Omat)
    
    ls=np.argsort(abs(val.imag))
    print val[ls[:10]]
    plt.plot(val[ls[:20]].real,val[ls[:20]].imag,plotsym[s])

plt.show()

