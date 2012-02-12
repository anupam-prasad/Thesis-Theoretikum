#! /usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt

from qmoperator import *
from axis import *

langle=0
rp=20.
np=60
order=10
kappa=0.05


def ho_eigen(ax,nval=2):
    hamiltonian=Operator('kinetic',ax)+Operator('kappa/2 |q^2|',ax)
    print 'n   eigenvalue'
    for n in range(nval):
        wf=hamiltonian.eigenstate(n)
        print n,(wf | (hamiltonian*wf)) / (wf|wf)

def hy_eigen(ax,nval=2):
    hamiltonian=Operator('kinetic',ax)+Operator('-1 |1/q|',ax)
    for n in range(nval):
        wf=hamiltonian.eigenstate(n)
        print n,(wf | (hamiltonian*wf)).real / (wf|wf).real


def cout_eigen(ax,nval=2):
    hamiltonian=Operator('kinetic',ax)+Operator('+1 cout(20.,1.)',ax)
    for n in range(nval):
        wf=hamiltonian.eigenstate(n)
        print n,(wf | (hamiltonian*wf)).real / (wf|wf).real

print 'fem'
#hy_eigen(Axis('r',200,0.,200.,kind='fem',order=20),8)
cout_eigen(Axis('r',300,0.,300.,kind='fem',order=10),10)


#ax1=Axis('r',np-order,0.,rp,kind='fem',order=order)
#ax2=Axis('r',np-ax1.n+1,0.,rp,kind='fem',order=np-ax1.n,axpar=['laguerre',kappa])
#ax=ax1+ax2
#print 'ax1+ax2'
#print 'hydrogen'
#hy_eigen(ax)
#print 'truncated'
#cout_eigen(ax,8)

