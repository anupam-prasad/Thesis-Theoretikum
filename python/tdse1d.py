#! /usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt

from qmoperator import *
from axis import *

# debug
import basisfunction as bf
np.set_printoptions(precision=7,suppress=True,linewidth=132)

def ho_eigen(ax,nval=5):
    hamiltonian=Operator('kinetic',ax)+Operator('+1 |q^2|',ax)
    print 'n   eigenvalue'
    for n in range(nval):
        wf=hamiltonian.eigenstate(n)
        print n,(wf | (hamiltonian*wf)) / (wf|wf)

ax=Axis.read()
