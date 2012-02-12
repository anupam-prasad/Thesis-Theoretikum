# a 2-dimensional time-dependent Schroedinger solver

import sys

from numpy import *
import matplotlib.pyplot as plt

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la

from axis import *
from tensor import *
from math import *

from matplotlib import colors, ticker

from scipy.integrate import odeint

# get the input list from input file
read_input_open(sys.argv[1])

# set up the axes
ax=axis.read()

# set up the laser field

# get an initial vector
wf=hamiltonian.eigenstate(time=tstart)

# time-propagate and plot
nt=int((tmax-tmin)/tinc+0.1)+1
for t in range(nt): 
    wf=hamiltonian.evolve(wf,wf.time+(tmax-tmin)/long(nt))
