#! /usr/bin/env python
import numpy as np
import findiff as fd

# set up some grid (equidistant)
N=12
L=float(N)
x=np.linspace(-L/2.,L/2,N+1) # N+1 points in interval, including end points
O=8

# control numpy matrix printing
# suppress ... set near-zero numbers = 0
np.set_printoptions(precision=5,linewidth=120,suppress=True)

print fd.general(x,O)*120*7
