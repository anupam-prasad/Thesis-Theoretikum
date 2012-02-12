#! /usr/bin/env python

# a 2-dimensional time-dependent Schroedinger solver

import sys
from read_input import read_input_open,read_input
import axis
from qmoperator import *
import physical_system as phys

# get the input list from input file
read_input_open(sys.argv[1])

# set up the axis
ax=Axis.read()

kinetic=Operator('kinetic',ax)
potential=Operator('+1 |q^2|',ax)
interaction=Operator('laser(t) |q|',ax)
hamiltonian0=kinetic+potential
hamiltonian=hamiltonian0+interaction

# read a time grid
tmin=read_input(0.,'.time grid',1,doc='start time for time-propagation')
tmax=read_input(1.,'.time grid',2,doc='end time for time-propagation')
tinc=read_input(0.1,'.time grid',3,doc='time-increment')
error=read_input(0.,'.time grid',4,doc='time-increment')

# time-propagate and plot
wf=WaveFunction.create(ax,'gauss',[1.,1.])
wf.t=tmin

# set time-grid and prpagate
t=np.arange(tmin+tinc,tmax+tinc,tinc)
wf0=hamiltonian.evolve(wf,t,error=error)

# show conservation of norms
print 'initial and final norms, wf0',abs(wf0[0]|wf0[0]),abs(wf0[-1]|wf0[-1])

# animated output
def animate():
    g,w,v0=wf0[0].grid(add=20)
    line,=pax.plot(g,v0)

    for n in range(len(wf0)):
        g,w,v0=wf0[n].grid(add=20)
        line.set_ydata(abs(v0)**2)
        fig.canvas.draw()
    
fig = plt.figure()
pax = fig.add_subplot(111)
win = fig.canvas.manager.window
fig.canvas.manager.window.after(10, animate)
plt.show()

