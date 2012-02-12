#! /usr/bin/env python

# demonstration of the FFT split step method

import sys
from read_input import read_input_open,read_input
import axis
from qmoperator import *
import physical_system as phys

# get the input list from input file
read_input_open(sys.argv[1])

tmin=read_input(0.,'.time grid',1,doc='start time for time-propagation')
tmax=read_input(1.,'.time grid',2,doc='end time for time-propagation')
tinc=read_input(0.1,'.time grid',3,doc='time-increment')

# set up the axis
ax=Axis.read()

kinetic=Operator('kinetic',ax)
potential=Operator('+1 |q^2|',ax)
hamiltonian=kinetic+potential

wf0=WaveFunction.create(ax,'gauss',[1.,1.])

# time-propagate and plot
wf0.t=tmin
t=np.arange(tmin+tinc,tmax+tinc,tinc)
wft=hamiltonian.evolve(wf0,t)

print 'initial and final norm',abs(wft[0]|wft[0]),abs(wft[-1]|wft[-1])

fig = plt.figure()
pax = fig.add_subplot(111)

def animate():
    g,w,v=wft[0].grid(add=20)
    line,=pax.plot(g,v)

    for wf in wft:
        g,w,v=wf.grid(add=20)
        line.set_ydata(abs(v)**2)
        fig.canvas.draw()
    
win = fig.canvas.manager.window
fig.canvas.manager.window.after(10, animate)
plt.show()

