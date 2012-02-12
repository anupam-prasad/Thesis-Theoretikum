# convenient timer functions
import sys
import time

# timer table and its maximal size
mytimer_max=1000
mytimer_table=[]

class MyTimer:
    def __init__(self): 
        self.time=0.
        self.calls=0
        self.started=True
        self.t0=None

    def __str__(self):
        if self.t0 is None: print 'timer not used'
        else: return ' time '+str(self.time)+' calls '+str(self.calls)+' name '+str(self.name)

    def start(self,name):
        """asign name etc. at first call
        start clock
        """
        # first call: initialize
        if self.started: 
            if self.t0 is None: self.name=name
            else: exit('timer "'+self.name+'" not stopped after previous start')
        self.calls+=1
        self.started=True
        self.t0=time.clock()

    def stop(self):
        """stop clock
        add elapsed time since "start()" """
        self.time=self.time+time.clock()-self.t0
        if not self.started: exit('timer stopped but not started: "'+self.name+'"')
        self.started=False

def table():
    """print table of all timers"""
    for t in mytimer_table:
        if t.t0 is not None: print t

def create(n):
    """add n new timers to table and return them
    timer instances must be created as "default variables" 
    in the argument list of the routine that uses them 
    (see usage() for an example) """
    for i in range(n): mytimer_table.append(MyTimer())
    if len(mytimer_table)>mytimer_max: exit('number of timers exceeds limit of '+str(mytimer_max))
    return mytimer_table[-n:] # return the newly created timers

def usage(tm=create(3)):
    """ MyTimer usage example:
    usage(...other variables,tm=create(3))
    as there are no static variables in python, 
    we misuse the "quasi-static" nature of the default arguments 
    to keep timers between subroutine calls"""
    import numpy as np
    import scipy.linalg as la
    N=400

    

    tm[2].start('both')

    tm[0].start('timer0')
    ev=la.eig(np.random.random((N,N)))[0]
    tm[0].stop()

    tm[1].start('timer1')
    ev=la.eig(np.random.random((2*N,2*N)))[0]
    tm[1].stop()

    tm[2].stop()

#usage()
#usage()
#table()


