import numpy as np
from read_input import read_input
from math import *

class Potential:
    """values, singularites etc. of potentials """
    
    def __init__(self,definition):
        self.name=definition.split('(')[0]
        print definition
        pars=definition.split('(')[0].strip(')').split(',')
        if self.name=='coulomb': 
            if pars != '': self.z=float(pars[0])
            else: self.z=1.

        elif self.name == 'cout': # "truncated coulomb"
            self.rp=float(pars.split(',')[0])
            self.z=float(pars.split(',')[1])

        elif self.name == 'cou1d': # "1d coulomb"
            print 'pars',pars
            self.z=float(pars[0])
            self.a=float(pars[1])

        else: sys.exit('potential not defined "'+name+'"')

        # check consistency of parameters
        if self.z<=0: exit('charge must be positive')

    @classmethod
    def read(cls):
        """create by reading in standard form"""
        definition=read_input('zero','.potential',force=True,doc='potential definition string')
        return Potential(definition)

    def v(self,q):  
            p=np.zeros(len(q))
            for i in range(len(q)): p[i]=self.pot(q[i])
            return p
        
    def pot(self,q):
        if self.name=='V=0': return 0.
        elif self.name == 'coulomb': return -self.z/q
        elif self.name == 'cout':                # "truncated coulomb"
            if x<self.rp: return -self.z/self.rp
            else: return (-self.z/x)
        elif self.name == 'cou1d': return -self.z/sqrt(q*q+self.a) # "1d coulomb"
        else: sys.exit('potential not defined: '+self.name)

