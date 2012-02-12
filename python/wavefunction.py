#! /usr/bin/env python

import matplotlib
import sys
from numpy import *
from axis import *

class WaveFunction:
    """wave function as elements of the Hilbert space
    all linear operations are defined"""

    def __init__(self,ax=None,time=None,model=None):
        """allocate by axis and time or form after model"""
        if ax is  not None:
            self.ax=ax
            self.time=time
        elif model is not None:
            self.ax=model.ax
            self.time=model.time
        else: sys.exit('must specify model= or ax=')

        self.S=False     # indicates that there is no overlap matrix multiplied on function
        self.c=zeros(self.ax.n,'complex')

    def __str__(self): return str(self.c)
            
    def __mul__(self,a): r=deepcopy(self);r.c*=a; return r   # multiply by scalar from the right
    def __add__(self,a): r=deepcopy(self);r.c+=a.c; return r # add two wavefunctions 
    def __sub__(self,a): r=deepcopy(self);r.c-=a.c; return r # subtract two wave fucntions
    def __imul__(self,a): self.c*=a; return self
    def __iadd__(self,a): self.c+=a.c; return self
    def __abs__(self): return abs(self.c)
 
    def __or__(self,wf):
        """ write scalar product in the form wf1|wf2"""
        if wf.S and self.S: sys.exit('both wave functions have been multiplied by operator, apply inverse')
        if self.ax != wf.ax: sys.exit('scalar product only for wave functions with same axes')
        elif self.S is not wf.S: # one of the functions carries overlapt already
            return vdot(self.c,wf.c)
        else:
            return vdot(self.c,dot(wf.ax.overlap(),wf.c))  

    def s_inv(self):
        """replace wf with S^-1 wf"""
        if not self.S: exit('cannot apply S^-1 on bare wave function')
        self.c=np.linalg.solve(self.ax.overlap(),self.c)
        self.S=False
        return self

    def grid(self,grid=None,add=0):
        """return grid, weights, and values of the wave function
        add ... add points beyond the minimum quadrature points
        grid... external grid (not implemented yet)
        """

        # which form of the wave function?
        if self.S: self.s_inv()

        base=[]
        weig=[]
        vals=[]
        for e in self.ax.e:
            if grid==None: x,w=e.basis.quadrature(add=add)
            else: exit('no external grid yet')
            for i in range(len(x)):
                v=e.val(x[i])[0]
                base.append(x[i])
                if w is not None: weig.append(w[i])
                vals.append(dot(v,self.c[e.i0:e.i0+e.n]))
        return array(base),array(weig),array(vals)
        
    @classmethod
    def create(cls,axes,name,pars):
        """given axes, name and parameters, return a wave function"""

        def gaussian(x,pars): return exp(-((x-pars[1])/pars[0])**2)

        # possible shapes on different coordinates
        shape={'gauss_x': gaussian,
               'gauss_y': gaussian
               }

        wf=WaveFunction(axes)
        iax=0
        func=shape[name+'_'+wf.ax.name]
        for e in wf.ax.e:
            x,w=e.basis.quadrature(add=0)
            for i in range(len(x)):
                v=e.val(x[i])[0]
                wf.c[e.i0:e.i0+e.n]+=v*w[i]*func(x[i],pars)

        # normalize
        wf.S=True  # we have calculated c_i = <i|wf>
        wf.s_inv() # remove the leading S
        wf.c/=sqrt(wf|wf)

        return wf

if __name__  == "__main__":
    """standalone: tests"""

    # generate wave function
    ax=Axis('x',100,-4.,4.,order=4)
    wf=WaveFunction.create(ax,'gauss',[1.,0.])

    # convert to grid
    g,w,v=wf.grid()

    plt.plot(g,abs(v))
    plt.plot(g,abs(np.exp(-g*g))/(myPi/2.)**0.25,'.')
    plt.show()

