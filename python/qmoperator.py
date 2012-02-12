import sys
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import scipy
from scipy.integrate import odeint
import scipy.linalg as la
from wavefunction import *
import physical_system as phys 
from odesolver import OdeSolver
  
definitions={
    # cartesian axis
    'x.overlap': ' 1  | ',
    'x.kinetic': '-1/2m |dd',
    # polar coordinates
    'r.overlap': ' 1  | ',
    'r.kinetic': '1/2m 1/qdq|1/qdq +l(l+1)/2m |1/q^2| '}

class VoxelOperator:
    """operator on a "voxel", i.e. a product of elements
    name ... name of the operator: kinetic
    coor ... (tuple of) coordinate(s)
    func ... (tuple of) basis funtions on voxel
    given as a list of matrices (for now)
    """
    
    def __init__(self,define,basis):
        
        # split into numerical prefactor and operator parts
        term=define.split()
        self.m=[(phys.parameter(term[0]),basis.matrix(term[1]))]
        self.i0=basis.i0          # starting index in global basis
        self.i1=basis.i0+basis.n  # end index in global basis
        for n in range(1,len(term)/2):
            self.m.append((phys.parameter(term[2*n]),basis.matrix(term[2*n+1])))

    def __add__(self,vo):
        """extend list of voxel operators"""
        vv=deepcopy(self)
        for m in vo.m: vv.m.append(m)
        return vv

class Operator:
    """consists of a set of voxel operators
    """

    def __init__(self,name,ax):

        self.name=name
        self.ax=ax
        self.vo=[]
        definitions.setdefault(ax.name+'.'+name,name) # return unaltered operator name as default
        self.definitions=definitions[ax.name+'.'+name]
        for n in range (len(ax.e)):
            self.vo.append(VoxelOperator(self.definitions,ax.e[n]))
            
    def __add__(self,op):
        """add by adding voxel operators """
        self.name+='+ '+op.name
        definitions.setdefault(self.ax.name+'.'+op.name,op.name) # return unaltered operator name as default
        res=deepcopy(self)
        for n in range(len(res.vo)): res.vo[n]=res.vo[n]+op.vo[n]
        return res

    def __mul__(self,app):
        """(overaload *: apply operator on wave function"""
        if isinstance(app,WaveFunction): return self.wf(app)
        else: sys.exit('cannot multiply operator on this object: '+str(app))
       
    def time(self,t):
        """update the operator time to t"""
        self.t=t
        return self
        

    def wf(self,wf): 
        """apply operator to wave function"""
        if wf.S: sys.exit('wave function is already multiplied by operator, apply inverse overlap first')
        op_wf=WaveFunction(model=wf)
        for v in self.vo:
            for m in v.m: 
                op_wf.c[v.i0:v.i1]+=m[0](self.t)*dot(m[1],wf.c[v.i0:v.i1])
        op_wf.S=True
        return op_wf 

    def matrix(self,t=0.,overlap=False):
        """matrix corresponding to operator"""
        if not hasattr(self,'m') or self.t != t: # setup matrix
            self.t=t
            self.m=np.zeros((self.ax.n,self.ax.n))
            for v in self.vo:
                for m in v.m:
                    self.m[v.i0:v.i1,v.i0:v.i1]=self.m[v.i0:v.i1,v.i0:v.i1]+m[0](t)*m[1]
        return self.m

    def eigensolver(self,time=None):
        """n'th eigenstate for operator at given time (if time-dependent)"""
        self.time=time
        (val,vec)=la.eig(self.matrix(),self.ax.overlap())
        ls=np.argsort(val.real)
        return val[ls],vec[:,ls]

    def eigenstate(self,n=0,time=None):
        """n'th eigenstate for operator at given time (if time-dependent)"""
        val,vec=self.eigensolver(time)
        wf=WaveFunction(self.ax,self.time)
        ls=np.argsort(val.real)
        wf.c=vec[:,ls[n]] 
        wf.e=val[ls[n]]
        return wf
    
    def evolve(self,wf,tgrid,method='runge kutta classical',error=0.):  # time-evolution
        """returns a time-evolution to a wave function"""

        # advanced version: apply operator
        def func(wf,t): return (self.time(t)*wf).s_inv()*complex(0.,-1.)

        solve1=OdeSolver(func,method,error,order=20) # a maximal order for the Chebyshev solver
        
        print 'initial time',wf.t
        wft=[wf]
        for t in tgrid:
            wft.append(WaveFunction(model=wf))
            wft[-1]=solve1.step(wft[-2].t,t,wft[-2])
            wft[-1].t=t
        return wft


