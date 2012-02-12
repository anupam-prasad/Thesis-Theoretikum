#!/usr/bin/env python
"""wrapper for general ode-solvers"""
import sys 
import numpy as np
import scipy.linalg as linalg
import scipy.linalg.flapack as lapack
import mytimer
import copy
import matplotlib.pyplot as plt
from math import *
from cmath import *


class RungeKutta:
    """general class for Runge-Kutta methods"""
    def __init__(self,name='classical'):
        self.order=2
        self.explicit=False
    def __str__(self):
        string='Crank Nicolson, 2nd order'

    def step(self,x0,x1,y0,deri,ideri):
        y1=y0+0.5*(x1-x0)*deri(x0,y0)

    # butcher tableaus for various methods
    global butcher_table
    butcher_table={
        'classical': (4,4,
                      [[0.,   0.,   0., 0.],
                       [1/2.,   0., 0., 0.],
                       [  0., 1/2., 0., 0.],
                       [  0.,   0., 1., 0.]],
                      [1/6.,1/3.,1/3.,1./6],
                      [0.,  1/2.,1/2.,1.],
                      ),
        'midpoint': (2,2,
                     [[0.,  0.],
                      [1/2.,0.]],
                     [0.,1.],
                     [0.,0.5]
                     ),
    'euler-forward': (1,1,
            [[0.,]],
            [1.],
            [1.]
            ),
    'euler-backward': (1,1,
            [[1.,]],
            [1.],
            [1.]
            )
        }

    def __init__(self,name='classical'):
        self.order=butcher_table[name][0]
        self.stage=butcher_table[name][1]
        self.a=np.array(butcher_table[name][2])
        self.b=np.array(butcher_table[name][3])
        self.c=np.array(butcher_table[name][4])
        self.explicit=True
        for j in range(np.shape(self.a)[1]):
            for i in range(j+1): 
                if self.a[i,j]!=0:
                    self.explicit=False
                    break
            if not self.explicit: break # no need to continue

    def __str__(self):
        string='Runge-Kutta: a,b,c'
        string+='\n'+str(self.a)
        string+='\n'+str(self.b)
        string+='\n'+str(self.c)
        return string

    def step(self,x0,x1,y0,deri,ideri=None):
        
        if self.explict:
            k=[]
            for j in range(self.stage):

                # next support vector kj
                k.append(copy.deepcopy(y0))
                for l in range(j):
                    ha=(x1-x0)*self.a[j,l]
                    if ha!=0: k[j]+=ha*k[l]

                tj=x0+(x1-x0)*self.c[j]
                k[j]=deri(k[j],tj)

        else:
            # treat only special cases here
            



        # y1 = y0 + h* sum[j] bj * kj
        y1=copy.deepcopy(y0)
        for j in range(self.stage):
            y1+=((x1-x0)*self.b[j])*k[j]
        return y1


class OdeSolver:
    """general solver for ODEs:
    y(t) = func(y,t)
    """

    def __init__(self,func,method='runge kutta classical',error=1.e-7,initial_step=None):
        """func ... derivative function
        method ... name of the method
        """
        self.func=func
        self.error=error
        self.hstep=initial_step
        
        self.nreject=0
        self.nsteps=0

        if   method=='runge kutta classical': self.method=RungeKutta('classical')
        elif method=='runge kutta midpoint':  self.method=RungeKutta('midpoint')
        elif method=='runge kutta euler-forward':  self.method=RungeKutta('euler-forward')
        elif method=='runge kutta euler-backward':  self.method=RungeKutta('euler-backward')
        else: exit('ODE method not implemented: '+method)

    def __str__(self):
        string=str(self.nreject)+' '+str(self.nsteps)
        return string

    def step(self,t0,t1,y0):
        """advance solution from t0 to t1"""
        
        if t0>t1: exit('time propatation only to larger times')

        # self.error=0 ... no error control
        if self.error==0: return self.method.step(t0,t1,y0,self.func)

        # guess step size if None given
        if self.hstep==None: self.hstep=(t1-t0)/64. 

        tm=t0
        ym=y0
        while tm<t1:
            self.nsteps=+1
            tn=min(tm+self.hstep,t1)

            # get step and half-step solutions over interval
            y1=self.method.step(tm,(tm+tn)/2,ym,self.func)
            y2=self.method.step((tm+tn)/2,tn,y1,self.func)
            y1=self.method.step(tm,       tn,ym,self.func)

            # get the error estimate and step size correction factor
            err_abs=abs(y1-y2).max()/(2.**(self.method.order+1)-2)
            err_abs=err_abs/max(self.error,abs(y1).max())
            correction_factor=(1/(self.error+1.e-12))*err_abs

             # accept step
            if correction_factor<1: 
                tm=tn
                ym=y2
            else:
                self.nreject+=1
                if self.nreject*10>self.nsteps and self.nsteps>100: 
                    print ' !!! many rejects - check your algorithm:',self.nreject,' out of',self.nstep 
           
            # estimate next step size
            self.hstep=min(self.hstep/((5.*correction_factor)**(1./self.method.order+1)),
                          1.5*self.hstep)
        # done
        return ym

    @classmethod
    def test(cls):

        print '\n\n *** d/dt y = -y [0,1] ***\n'
        try:
            eps=float(sys.argv[3])
            N=int(sys.argv[4])
        except:
            sys.exit('supply eps, N (command line arguments): accuracy, matrix dimension N')

        def func0(y,t): return -y
        expt=OdeSolver(func0,'runge kutta classical',eps)
        print 'method:\n',expt.method

        # set initial values, solve and compare to exact
        y=1.
        t1=1.
        y=expt.step(0.,t1,y)
        if eps>0.: print '\ny*exp(t1)=1?:',(y*exp(t1)-1).real,'\ntrue error/eps:',((y*exp(t1)-1)/eps).real
        else: print '\ny*exp(t1)=1?:',y*exp(t1),'\ntrue error:',(y*exp(t1)-1)

        print '\n\n *** time-propagation by a random symmetric matrix ***'

        # generate a random hermitian matrix
        np.random.seed(1)
        ham=np.matrix(np.random.random((N,N)))
        ham=(ham+ham.T)/2 

        # define d/dt y = -i H y
        def func1(y,t): return -1j*ham*y

        expm=OdeSolver(func1,'runge kutta classical',eps)

        # get a normalized random vector
        y0=np.matrix(np.random.random((N,1))+1j*np.random.random((N,1)))
        y0=y0/sqrt((y0.T*y0).real)

        # solve by full diagonalization 
        # y(t) = V exp(-i*t*eigvals) V.T y(0)
        val,vec=linalg.eig(ham)  
        vec=np.matrix(vec)
        cvy=vec.T*y0
        for i in range(len(cvy)): cvy[i]*=exp(-1j*val[i]*t1)
        yfd=vec*cvy
        print 'maximal eigenvalue: ',abs(val).max()

        # solve
        yrk=expm.step(0.,t1,y0)
        print 'solver info\n',expm

        # compare results
        np.set_printoptions(precision=1,suppress=True,linewidth=132)
        print '\nerror/eps'
        if eps>0.: print abs(yrk.T/yfd.T-1)/eps
        else: print 'errors\n ',abs(yrk.T/yfd.T-1)
 
    def consistency(self):

        y0=1.
        t1=1.
        tt=[]
        ee=[]
        print self.method
        print '\nt,y,y-exp(-t),log(|y-exp(-t)|)'
        for t in np.arange(t1/10,t1,t1/10.):
            y=self.step(0.,t,y0)
            print "%f" % t,(y-exp(-t)).real
            tt.append(t)
            ee.append(abs(y-exp(-t)))

        print '\nlog/log slope-1 ~ order: ',((log(abs(ee[-1]))-log(abs(ee[0])))/(log(abs(tt[-1]))-log(abs(tt[0])))-1).real
        print '\n'
        plt.loglog(tt,ee)
        plt.show()

    def convergence(self):
        y0=1.
        t1=1.
        tt=[]
        ee=[]
        
        for N in range(1,11):
            y=y0
            emax=0
            for t in np.arange(0,t1,t1/N):
                y=self.step(t,t+t1/N,y)
                emax=max(emax,abs((y*exp(t+t1/N))-1))
            tt.append(N)
            ee.append(emax)
            
        print '\nlog/log slope-1 ~ order: ',((log(abs(ee[-1]))-log(abs(ee[0])))/(log(abs(tt[-1]))-log(abs(tt[0])))+1).real
        print '\n'
        print ee
        plt.loglog(tt,ee)
        plt.show()

    def stability(self):
        lamax=3.
        N=50
        y0=1.
        nplt=10

        print 'hey'

        ll=[]
        em=[]
        ym=[]
        la=0.
        for p in range(nplt):
            la=la+lamax/nplt
            tt=[]
            ee=[]
            yy=[]
            y=y0
            for n in range(N):
                y=self.step(0,la,y)
                tt.append(float(n)/N)
                yy.append(abs(y))
                ee.append(abs(y*exp(la*(n+1)/N)-1))
            ll.append(la)
            ym.append(yy[-1])
            em.append(ee[-1])
            if p%(nplt/10)==0: 
                plt.plot(tt,np.log10(np.array(yy)))
        
        print '\n lambda  error  value'
        for i in range(len(ll)):
            print '%6.1f' % ll[i],'%7.1e' %em[i],'%7.1e' % ym[i]
        plt.show()


if __name__  == "__main__":
    """standalone: tests"""

    global func0
    def func0(y,t): return -y

    try:
        what=sys.argv[1]
        method=sys.argv[2]
    except:
        sys.exit('supply  what to do and method (command line arguments)')

    print '\n\n e(-1)-y_N for d/dt y = - y on [0.1]\n'

    if method=='classical': expt=OdeSolver(func0,'runge kutta classical',0.)
    elif method=='midpoint':  expt=OdeSolver(func0,'runge kutta midpoint',0.)
    elif method=='euler-forward':  expt=OdeSolver(func0,'runge kutta euler-forward',0.)
    elif method=='euler-backward':  expt=OdeSolver(func0,'runge kutta euler-backward',0.)
    else: print 'not set up for method:', method

    if what == 'consistency': expt.consistency()
    elif what == 'convergence': expt.convergence()
    elif what == 'control': OdeSolver.test()
    elif what == 'stability': expt.stability()
    else: exit('uknown input: '+what+'. Enter control, consistency, convergence, stability as first argument')
    
   
