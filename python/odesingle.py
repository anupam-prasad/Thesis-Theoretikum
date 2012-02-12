import sys 
import numpy as np
import scipy.linalg as linalg
import scipy.linalg.flapack as lapack
import mytimer
import copy
from math import *
from cmath import *
from orthopol import *
from my_constants import myPi

class Chebyshev:
    def __init__(self,eps=1.e-12,maxorder=15):
        """Chebyshev exponentiation to maximal order"""
        if not isinstance(maxorder,int): exit('non-integer order for Chebychev')
        self.omax=maxorder
        self.explicit=True
        self.eps=eps
        x,w=ChebyshevPolynomial().quadrature(maxorder+10)

        v=np.zeros((len(x),maxorder))
        for i in range(len(x)): v[i,:]=ChebyshevPolynomial().val(x[i],maxorder)[0]

        def bracket(cv,v,f):
            m=np.zeros((np.shape(cv)[1],np.shape(v)[1]))
            for i in range(np.shape(cv)[1]):
                for j in range(np.shape(v)[1]):
                    m[i,j]=m[i,j]+np.dot(cv[:,i]*f,v[:,j])
            return m

        # modify weight, scale, shift
        np.set_printoptions(precision=2,suppress=True,linewidth=132)
        self.a=bracket(v,v[:,:1],np.exp(x)*w)*2/myPi
        self.a[0]*=0.5 # first Chebyshef is normalized to pi

    def __str__(self):
        string='Chebyshev, max order= '+str(self.omax)
        return string
        
    def step(self,x0,x1,y0,deri):
        t1=np.zeros((1,))
        tn=y0
        yn=self.a[0]*tn
        for n in range(1,self.omax):
            t2=t1
            t1=tn
            tn=ChebyshevPolynomial().b(n)*(x1-x0)*deri(t1,0.)+ChebyshevPolynomial().c(n)*t2
            t2=self.a[n]*tn
            yn+=t2
            if abs(t2).max() < abs(yn).max()*self.eps: break
        return yn

class CrankNicolson:

    def __init__(self,deri,solve):
        """deri(y,t)      ...derivative
        solve(y,t0,t1) ...return ys: y=ys-0.5*(t1-t0)*deri(ys,0.5*(t1+t0))
        """
        self.order=2
        self.explicit=False
        self.deri=deri
        self.solve=solve
        if self.solve==None: exit('for Crank-Nicolson must specify solver')

    def __str__(self):
        string='Crank-Nicolson'
        return string

    def step(self,x0,x1,y0,deri):
        ym=y0+0.5*(x1-x0)*deri(y0,0.5*(x1+x0))
        return self.solve(ym,x0,x1)

class RungeKutta:
    """general class for Runge-Kutta methods"""
    def __init__(self,name='classical'):
        self.order=2
        self.explicit=False
    def __str__(self):
        string='Crank Nicolson, 2nd order'

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
        
        if not self.explicit: exit('only explict methods implemented')
        k=[]
        for j in range(self.stage):

            # next support vector kj
            k.append(copy.deepcopy(y0))
            for l in range(j):
                ha=self.a[j,l]*(x1-x0)
                if ha!=0: k[j]+=k[l]*ha

            tj=x0+self.c[j]*(x1-x0)
            k[j]=deri(k[j],tj)

       # y1 = y0 + h* sum[j] bj * kj
        y1=copy.deepcopy(y0)
        for j in range(self.stage):
            y1+=k[j]*((x1-x0)*self.b[j])
        return y1

