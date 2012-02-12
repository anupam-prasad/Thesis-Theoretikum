import sys
import numpy as np
import matplotlib.pyplot as plt

class MyFunction:

    # draw the functions on the interval
    def plot(self,x0=None,x1=None,n=2,pts=100,normalization=None): 
        """Plot all functions on element"""
        if x0 is None: x0=self.range()[0]
        if x1 is None: x1=self.range()[1]
        if hasattr(self,'n'): 
            if self.n is not None: n=self.n

        x=np.linspace(x0,x1,pts)
        f=np.zeros((len(x),n))
        for i in range(len(x)): f[i,:]=self.val(x[i],n)[0]

        if normalization == 'max=1':
            for i in range(n):
                f[:,i]=f[:,i]/abs(f[:,i]).max()
                if f[:,i].max()<-f[:,i].min(): f[:,i]=-f[:,i]


        for i in range(n): plt.plot(x,f[:,i])
