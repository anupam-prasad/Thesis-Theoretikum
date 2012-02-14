#! /usr/bin/env python
import sys
import matplotlib
matplotlib.use('PS') 
import traceback
from copy import *
from basisfunction import *
from finiteelement import *
from read_input import *
from my_constants import *
import matplotlib.pyplot as plt


class Axis:
    def __init__(self,name,n,lb,ub,kind='fem',order=2,axpar=None):
        """Setup elements with given order, equidistant 
        name... r,z,etc. (determines axis boundary conditions)
        n ..... number of discretization coefficients
        lb..... lower boundary of the axis
        ub..... upper boundary of the axis
        kind... 'quadrature' for quadrature grid, 'fem' for finite elements 
        order.. order of the fem method
        axpar.. further axis parameters
        """    

        # boundary conditions, jacobian
        ax_prop ={'r':     (False,True,'q^2',  0.,myInf),
                  'rho':   (False, True,'q',   0.,myInf),
                  'phi':   (True,True,'1',     0.,2.*myPi),
                  'x':     (True, True,'1',-myInf,myInf),
                  'y':     (True, True,'1',-myInf,myInf),
                  'z':     (True, True,'1',-myInf,myInf),
                  'cos':   (False,False,'1',-1.,1.)}

        self.name=name
        self.lb=lb
        self.ub=ub

        # finite element style axis
        self.e=[]
        self.n=n # will be overwritten for fem type

        if   kind == 'monomial': self.e.append(Monomial(n,lb,ub,jc=ax_prop[name][2]))  
        elif kind == 'legendre': self.e.append(LegendreScaled(n,lb,ub,jc=ax_prop[name][2]))
        elif kind == 'laguerre':
            if axpar is None: exit('for LaguerreExpon basis need to supply parameter kappa') 
            self.e.append(LaguerreExpon(n,float(axpar[0]),jc=ax_prop[name][2]))
        elif kind == 'fd': self.e.append(FDGrid(n,lb,ub,jc=ax_prop[name][2],order=order))
        elif kind == 'trigon': self.e.append(CosSin(n,jc=ax_prop[name][2]))
        elif kind == 'fem':
            ne=n/(order-1)
            i0=0
            for i in range(ne):
                self.e.append(FiniteElement(lb+(ub-lb)*i/ne,lb+(ub-lb)*(i+1)/ne,i0,order,
                                            lb0=ax_prop[name][0] and i==0,
                                            rb0=ax_prop[name][1] and i==ne-1,
                                            jc =ax_prop[name][2],fempar=axpar))
                i0=i0+self.e[i].n-1
            self.n=self.e[-1].i0+self.e[-1].n

        elif kind != 'dum':
            traceback.print_exc(file=sys.stdout)
            sys.exit('axis kind not defined "'+kind+'"')
            
        self.e[0].i0=0

        # check
        if self.lb<ax_prop[name][3]: exit('lower limit of '+self.name+'-axis < '+str(ax_prop[name][3]))
        if self.ub>ax_prop[name][4]: exit('upper limit of '+self.name+'-axis > '+str(ax_prop[name][4]))
        
    def __str__(self,text=None):
        if text is not None:
            string = '\n'+text
        else:
            string = '\n'
        for el in self.e: string+= '\n'+str(el)
        return string

    def show(self,text):
        print text
        for el in self.e: print el
        
    def plot(self,normalization=None):
        """Add to plot all element functions on axis"""
        for e in self.e: e.plot(normalization=normalization)

    def type(self):
        """axis types: fem, legendre, laguerre etc."""
        if hasattr(self.e[0],'basis'): return 'fem'
        else: return self.name

    def __add__(self,ax1):
        """append axis to axis"""
        
#        if isinstance(self.e[-1].basis,LaguerreExpon): 
#            exit('cannot append to axis with basis LaguerreExpon')

        if self.type()!='fem' or ax1.type()!='fem': exit('append only fem to fem')

        # redo first axis with open boundary on the right
        ax=deepcopy(self)
        ax.e=[]
        for e in self.e:
            ax.e.append(FiniteElement(e.x0,e.x1,e.i0,np.shape(e.tr)[0],e.jc,e.lb(),False,e.par))

        # concatenate
        for e in ax1.e:
            nf=np.shape(e.tr)[0]
            ax.e.append(FiniteElement(ax.e[-1].x1,
                                      ax.e[-1].x1+e.x1-e.x0,
                                      ax.e[-1].i0+ax.e[-1].n-1,
                                      nf,e.jc,False,e.ub(),e.par))
        ax.n=ax.e[-1].i0+ax.e[-1].n
        if hasattr(ax,'o'): del ax.o # overlap matrix is not valid
        return ax

    def len(self): return self.e[-1].i0+self.e[-1].n

            
    def overlap(self):
        """return overlap matrix for axis (create if needed)"""
        if not hasattr(self,'o'): # setup matrix
            self.o=np.zeros((self.n,self.n))
            for n in range(len(self.e)):
                i0=self.e[n].i0
                i1=i0+self.e[n].n
                self.o[i0:i1,i0:i1]=self.o[i0:i1,i0:i1]+self.e[n].matrix('|')
        return self.o

    def overlap_inv(self):
        """return inverse of overlap matrix (create if needed)"""
        if not hasattr(self,'i'): # setup inverse
            self.i=np.linalg.inv(self.overlap())
        return self.i

    def order(self):
        o=0
        for e in self.e: o=max(o,np.shape(e.tr)[0])
        return o

    @classmethod
    def read(cls,l=0):
        """read axis parameters from file as below and set up, and return axis
        .coordinate axis
        str name,int n,float lb,float up,str type,int order
        """
        name=read_input('none','.coordinate axis',1,l+1,'axis name: r,x,z,...',True)
        n=read_input(0,'.coordinate axis',2,l+1,'numer of discretization coefficients',True)
        lb=read_input(0.,'.coordinate axis',3,l+1,'lower boundary of axis',True)
        ub=read_input(0.,'.coordinate axis',4,l+1,'upper boundary of axis',True)
        typ=read_input('fem(legendre)','.coordinate axis',5,l+1,'type: fem,legendre...',False)
        order=read_input(4,'.coordinate axis',6,l+1,'finite element order',False)

        kappa=typ.rpartition('(')[2].partition(')')[0]
        typ=typ.partition('(')[0]
        typ=typ.strip('(')

        ax=Axis(name,n,lb,ub,typ,order,axpar=[kappa])
        return ax

    #Added Thursday 02.02.2012 - Calculates inner product of two vectors in FEM
    #representation
    def FEM_InnerProduct(self,v1,v2):
	vconj1=v1.conjugate()
	tempv=np.dot(vconj1,self.overlap())
	result=np.dot(tempv,v2)
	return result

    def FEM_Normalize(self,v):
	normsq=np.abs(self.FEM_InnerProduct(v,v))
	norm=np.sqrt(normsq)
	vout=v/norm
	return vout

    #Added Tuesday 31.01.2012 by Anupam - Returns vector corresponding to fem representation of
    #a momentum eigenstate for given axis and a specified momentum eigenstate - Have to change
    def FEM_MomentumEigenstate(self,momentum):

        nelem=self.n

        #Convert momentum eigenstate into finite element representation
        fem_vec=np.zeros(int(nelem))+0j
        invert_vec=np.zeros(int(nelem))+0j
        iter1=0

	
        for e in self.e:
                (x,w) = e.quadrature(add=10)
                (temp_a,temp_b)=e.val(x[0])
                iter2=len(temp_a)
                quad_iter=len(x)
                v=np.zeros([quad_iter,iter2])+0j
                
                for k in range(0,quad_iter):
                        (a,b) = e.val(x[k])+0j
                        c=a * np.exp(-momentum*x[k]*1j)
                        for l in range(0,iter2):
                                v[k][l]=c[l]*w[k]

                integral=v.sum(axis=0)
                for k in range(0,iter2):
                        invert_vec[iter1+k]=invert_vec[iter1+k]+integral[k]
                iter1=iter1+iter2-1
	
	fem_vec=np.dot(invert_vec,self.overlap_inv())
	return fem_vec


if __name__ == "__main__": 
    ax=Axis('x',4,0.,4.,'fem',3)+Axis('x',2,0.,4.,'fem',2)+Axis('x',1,0.,1.,'fem',2,axpar=['laguerre',1.])
    ax=Axis('x',9,0.,8.,'fem',3)+Axis('x',1,0.,1.,'fem',2,axpar=['laguerre',1.])
    ax=Axis('x',13,0.,12.,'fem',3)
    plt.subplot(311)
    ax.plot(normalization="max=1")
    for e in ax.e:
        plt.plot([e.x0,e.x0],[-0.05,1.05])
    plt.show()


