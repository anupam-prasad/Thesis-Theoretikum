#! /usr/bin/env python

# an L=0 interparticle He code

"""He in interparticle coordinates"""

import sys

from copy import *
from numpy import *
import matplotlib.pyplot as plt
import scipy.linalg as la

from math import *
from cmath import *
from read_input import *

from matplotlib import colors, ticker
from mytimer import *
from my_constants import He_ground

tm=create(20)

# inputs
read_input_open(sys.argv[1])
chrg=read_input(2.,  ".nuclear charge",doc="nuclear charge")
xsym=read_input(1,   ".exchange symmetry",doc="exchange symmetry = +1 or -1")
alfa=read_input(1.,  ".interparticle basis",1,doc="r1-exponent",force=True)
beta=read_input(alfa,".interparticle basis",2,doc="r2-exponent")
psum=read_input(0,   ".interparticle basis",3,doc="maximal sum of powers",force=True)
mmax=read_input(psum,".interparticle basis",4,doc="maximal r1 power")
nmax=read_input(psum,".interparticle basis",5,doc="maximal r2 power")
kmax=read_input(psum,".interparticle basis",6,doc="maximal r3 power")
thet=read_input(0.,".complex scaling angle",doc="complex scaling angle in rad")
read_input_finish()

if chrg <=0: exit("need positive nuclear charge, is: "+str(chrg))

integral_table=[]

def factorial(n):
    f=1
    for i in range(n): f*=i+1
    return f


def integrals_extend(tab,mmax,nmax,kmax,a,b):
    """extend existing integral tables J,A to maximal indices mmax,nmax,kmax"""

    tm[10].start('extend tables')
    def extend_part(J,A,m0,m1,n0,n1,k0,k1):
        """partial extend for index ranges - J,A must be filled up to m0,n0,k0"""
        for m in range(m0-1,m1-1):
            for n in range(n0-1,n1-1):
                for k in range(k0-1,k1-1):
                    A[m+1,n+1,k+1]=(A[m,n+1,k+1]*m+
                                    A[m+1,n,k+1]*n+
                                    A[m+1,n+1,k]*k*2)
                    if m==0: A[m+1,n+1,k+1]+=2*(b**(-(n+k+1)))*factorial(n+k)
                    if n==0: A[m+1,n+1,k+1]+=2*(a**(-(m+k+1)))*factorial(m+k)
                    A[m+1,n+1,k+1]/=a+b
                    
                    J[m+1,n+1,k+1]=(A[m+1,n+1,k+1]+
                                    J[m,n+1,k+1]*m+
                                    J[m+1,n,k+1]*n)/(a+b)

    # new size arrays
    J=zeros((mmax+1,nmax+1,kmax+1))
    A=zeros((mmax+1,nmax+1,kmax+1))
    l=shape(tab[0])
    m=shape(A)
    J[:l[0],:l[1],:l[2]],A[:l[0],:l[1],:l[2]]=tab

    # fill in missing parts
    extend_part(J,A,   1,l[0],1,   l[1],l[2],m[2])
    extend_part(J,A,   1,l[0],l[1],m[1],   1,m[2])
    extend_part(J,A,l[0],m[0],   1,m[1],   1,m[2])

    tm[10].stop()
    return (J,A)
        
class InterBas:
    """monomial inter-particle coordinate (r1,r2,r3) basis functions (Hylleras)"""
    def __init__(self,f,m,n,k,a,b,c=0.):
        self.f=f
        self.p=array([m,n,k])
        self.e=array([a,b,c])

    def __str__(self): return str(self.f)+" "+str(self.p)+" "+str(self.e)
    def __ne__(self,o): return self.f!=o.f or  any(self.p!=o.p) or  any(self.e!=o.e)
    def __eq__(self,o): return self.f==o.f and all(self.p==o.p) and all(self.e==o.e)

    def integral(self,b):
        """get integral from table, extend table if needed"""
        tm[12].start('all integrals')

        # NOTE: the operations below are ridiculously costly
        #       seems to be a serious shortcoming of the implementation of the "objects"
        tm[6].start('powers and exponents')
        p=self.p+b.p
        e=self.e+b.e
        tm[6].stop()

        # select table
        tm[13].start('select table entry')
        s=(max(self.s,b.s)*(max(self.s,b.s)+1))/2+min(self.s,b.s)
        tm[13].stop()
        tm[8].start('get J from table')
        J=integral_table[s][0]
        tm[8].stop()

        # check whether table must be extended
        # NOTE: these checks are surprisingly costly
        tm[9].start('check shape')
        extab=shape(J)[0]<p[0]+3 or shape(J)[1]<p[1]+3 or shape(J)[2]<p[2]+3
        tm[9].stop()
        if extab: 
            integral_table[s]=integrals_extend(integral_table[s],
                             max(shape(J)[0],p[0]+3),
                             max(shape(J)[1],p[1]+3),
                             max(shape(J)[2],p[2]+3),e[0],e[1])
            J=integral_table[s][0]

        # get integral from table
        tm[5].start('get from table')
        val=J[p[0]+2,p[1]+2,p[2]+2]*self.f*b.f
        tm[5].stop()
        tm[12].stop()
        return val

    def __or__(self,b):
        """overlap between two basis functions"""
        return self.integral(b)

    def x12(self):
        """return basis with r1 and r2 exchanged"""
        b=InterBas(self.f,self.p[1],self.p[0],self.p[2],self.e[1],self.e[0],self.e[2])
        if hasattr(self,'s'): b.s=self.x
        if hasattr(self,'x'): b.x=self.s
        return b

class InterTerm:
    """a single term of a interparticle operator"""
    def __init__(self,def_term):
        """factors, powers, derivatives"""
        self.f=def_term[0]
        self.p=array(def_term[1:4])
        self.d=array(def_term[4:7])

    def __ne__(self,o): return self.f!=o.f or  any(self.p!=o.p) or  any(self.d!=o.d)
    def __eq__(self,o): return self.f==o.f and all(self.p==o.p) and all(self.d==o.d)
    def __str__(self): return str(self.f)+' '+str(self.p)+' '+str(self.d)

    def apply(self,bas):
        """applying to single basis function results in list of basis functions"""

        t=[]
        t.append(deepcopy(bas))
 
        # derivatives
        for l in range(3):
            for i in range(self.d[l]):
                for n in range(len(t)): 
                    t.append(deepcopy(t[n]))
                    # d/dr exp(-f r)   
                    t[-1].f*=-t[-1].e[l] 
                    # d/dr r^p       
                    t[ n].f*= t[ n].p[l]
                    t[ n].p[l]-=1

        # power changes and factors
        for tt in t: 
            tt.p+=self.p
            tt.f*=self.f
        # remove zero terms
        tn=[]
        for tt in t : 
            if tt.f!=0.: tn.append(tt)
        return tn

    # merge same basis functions terms
    # (not implemented)

    def sym(self):
        """exchange symmetric operator term"""
        return InterTerm([self.f,self.p[1],self.p[0],self.p[2],self.d[1],self.d[0],self.d[2]])

class InterOper:
    """interparticle coordinate (r1,r2,r3) operator is a set of terms, consisting of 
    factors, power changes, and derivatives for each coordinate"""
    def __init__(self,op,sym=True):
        """converts a list 'op' of [factors,values,derivatives] into an operator
        if sym==True, r1<->r2 terms will be supplemented to the operator
        """
        self.name=op[0]
        self.term=[]
        # convert all list entries to InterTerm type, supplement by exchanged term where needed
        for o in op[1:]: 
            self.term.append(InterTerm(o))
            if sym and InterTerm(o).sym()!=InterTerm(o): self.term.append(InterTerm(o).sym())

    def __str__(self): 
        s=self.name
        for t in self.term: s+='\n '+str(t)
        return s

    def apply(op,basi):
        """apply an operator term to a single basis function
        returns a list of basis functions
        """
        opt=[]
        # apply each term
        for t in op.term: 
            ta=t.apply(basi)
            for a in ta: opt.append(a)
        return opt

    def matrix(self,bas,xsym):
        """construct a matrix from an operator for a given basis
        bas  ... list of InterBas
        xsym ... exchange symmetry 0, +1 or -1
        """
        Omat=zeros((len(bas),len(bas)))
        basx=[]
        asym=[]
        # loop through left hand basis functions
        for i in range(len(bas)):

            # list of 1-2 exchanged basis functions (if needed)
            if xsym!=0: basx.append(bas[i].x12())
            # indicate whether exchanded function term should be added
            asym.append(basx[i]!=bas[i] and xsym !=0)

            # generate a list of terms by applying the oberator to the basis function
            tm[2].start('apply operators')
            Obasi=self.apply(bas[i])
            tm[2].stop()

            # loop through right hand basis functions
            tm[3].start('compute integrals')
            for j in range(i+1):
                for term in Obasi:
                    if asym[j]:
                        Omat[i,j]+=(term|bas[j])+xsym*(term|basx[j])
                    else:
                        Omat[i,j]+=(term|bas[j])*2

                # assume symmetric matrix
                Omat[j,i]=Omat[i,j] 
            tm[3].stop()

        return Omat
                        
def basis_set(psum,mmax,nmax,kmax,alfa,beta,xsym):
    """ generate basis
    r1^m r2^n r3^k exp(-alfa r1 - beta r2)
    with the constraint on the powers
    m+n+k <= psum,   m <= mmax, n <= nmax, k <= kmax
    if xsym=+1 or -1
    include only functions that are not related to previous by exchange symmetry
    """
    tm[0].start('basis setup')
    bas=[]
    for m in range (mmax+1):
        for n in range (nmax+1):
            for k in range (kmax+1):
                if m+n+k>psum: break
                bn=InterBas(1.,m,n,k,alfa,beta)

                # add only if not related by symmetry
                app=True
                bs=bn.x12()
                for b in bas: 
                    if b==bs: app=xsym==0; break

                # accept basis function (except exchange symmetric for anti-symmetric case)
                if app and (bs!=bn or xsym!=-1): bas.append(bn)

    # assign numbers to sets of exponents
    nset=0                
    for i in range(len(bas)):
        for a in bas[:i]:
            if all(a.e==bas[i].e): # exponents match
                bas[i].s=a.s       # use number of matching set
                break
        if not hasattr(bas[i],'s'):
            bas[i].s=nset # new set
            nset+=1       # increment set counter

    # exchanged exponent numbering
    for i in range(len(bas)):
        for a in bas[:i]:
            if a.e[1]==bas[i].e[0] and a.e[0]==bas[i].e[1]: 
                # exchanged exponents match
                bas[i].x=a.s       # use number of matching set
                break
            if hasattr(a,'x') and all(bas[i].e==a.e): 
                    bas[i].x=a.x
                    break
        if not hasattr(bas[i],'x'):
            bas[i].x=nset # new set
            nset+=1       # increment set counter
        
    tm[0].stop()

    # create emtpy integral tables
    nset*=2
    for i in range((nset*(nset+1))/2+1): integral_table.append((zeros((1,1,1)),zeros((1,1,1))))       
    return bas

def rescale(A,B):
    d=zeros((shape(B)[0]))
    for i in range(len(d)):
        d[i]=B[i,i]
        for j in range(i+1):
            A[i,j]=A[i,j]/sqrt(d[i]*d[j])
            A[j,i]=A[i,j]

def test():

    # hamiltonian (unsymmetrized)
    # factor, power-changes, derivatives 
    # (exchange-symmetric part will be supplemented)
    ham=['hamiltonian',
         [-2.,   0, 0,-1,0,0,1],
         [-1.,   0, 0, 0,0,0,2],
         [-1.,  -1, 0, 0,1,0,0],
         [-0.5,  0, 0, 0,2,0,0],
         [-0.5,  1, 0,-1,1,0,1],
         [ 0.5, -1, 2,-1,1,0,1],
         [-0.5, -1, 0, 1,1,0,1],
         [-chrg,-1, 0, 0,0,0,0],
         [ 1. ,  0, 0,-1,0,0,0],
         ]
    kin=['kinetic energy',
         [-2.,   0, 0,-1,0,0,1],
         [-1.,   0, 0, 0,0,0,2],
         [-1.,  -1, 0, 0,1,0,0],
         [-0.5,  0, 0, 0,2,0,0],
         [-0.5,  1, 0,-1,1,0,1],
         [ 0.5, -1, 2,-1,1,0,1],
         [-0.5, -1, 0, 1,1,0,1],
         ]
    cou=['coulomb potentials',
         [-chrg,-1, 0, 0,0,0,0],
         [ 1. ,  0, 0,-1,0,0,0],
         ]

    # overlap
    ovr=['overlap',[1.,0,0,0,0,0,0]]
    

    # set up the basis
    bas=basis_set(psum,mmax,nmax,kmax,alfa,beta,xsym)

    # get operator matrices
    tm[1].start('operator setup and matrix formation')
    Ovr=InterOper(ovr)
    Ham=InterOper(ham)
    Kin=InterOper(kin)
    Pot=InterOper(cou)
    Hmat=exp(-2j*thet)*Kin.matrix(bas,xsym)+exp(-1j*thet)*Pot.matrix(bas,xsym)
    Omat=Ovr.matrix(bas,xsym)
    tm[1].stop()

    # rescale to get nicer condition numbers
    rescale(Hmat,Omat)
    rescale(Omat,Omat)

    # diagonalize
    set_printoptions(precision=9,suppress=True,linewidth=132)
    print '\n eigenvalues'
    evals=la.eigvals(Hmat,Omat)
    ls=argsort(evals.real)
    print evals[ls[:20]]
    val=sort(la.eigvals(Omat.real)).real
    print 'He ground state energy', He_ground,'(au), best known value'
    print 'Overlap condition number',val[-1]/val[0]
    table()

test()
