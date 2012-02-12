#! /usr/bin/env python

"""He in independent particle coordinates"""

import sys

from copy import *
from numpy import *
from orthopol import *
from basisfunction import *
import matplotlib.pyplot as plt
import scipy.linalg as la

from math import exp
from read_input import *

from matplotlib import colors, ticker
from mytimer import *
from my_constants import *
from axis import Axis
from coulomb_integral import *

class GauntCoeffTable:
    def __init__(self,lmax):
        """set up tables of associated Legendre polynomials up to maximal l-values"""
        self.quadx,self.quadw=LegendrePolynomial().quadrature(max(2,(3*lmax+1)/2))
        self.Passoc=[]
        for m in range(lmax+1):
            self.Passoc.append([])
            pa=LegendreAssociated(lmax,m)
            for l in range(m): self.Passoc[m].append([]) # append empty lists for more convenient indexing 
            v=zeros((len(self.quadx),lmax+1))
            for i in range(len(self.quadx)):
                v[i,m:]=pa.val(self.quadx[i])[0]
            for l in range(m,lmax+1):
                # append normalized sets of values
                self.Passoc[m].append(array(v[:,l]*pa.normY(l,m)))

    def coeff(self,lc,l1,l2,mc,m1,m2):
        """return Gaunt coefficient"""
        if mc!=m1+m2: return 0
        if abs(l1-l2)>lc>l1+l2: return 0
        if abs(mc)>lc or abs(m1)>l1 or abs(m2)>l2: return 0
        pc=self.Passoc[abs(mc)][lc]
        p1=self.Passoc[abs(m1)][l1]
        p2=self.Passoc[abs(m2)][l2]
        return sum(pc*p1*p2*self.quadw)*2*myPi

    def test(self):
        """simple test with one of the L's==0"""
        for m in range(len(self.Passoc)):
            for l in range(m,len(self.Passoc[m])):
                if abs(self.coeff(l,l,0,m,m,0)-sqrt(0.25/myPi))>1.e-10:
                    exit('failed test1')
                if abs(self.coeff(l,l,0,-m,-m,0)-sqrt(0.25/myPi))>1.e-10:
                    exit('failed test2')
                if abs(self.coeff(0,l,l,0,m,-m)-sqrt(0.25/myPi))>1.e-10:
                    exit('failed test3')
                if abs(self.coeff(0,l,l,0,-m, m)-sqrt(0.25/myPi))>1.e-10:
                    exit('failed test4')
        print 'passed GauntCoeff test'

def ClebschGordan(L,M,l1,m1,l2,m2):
    """dummy, to be replaced by full function"""

    if L!=0: exit('only implemented for L=M=0')
    
    if M!=m1+m2: return 0
    if L>l1+l1: return 0
    if L<abs(l1-l2): return 0

    if m1>l1: exit('illegal: m1>l1')
    if m2>l2: exit('illegal: m2>l2')
    if M > L: exit('illegal: M > L')

    return 1.

class BasTwoAngle:
    """angular function with its list of radial blocks"""
    def __init__(self,L,M,l1,l2):
        """List of Clebsch-Gordan coefficients and Ylm's for 
        |L,M,l1,l2> = sum[m1,m2] C(L,M,l1,m1,l2,m2) Ylm(l1,m1) Ylm(l2,m2)
        """
        self.L=L
        self.M=M
        self.l1=l1
        self.l2=l2
        self.CG=[]
        self.m1=[]
        self.m2=[]
        for m in range(-l1,l1+1):
            cg=ClebschGordan(L,M,l1,m,l2,M-m)
            if cg!=0: 
                self.CG.append(cg)
                self.m1.append(m)
                self.m2.append(M-m)
        self.brad=[]

    def __str__(self):
        string='L,M,l1,l2= '+str(self.L)+','+str(self.M)+','+str(self.l1)+','+str(self.l2)
        for br in self.brad:
            string+=str(br)
        return string

    def x12(self):
        """1-2 exchanged angular function"""
        b=copy(self)
        b.L =self.L
        b.M =self.M
        b.l1=self.l2
        b.l2=self.l1
        b.CG=self.CG
        b.m1=self.m2
        b.m2=self.m1
        return b

class BasTwoRadial:
    """radial elements with a list of function indices"""

    def __init__(self,e1,e2):
        """|i>=|L,M,l0,l1> |e.k0> |e.k1>
        |L,M,l0,l1>
        |e.k> ...radial function from FE basis, multiplied by r^l if lower boundary =0
        """
        self.e1=e1
        self.e2=e2
        self.k1=[]
        self.k2=[]
        self.i=[]

    def __str__(self):
        string='\n'+str(self.e1)+' | '+str(self.e2)+'\n' 
        for i in range(len(self.i)):
            string+=str(self.i[i])+','+str(self.k1[i])+','+str(self.k2[i])+'\n'
        return string

    def x12(self):
        """1-2 exchanged radial block"""
        b=copy(self)
        b.e1=self.e2
        b.e2=self.e1
        b.k1=self.k2
        b.k2=self.k1
        b.i =self.i
        return b
    
class BasTwo:
    """list of angular functions, each with its list of radial blocks
    """

    def __init__(self,xsym,ax,Lmax,Mmax,lmax,parity='natural',ax2=None,psum=None,Lmin=None,Mmin=None,lsum=None):
        """two particle basis
        xsym=+1,-1,0 ... exchange symmetry
        parity=natural or unnatural
        ax  ... coordinate axis for radial functions
        constraints:
        L in [Lmin,Lmax], default Lmin=Lmax
        M in [Mmin,Mmax], default Mmin=Mmax
        l1 in [0,lmax]
        l2 in [0,lsum-l], default lsum=lmax
        ax2 ... second, possibly shorter axis, default ax2=ax
        psum... inner element basis functions contstrained to k1+k2<=psum
        """                                  

        #  set the defaults
        def default_if_None(val,deflt): 
            if val is None: 
                return deflt
            else:
                return val

        self.xsym=xsym
        self.ax  =ax
        self.bx  =default_if_None(ax2,ax)
        self.ax.show('radial axis 1')
        self.bx.show('radial axis 2')
        L0=default_if_None(Lmin,Lmax)
        M0=default_if_None(Mmin,Mmax)
        ls=default_if_None(lsum,2*lmax)+1
        ks=default_if_None(psum,self.ax.order()+self.bx.order())

        self.gaunt=GauntCoeffTable(2*lmax)
        self.bang=[]
        self.len=-1
        block_i0=0
        count=0
        for L in range(L0,Lmax+1):
            for M in range(M0,Mmax+1):
                for l1 in range(lmax+1):
                    for l2 in range(lmax+1):
                        if l1+l2>ls: continue
                        if parity=='natural'   and (L+l1+l2)%2==1: continue
                        if parity=='unnatural' and (L+l1+l2)%2==0: continue
                        if xsym!=0 and l1<l2: continue # skip exchange symmetric angular part
                        self.bang.append(BasTwoAngle(L,M,l1,l2))
                        ba=self.bang[-1]

                        # generate product basis
                        for e1 in self.ax.e:
                            for e2 in self.bx.e:
                                ba.brad.append(BasTwoRadial(e1.centrifugal(l1),e2.centrifugal(l2)))
                                br=ba.brad[-1]
                                for k1 in range(e1.n):
                                    for k2 in range(e2.n):
                                        count+=1
                                        br.k1.append(k1)
                                        br.k2.append(k2)
                                        itotal=block_i0+e1.i0+k1+self.ax.len()*(e2.i0+k2)
#                                        print 'block',L,M,l1,l2,itotal,block_i0
                                        self.len=max(self.len,itotal+1)
                                        br.i.append(itotal)
                        block_i0=block_i0+self.ax.len()*self.bx.len()        
        print 'total',self.len

    def __str__(self,text=None):
        string='exchange '+str(self.xsym)+'\n'
        string+='axes: '+self.ax.__str__('1st')+self.ax.__str__('2nd')
        for ba in self.bang:
            string+=str(ba)
        return string

    def hamiltonian(self):
        """hamiltonian matrix for basis"""

        # compute the radial two-particle integrals
        class IntegTwoTab:
            """table to two-electron integrals"""

            def __init__(self,e1,e2,lmax):
                """basic two-particle integrals for multipole expansion"""
                self.basic=[]
                for e in e1:
                    for f in e2:
                        self.basic.append([e,f,coulomb_integrals_radial(e,f,lmax)])


            def compose(self,bai,baj,e1,e2,gaunt):
                """compose specific integral table from basic tables"""
                for tab in self.basic: 
                    if [e1,e2]==tab[:2]: break
                else:
                    exit('no basic table for radial sets'+str(e1)+str(e2))

                ints=zeros((shape(tab[2])[1]))
                lmin=max(abs(bai.l1-baj.l1),abs(bai.l2-baj.l2))
                lmax=min(bai.l1+baj.l1,bai.l2+baj.l2)
                for i in range(len(bai.CG)):
                    l1i,l2i,m1i,m2i=bai.l1,bai.l2,bai.m1[i],bai.m2[i]
                    for j in range(len(baj.CG)):
                        l1j,l2j,m1j,m2j=baj.l1,baj.l2,baj.m1[j],baj.m2[j]
                        cgij=bai.CG[i]*baj.CG[j]
                        for l in range(lmin,lmax+1):
                            gf1=gaunt.coeff(l1j,l,l1i,m1j,m1j-m1i,m1i)
                            if gf1==0: continue
                            gf2=gaunt.coeff(l2i,l,l2j,m2i,m2i-m2j,m2j)
                            if gf2==0: continue
                            ints+=4*myPi/(2*l+1)*gf1*gf2*cgij*tab[2][l,:]

                # convert to convenient format
                Vee=zeros((e1.n,e2.n,e1.n,e2.n))
                ijkl=-1
                for i in range(e1.n):
                    for j in range(i+1):
                        for k in range(e2.n):
                            for l in range(k+1):
                                ijkl+=1
                                Vee[i,k,j,l]=ints[ijkl]
                                Vee[j,k,i,l]=ints[ijkl]
                                Vee[i,l,j,k]=ints[ijkl]
                                Vee[j,l,i,k]=ints[ijkl]
                return Vee

        def addHmatSmat(Hmat,Smat,bai,baj,bri,brj,xsym=1):
            """add to global matrix for given radial block"""

            # shorter notation 
            e1i,e2i,e1j,e2j=bri.e1,bri.e2,brj.e1,brj.e2

            # return if same particle matrix elements do not overlap
            if not e1i.overlap(e1j) or not e2i.overlap(e2j): return
                
            # compose two-particle integrals for current angular and radial blocks
            Vee=int2table.compose(bai,baj,e1i,e2i,gaunt)

            Hmat1=0.5*e1i.matrix('1/qdq|1/qdq',e1j)+0.5*(bai.l1*(bai.l1+1))*e1i.matrix('|1/q^2|',e1j)-charg1*e1i.matrix('|1/q|',e1j)
            Hmat2=0.5*e2i.matrix('1/qdq|1/qdq',e2j)+0.5*(bai.l2*(bai.l2+1))*e2i.matrix('|1/q^2|',e2j)-charg2*e2i.matrix('|1/q|',e2j)
#            Hmat1=0.5*e1i.matrix('1/qdq|1/qdq',e1j)+0.5*(bai.l1*(bai.l1+1))*e1i.matrix('|1/q^2|',e1j)+0.5*e1i.matrix('|q^2|',e1j)
#            Hmat2=0.5*e2i.matrix('1/qdq|1/qdq',e2j)+0.5*(bai.l2*(bai.l2+1))*e2i.matrix('|1/q^2|',e2j)+0.5*e2i.matrix('|q^2|',e2j)
                
            Smat1=e1i.matrix('|',e1j)
            Smat2=e2i.matrix('|',e2j)

            def vee_idx(k1i,k2i,k1j,k2j,nd1_sq):
                kk1=(max(k1i,k1j)*(max(k1i,k1j)+1))/2+min(k1i,k1j)
                kk2=(max(k2i,k2j)*(max(k2i,k2j)+1))/2+min(k2i,k2j)
                return kk1+nd1_sq*kk2

            ikjk=-1
            for jk in range(len(brj.i)):
                k1j,k2j,j=brj.k1[jk],brj.k2[jk],brj.i[jk]
                for ik in range(len(bri.i)):
                    k1i,k2i,i=bri.k1[ik],bri.k2[ik],bri.i[ik]
                    Hmat[i,j]+=xsym*Vee[k1i,k2i,k1j,k2j]
                    Hmat[i,j]+=xsym*Hmat1[k1i,k1j]*Smat2[k2i,k2j]
                    Hmat[i,j]+=xsym*Smat1[k1i,k1j]*Hmat2[k2i,k2j]
                    Smat[i,j]+=xsym*Smat1[k1i,k1j]*Smat2[k2i,k2j]
            for x in la.eigvals(Smat1):
                if x<1.e-7: exit('ill-conditioned single particle matrix')
            for x in la.eigvals(Smat2):
                if x<1.e-7: exit('ill-conditioned single particle matrix')
        
        # get all elements on first and second coordinate
        print 'computing basic 2e integral tables...'
        int2table=IntegTwoTab(self.ax.e,self.bx.e,2*lmax)
        Hmat=zeros((self.len,self.len))
        Smat=zeros((self.len,self.len))
        for li in range(len(self.bang)):
            bai=self.bang[li]
            for lj in range(li+1):
                baj=self.bang[lj]

                # rotationally symmetric - skip if angular momenta differ
                if bai.L!=baj.L: continue
                if bai.M!=baj.M: continue
                print 'block l1,l2',bai.l1,bai.l2,baj.l1,baj.l2
                for ri in range(len(bai.brad)):
                    bri=bai.brad[ri]
                    for rj in range(len(baj.brad)):
                        if li==lj and ri<rj: continue #  skip symmetric part of the basis
                        brj=baj.brad[rj]

                        addHmatSmat(Hmat,Smat,bai,baj,bri,brj)
                        if self.xsym!=0: # exchange symmetric terms
                            addHmatSmat(Hmat,Smat,bai,baj.x12(),bri,brj.x12(),self.xsym)

        return Hmat,Smat

lmax=8
gaunt=GauntCoeffTable(lmax)
gaunt.test()
charg1=2
charg2=2
ax=Axis('r',20,0.,5.,order=2)
print 'axis length',ax.len()
bas=BasTwo(0,ax,0,0,2)
#print bas
set_printoptions(precision=4,suppress=False,linewidth=132)
Hmat,Smat=bas.hamiltonian()

plt.matshow(Hmat)
plt.show()

print 'solving eigenproblem...'
val=la.eigvals(Hmat,Smat)
ls = argsort(val.real)
print 'eig',val[ls[:20]].real
