#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from read_input import *
from axis import *
from tensor import *
from myeigen import inversiter

# open the first string after the program name as input file 
read_input_open(sys.argv[1])

# read system parameters
field=read_input(0.,'.system parameters',1,doc='field strength (a.u.)')
lmin =read_input(0, '.system parameters',2,doc='minimal angular momentum')
lmax =read_input(0, '.system parameters',3,doc='maximal angular momentum')
mqn  =read_input(0, '.system parameters',4,doc='m quantum number')
steps=read_input(0,'.field steps',1,doc='increase to field strength in steps')
nplot=read_input(0,'.field steps',2,doc='number of energies to plot')
metho=read_input('full','.field steps',3,doc='how to find roots: full,invit')
theta=read_input(0.,'.complex scaling angle',1,doc='how to find roots: full,invit')


lmin=max(mqn,lmin) # minimal angular momentum cannot be below m-quantum number

ax=axis_read(0) # read the axis

# check input, write docu-file
read_input_finish()

# hydrogen atom (L=0)
print 'hydrogen atom'
hamr=np.zeros((ax.n,ax.n))
kinr=np.zeros((ax.n,ax.n))
potr=np.zeros((ax.n,ax.n))
angr=np.zeros((ax.n,ax.n))
dipr=np.zeros((ax.n,ax.n))
ovrr=np.zeros((ax.n,ax.n))
for e in ax.e:
    kinr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('1/qdq|1/qdq')
    potr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('coulomb(1.)')
    angr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=0.5*e.matrix('|1/q^2|')
    dipr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=field*e.matrix('|q|')
    ovrr[e.i0:e.i0+e.n,e.i0:e.i0+e.n]+=e.matrix('|')
hamr=kinr+potr

ovro=np.eye((lmax+1-lmin))
ango=np.zeros((lmax+1-lmin,lmax+1-lmin))
dipo=np.zeros((lmax+1-lmin,lmax+1-lmin))
ango[0,0]=lmin*(lmin+1)
for l in range(lmin,lmax+1):
    ango[l-lmin,l-lmin]=l*(l+1)
    dipo[l-lmin-1,l-lmin  ]=sqrt((l*l-mqn*mqn)/float(4*l*l-1))
    dipo[l-lmin,  l-lmin-1]=dipo[l-lmin-1,l-lmin]

Ham=Tensor(TensorProd(hamr,ovro))+Tensor(TensorProd(angr,ango))+Tensor(TensorProd(dipr,dipo))
Ovr=Tensor(TensorProd(ovrr,ovro))

(val,vec)=la.eig(Ham.array().real,Ovr.array().real)      # solve eigenproblem

np.set_printoptions(precision=9,suppress=True,linewidth=132)
print np.sort(val.real)[:10]   # show results (real part of eigenvalues, sorted !)

Tmat=exp(-2j*theta)*(Tensor(TensorProd(kinr,ovro))+Tensor(TensorProd(angr,ango))).array()
Vmat=exp(-1j*theta)*Tensor(TensorProd(potr,ovro)).array()
Hmat=Tmat+Vmat
Dmat=exp( 1j*theta)*Tensor(TensorProd(dipr,dipo)).array()/field
Omat=Ovr.array().real


en1=zeros((steps+1))
en2=zeros((steps+1,2))
fiel=[]
ener=[]
for j in range(nplot): ener.append([]) # list of empty lists 

if metho == 'full':
    for i in range(steps+1):
        fiel.append(i*field/steps)
        val=la.eigvals(Hmat+(i*field/steps)*Dmat,Omat)
        ls=np.argsort(val.real)
        for j in range(nplot):
            ener[j].append(val[ls[j]])
        
        # find energies nearest to -0.5
        dif1=0.5
        for j in range(nplot): 
            if abs(ener[j][i]+0.5)<dif1:
                dif1=abs(ener[j][i]+0.5)
                en1[i]=ener[j][i]

        eguess=-0.125
        dif2=0.125
        dif3=0.125
        for j in range(nplot): 
            if i>3: eguess=-0.125+float(i)/3.*(en2[3,0]+0.125)
            if abs(ener[j][i]-eguess)<dif2:
                dif2=abs(ener[j][i]-eguess)
                en2[i,0]=ener[j][i]
        
elif metho=='invit':
    val=np.sort(la.eig(Hmat+field/steps*Dmat,Omat)[0].real)[:nplot]
    print 'initial',val
    for j in range(nplot):
        wf=np.random.random((np.shape(Hmat)[0]))
        ee=val[j]
        for i in range(1,steps+1):
            if j==0: fiel.append(i*field/steps)
            ee,wf=inversiter(ee,wf,Hmat+i*field/steps*Dmat,Omat)
            ener[j].append(ee)
    
else: exit('root finding undefined: '+metho)

print en1[:10]
print en2[:10,0]

ls=np.argsort(abs(val.imag))
print val[ls[:20]]
plt.plot(val[ls[:100]].real,val[ls[:100]].imag,'o')
plt.show()

for j in range(nplot): 
    plt.plot(np.array(fiel),np.array(ener[j]))
plt.plot(np.array(fiel),en2[:,0])
plt.show()

plt.plot(np.array(fiel),en1)
plt.plot(np.array(fiel),en2[:,0])
plt.show()

plt.loglog(np.array(fiel[1:]),-0.5-en1[1:])
#plt.loglog(np.array(fiel[1:]),abs(-0.125-en2[1:,0]))
plt.show()

