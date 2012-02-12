#!/usr/bin/env python
"""a general abstract matrix class with matrix operations implemented using LAPACK/BLAS"""
import sys 
import numpy as np
import scipy.linalg.flapack as lapack
import mytimer

class GMatrix:
    """general abstract matrix class
    handles: real full symmetric, real banded symmetric
    extendible to any type of matrics
    """
    def __init__(self,mat,rows=None,cols=None,Blo=None,Bup=None,storage='full'):
        """a matrix has 
        rows,cols ...numbers of rows,columns
        bup,blo   ...upper and lower band width 
        store     ...storage format: full, upper banded (... to be extended)
        symm      ...symmetry: True or False
        (hermit)  ...hermitian, to be implemented
        """

        self.m=mat
        self.store=storage
        if self.store is None: self.store='full'

        # fill in the matrix parameters
        self.rows=rows
        self.cols=cols
        if self.store is 'full':
            self.symm=None
            if rows is None: self.rows=np.shape(mat)[0]
            if cols is None: self.cols=np.shape(mat)[1]
            self.blo=None
            self.bup=None
                    
        elif self.store is 'upper_banded':
            self.symm=True
            self.rows=self.cols=np.shape(mat)[1]
            self.bup=self.blo =np.shape(mat)[0]-1

        else: sys.exit('undefined storage type: '+storage)

        # mapping from matrix indices to storage
        def map_full((i,j)):         
            return (i,j)
        def map_upper_banded((i,j)):
            if i<j: return (self.bup-j+i,j)
            else:   return (self.bup-i+j,i)

        if   self.store=='full': self.map=map_full
        elif self.store=='upper_banded': self.map=map_upper_banded
        else: sys.exit('map not defined for storage= '+self.store)

    def __getitem__(self,ij):
        """return matrix element self[i,j]"""
        return self.m[self.map(ij)]

    def __setitem__(self,ij,val): 
        """assign to matrix element self[i,j]=val"""
        self.symm=None # assigement may destroy symmetry
        self.m[self.map(ij)]=val

    def __str__(self): 
        string=str(self.rows)+' '+str(self.cols)+' '+str(self.bup)+' '+str(self.blo)+' '+self.store
        string+='\n'+str(self.m)
        return string

    def irow(self,j):
        """range of row indices in j'th column"""
        if self.store=='full': return range(self.rows)
        elif self.store=='upper_banded': return range(max(0,j-self.bup),min(self.rows,j+self.bup+1))
        else: sys.exit('irow() not defined for storag: '+self.store)

    def bandwidth(self):
        """determine upper and lower band width of a matrix
        = number of super- and sub-diagonals"""

        if self.store!='full': sys.exit('band-width determinantion only for storage = "full"')

        # upper band width
        if self.bup is None:
            self.bup=self.cols
            while self.bup>0:
                self.bup-=1
                for i in range(self.cols-self.bup):
                    if self[i,i+self.bup]!=0: nonzero=True
                    else: nonzero=False;break
                if nonzero: break
                    
                if self.symmetric():
                    self.blo=self.bup
                    return
                    
        # lower band width
        if self.blo is None:
            self.blo=self.rows
            while self.blo>0:
                self.blo-=1
                for i in range(self.rows-self.blo):
                    if self[i+self.blo,i]!=0: nonzero=True
                    else: nonzero=False;break
                if nonzero: break

    def symmetric(self):
        """True if matrix is symmetric, i.e. M[i,j]=M[j,i] for all i,j """
        if self.symm != None: return self.symm
        for store in ['upper','lower','upper_banded','lower_banded']:
            if self.store==store: return True
        if self.rows!=self.cols: return False
        for i in range(self.cols):
            for j in range(i):
                if self[i,j]!=self[j,i]: return False
        else:
            return True

    def restore(self,storage):
        """copy matrix into a new matrix with the given storage format
        """
        # NOTE: a lot of checks should be introduced here!

        if storage=='full': m=np.zeros((self.rows,self.cols))
        elif storage=='upper_banded':m=np.zeros((self.bup,self.cols))
        else: sys.exit('undefined storage type: '+storage)

        a=GMatrix(m,self.rows,self.cols,self.blo,self.bup,storage)

        for j in range(self.cols):
            for i in self.irow(j):
                a[i,j]=self[i,j]
        return a

    def eigvals(self): 
        """eigenvalues of the matrix"""
        return self.eig(vectors=False) # solve eigenproblem without computing eigenvectors 

    def eig(self,vectors=True):
        """solve matrix eigenproblem
        vectors ...True = eigenvalues and eigenvectors, False = eigenvalues only
        """
        if vectors: right_v=1
        else: right_v=0

        # note: compute_v=0 means now eigenvectors will be calculated

        if self.store is 'full':
            if self.symmetric(): w,v,info=lapack.dsyev(self.m,compute_v=right_v)
            else: sys.exit('no diagonalization defined for non-symmetric matrix')

        elif self.store is 'upper_banded':
            w,v,info=lapack.dsbev(self.m,compute_v=right_v)

        else: sys.exit('eigvals not defined for storage type: '+self.store)

        if info!=0: sys.exit('error in lapack eig')
        return w


def compare_eigen_methods():
    """ timing of different eigensolver methods"""
    import scipy.linalg as linalg 
    print '\n *** diagonalize a real symmetric band matrix by different methods ***\n'
    mt=mytimer.create(10)    
    try:
        N=float(sys.argv[1])
        S=float(sys.argv[2])
    except:
        sys.exit('supply N S (command line arguments): matrix dimension N and and number of super-diagonals S ')

    np.random.seed(1)
    a=np.random.random((N,N))

    ab=GMatrix(np.random.random((S+1,N)),storage='upper_banded')
    mt[0].start('lapack symmetric upper banded')
    print ab.store+'\n',ab.eigvals()[:5]
    mt[0].stop()

    a=ab.restore('full')
    mt[1].start('lapack symmetric full')
    print a.store+'\n',a.eigvals()[:5]
    mt[1].stop()

    print 'lingalg'
    mt[2].start('linalg general')
    print np.sort(linalg.eigvals(a.m).real)[:5]
    mt[2].stop()

    mytimer.table()


if __name__ == "__main__":
    """standalone: tests"""
    compare_eigen_methods()

    
