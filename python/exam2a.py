#!/usr/bin/env python
import sys 
import numpy as np

class BlockMatrix:
    """block matrix class"""
    
    def __init__(self,iblocks,jblocks,mat=None):
        """block matrix
        iblocks ... block sizes on the row-index
        jblocks ... block sizes on the column-index"""
        self.I=iblocks
        self.J=jblocks
        if mat is not None:
            if (sum(self.I)!=np.shape(mat)[0] or
                sum(self.J)!=np.shape(mat)[1]): exit ('block sizes do not add up to dimensions')

        # create list of matrix blocks
        self.b=[]
        i0=0
        for i in self.I:      # loop through row blocks
            i1=i0+i
            j0=0
            for j in self.J:  # loop through column blocks
                j1=j0+j
                if mat is None: self.b.append(np.matrix(np.zeros((i,j))))
                else:
                    self.b.append(np.matrix(mat[i0:i1,j0:j1]))
                j0=j1
            i0=i1

    def ij(self,i,j): 
        """convert a double index in to linear index in list of blocks"""
        return i*len(self.J)+j

    def __add__(self,b):
        """returns sum of two matrices"""
        if any(self.I!=b.I): exit('left hand blocking does not match')
        if any(self.J!=b.J): exit('right hand blocking does not match')
        c=BlockMatrix(self.I,self.J)
        for i in range(len(self.b)):
            c.b[i]=self.b[i]+b.b[i]
        return c
    
    def __mul__(self,b):
        """matrix product of two block-matrices
        returns block-matrix"""
        if any(self.J!=b.I): exit('joint index does not match, cannot multiply')
        c=BlockMatrix(self.I,b.J)
        ij=-1
        for i in range(len(self.I)):
            for j in range(len(b.J)):
                ij+=1
                for k in range(len(self.J)):
                    c.b[ij]+=self.b[self.ij(i,k)]*b.b[b.ij(k,j)]
        return c

    def trans(self):
        c=BlockMatrix(self.J,self.I)
        for i in range(len(self.I)):
            for j in range(len(self.J)):
                c.b[c.ij(j,i)]=self.b[self.ij(i,j)].T
        return c

    def full(self):
        """merge block matrix into a full matrix"""
        a=np.matrix(np.zeros((sum(self.I),sum(self.J))))
        ij=-1
        i0=0
        for i in self.I:
            i1=i0+i
            j0=0
            for j in self.J:
                j1=j0+j
                ij+=1
                a[i0:i1,j0:j1]=self.b[ij]
                j0=j1
            i0=i1
        return a

    @classmethod
    def check(cls):
        """check routines for random matrices"""
        I=np.array([3,4,7])
        K=np.array([2,4])
        J=np.array([5,6])
        
        
        # test full
        A=np.matrix(np.random.random((sum(I),sum(K))))
        bA=BlockMatrix(I,K,A)
        if (abs(A-bA.full())>1.e-10).any(): exit('full() failed')
        print 'full() passed'

        # test +
        B=np.matrix(np.random.random(np.shape(A)))
        bB=BlockMatrix(I,K,B)
        if (abs((bA+bB).full()-(A+B))>1.e-10).any(): exit('+ failed')
        print "+ test passed"
        
        # test *
        B=np.matrix(np.random.random((sum(K),sum(J))))
        bA=BlockMatrix(I,K,A)
        bB=BlockMatrix(K,J,B)
        if abs((bA*bB).full()-(A*B)>1.e-10).any(): exit('* failed')
        print "* test passed"

        # test trans
        if (abs(bA.trans().full())-A.T>1.e-10).any(): exit('trans failed')
        print 'trans() test passed'

if __name__ == "__main__": 
    print "test"
    BlockMatrix.check()
