import numpy as np

class MyMatrix:
    def __init__(self,m,store=None,symm=None,spars=None,copy=False):
        """matrix plus its structural information
        m ...data storage
        symm ...'general','symmet','hermit','utriang','ltriang'
        nsup ...number of super-diagonals
        nsub ...number of sub-diagonals
        store...form of storage: 'full','triangular','banded',
        """
        self.m=m
        self.type='unknown'
        self.ub=min(m.shape())
        self.lb=max(m.shape())
        self.store='full'

    def QR(self):
        """overwrite matrix with its QR decomposition"""
        qr=np.matrix(self.m)
        for i in range(qr.shape()[0]):
            w=qr[i:,i]
            if w[0]>0: w[0]+=abs(w.norm())
            else:      w[0]-=abs(w.norm())
            qr[i:,i+1:]-=2.*w*w.transpose()*qr[i:,i:]/(w.transpose*w)
            qr[i+1:,i]=w[1:]/w[0] # we do not store the exact Q factors, rather scaled

    def OrthoIter(A,r):
        """determine r largest eigenvalues by orthogonal iteration"""
        
