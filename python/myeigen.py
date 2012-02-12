#! /usr/bin/env python

# inverse iteration 

import sys
import numpy as np
from math import * 
import scipy.linalg as la
import mytimer as mt

def inversiter(eguess,wfguess,A,S,Klu=100,Kpower=1):

    vk=np.matrix(np.random.random((wfguess.size,1)))
    olu,opi=la.lu_factor(S)
    
    for L in range(Klu):
        # get the LU decompostion of A-e*S
        lu,piv=la.lu_factor(A-eguess*S)
        
        for K in range(Kpower):
            vk=np.matrix(S)*vk
            vk=np.matrix(la.lu_solve((lu,piv),vk,overwrite_b=False))
            ek=(vk.T*A*vk).diagonal()/(vk.T*S*vk)
        
        if (eguess-ek[0,0])/eguess<1.e-9: return ek[0,0],vk[:,0]
        eguess=ek[0,0]
    
