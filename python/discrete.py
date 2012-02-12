import sys

import numpy as np
import matplotlib.pyplot as plt

import scipy
import scipy.special as ss
import scipy.special.orthogonal as so
import scipy.linalg as la

import cmath
import myfunction as mf
import orthopol as op



class Grid:
    
    def __init__(self,x0,x1,n):
        
