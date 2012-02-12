import sys
import numpy as np
import scipy.special.orthogonal as so

quad_table={}     # keep quadrature rules in this dictionary

def rule(n,a=0.,b=1.,kind='legendre'):
    # if table entry does not exist, generate it
    try:
        x,w=quad_table[kind+'_'+str(n)]
    except:
        if kind=='legendre': x,w=so.p_roots(n)
        else: sys.exit('undefined quadrature kind "'+kind+'"')

        # append to table
        quad_table.setdefault[kind+'_'+str(n),(x,w)]


    # return shifted and scaled rule
    return a+x*(b-a),w*(b-a)
