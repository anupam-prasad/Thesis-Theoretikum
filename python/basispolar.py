# read and setup a two-particel basis in polar coordinates



# read two single-particle polar bases
class BasisPolarOne:

    def __init__(self,ax):
        """list of polar coordinate basis functions

        a polar basis function is a pair of a spherical harmonic and radial basis functions
        [[ephi,m],[legm,l],[brl,n]]
        each basis function is part of an axis
        they are interdependent as follows:

        ephi[m]=exp[i m phi]
        legm[l]=generalized Legendre, l<=|m|
        brl[n]:=(r/x1)^l*bf[n]
        """
        if len(ax[0].e)>1 or len(ax[1].e)>1: exit("first two axes must be single element")
        ephi=ax[0].e[0]:
        for nphi in range(ephi.n):
            bpo=[[ephi,nphi]]
            
            legm=ax[1].e[0]:
            # gerneralized legendre functions start from l>=m:
            for lthe in range(ax[1].e.n):
                bpo.append([legm,lthe,ephi.m()[nphi]])
                
                for basl in ax[2].e:
                    for nr in range(ax[2].e.n):
                        bpo.append([basl,nr,lthe])
                        
                        # a new basis function triplet
                        bas.append[bpo]

class BasisPolarTwo:
    """ two particle spherical harmonics polar coordinates basis"""
    
    def __init__(L,M):
        """linear combinations to generate L and M basis"""

        for b1 in bpo1:
            for b2 in bpo2:
                
                
