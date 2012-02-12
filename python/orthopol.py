#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.special.orthogonal as so # has the quadratur rules
import myfunction as mf               # separate home-made module with class MyFunction

class OrthogonalPolynomial(mf.MyFunction):
    """
    values and derivatives of orthogonal polynomials
    mf.MyFunction has method "plot()"
    """
    def val(self,q,n=0):
        """values and derivatives up to ORDER n (degree n-1)
        """
        # values and derivatives
        v=np.zeros((n))
        d=np.zeros((n))
        v[0]=1.
        d[0]=0.
        if n>0:
            for i in range(1,n):
                v[i]=(self.a(i)+self.b(i)*q)*v[i-1]                 +self.c(i)*v[i-2]
                d[i]=(self.a(i)+self.b(i)*q)*d[i-1]+self.b(i)*v[i-1]+self.c(i)*d[i-2]
        return v,d

    def test(self,n):
        """check orthogonality up to order n
        """
        # get the quadrature rule 
        x,w=self.quadrature(n)

        # print it - should be diagonal
        np.set_printoptions(precision=5,suppress=True)

        # compute the overlap matrix
        o=np.zeros((n,n))
        for j in range(len(x)):
            v=self.val(x[j],n)[0]
            for i in range(n):
                o[i,:]=o[i,:]+v[i]*v[:]*w[j]
        print o

        # plot the functions
        self.plot(-1.,1.,n)
        plt.show()

class LegendrePolynomial(OrthogonalPolynomial):
    """recurrence coefficients and quadrature
    inherits from OrthogonalPolynomial:
    val()... return values and derivatives
    plot().. plot on reasonable interval
    """
    def a(self,i): return 0.
    def b(self,i): return (2.*float(i)-1.)/float(i)
    def c(self,i): return (-float(i)+1.)/float(i)
    def quadrature(self,n): return so.p_roots(n)

class LaguerrePolynomial(OrthogonalPolynomial):
    """recurrence coefficients and quadrature
    """
    def a(self,i): return (2*float(i)-1.)/float(i)
    def b(self,i): return -1./float(i)
    def c(self,i): return (-float(i)+1.)/float(i)
    def quadrature(self,n): return so.la_roots(n,0)

class ChebyshevPolynomial(OrthogonalPolynomial):
    """recurrence coefficients and quadrature
    """
    def a(self,i): return 0.
    def b(self,i): 
        if i==1: return 1. # Chebyshevs are normalized inconsistently
        else: return 2.
    def c(self,i): return -1.
    def quadrature(self,n): return so.t_roots(n)

class HermitePolynomial(OrthogonalPolynomial):
    """recurrence coefficients and quadrature
    """
    def a(self,i): return 0.
    def b(self,i): return 2.
    def c(self,i): return -2.*float(i-1)
    def quadrature(self,n): return so.h_roots(n)

#print ''
#print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
#print '!!! WARNING: OrthogonalPolynomial was redefined for consistency with BasisFunction'
#print '!!! now n is the order = maximal degree + 1 '
#print '!!! i.e. n is the number of functions (n was the degree before) '
#print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
#print ''
if __name__  == "__main__":
    """standalone: tests"""

    np.set_printoptions(linewidth=132)
    print 'Gauss-Hermit quadratur rules'
    for J in range(2,7):
        print 'x[',J,']=',so.h_roots(J)[0]
        print 'w[',J,']=',so.h_roots(J)[1]

    print 'Chebyshev'
    ChebyshevPolynomial().test(5)
    print 'Laguerre'
    LaguerrePolynomial().test(5)
    print 'Legendre'
    LegendrePolynomial().test(5)
    print 'Hermite'
    HermitePolynomial().test(5)
