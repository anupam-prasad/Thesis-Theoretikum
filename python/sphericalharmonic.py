# spherical harmonics

class Ysphhar:
    def __init__(self,l,m):
        if abs(m)>l: exit('ILLEGAL |m|>l in in Y_lm')
        self.l=l
        self.m=m

    def val(phi,cth):
        for m in range(-self.m,self.m+1):
            mv,md=
            lv,ld=LegendreGeneralized(self.m,self.l)
            
