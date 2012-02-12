class ODE:
    
    def runge_kutta(deriv,y0,t,h,stage=4):
        """runge-kutta solver for objects whose derivatives form a linear space"""
        for i in range(stage):
            y1=y0

            for j in range(i-1):
                y1=y1+ark[i,j]*yrk[j]
                
            yrk[i]=h*deriv(y1,t+h*crk[i])
        

    def solve(func,y,tend,tstart=None,meth='runge-kutta'):
        
        if hasattr(y,time): 
            tt=y.time
        else: 
            if tstart is None: sys.exit('object has no time of its own, specify tstart')
            tt=tstart

        if meth == 'runke-kutta': runge_kutta(func,y,tt,tend-tt)
        else: sys.exit('undefined ODE method')
