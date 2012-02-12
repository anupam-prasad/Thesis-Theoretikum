#!/usr/bin/env python
"""laser field and vector potential generator"""
import sys 
import numpy as np
from math import *
from read_input import *
import my_pyplot # pyplot with include guard

# import the unit conversion
from my_constants import *
from units import *

# to do:
# class LaserPulse:
# a list of single laser pulses, using SinglePulse.read()
# whose fields/vector potentials are added

# the laser class
class SinglePulse:
    """complete definition of a laser pulse
    shape    ... pulse shape
    intensity... in W/cm2
    FWHM    .... full width half maximum of intensity
    omega   .... carrier frequency(a.u.)
    """

    def __init__(self,wavelength,duration_fwhm,peak_intensity,shape,phase=0,delay=0.):

        def envelope_gauss(Omega_t): return exp(-Omega_t**2),-2*exp(-Omega_t**2)
        def envelope_cossq(Omega_t):
            if Omega_t<-myPi/2 or Omega_t>myPi/2: return 0.,0.
            return cos(Omega_t)**2,-2*cos(Omega_t)*sin(Omega_t)

        self.shape=shape
        self.phase=phase
        self.delay=delay
        self.fwhm=duration_fwhm
        self.omega=SI().velocity(speed_of_light)/wavelength*2*myPi
        self.apeak=sqrt(peak_intensity)/self.omega
        
        if   shape=='gauss':
            self.Omega=sqrt(2.*log(2.))/duration_fwhm
            self.envelope=envelope_gauss
        elif shape=='cossq':
            self.Omega=0.5*myPi/duration_fwhm
            self.envelope=envelope_cossq
        else: exit('undefined pulse shape: '+shape)
        print 'Omega',self.Omega,duration_fwhm


    def field_pot(self,t): 
        e,d=self.envelope(self.Omega*(t-self.delay))
        arg=self.omega*(t-self.delay)+self.phase
        sin_arg=sin(arg)
        cos_arg=cos(arg)
        ap =self.apeak*(e*sin_arg)
        ef =self.apeak*(d*sin_arg*self.Omega+e*cos_arg*self.omega)
        return ef,ap

    def efield(self,t): return self.field_pot(t)[0]
    def vecpot(self,t): return self.field_pot(t)[1]

    # to do:
    # def envelope_trapez(Omega_t):
    # one more envelope  
    #      !     ^ E(t)
    #      !     |
    #      !     |  /|-----|\
    #      !     | / |     | \
    #      !      /--|-----|--\----> t
    #      !       a    b    c
    # NOTE: for the FIELD! needs integration
         
    # def field_pot_file():
    # in __init__:
    # - get a table for the FIELD from file
    # - integrate to also get VECTOR POTENTIAL from file
    # - set up an interpolation table: QUADRATIC interpolation for vector potential

    # to do:
    # def check(self): check the relation A(t)= -int[-\infty,t] E(t') dt'
    #                  by explicit numerical integration (to detect errors in the pulse definitions)


    def show(self):
        """plot the pulse with t = t_peak+-3*FWHM"""
        tt=[]
        ff=[]
        aa=[]
        for t in np.arange(-2*self.fwhm,2*self.fwhm,4*self.fwhm/401):
            f,a=self.field_pot(t)
            tt.append(t)
            ff.append(f)
            aa.append(a)
            
        plt.plot(tt,ff)
        plt.show()
        plt.plot(tt,aa)
        plt.show()

    @classmethod
    def read(cls):
        """read definition of one or several pulses"""
        i=0
        wave_length=   read_input(800.,   '.laser pulse',1,i+1,doc='carrier wave length (nm)')        
        duration=      read_input(1.,     '.laser pulse',2,i+1,doc='FWHM of pulse duration (optical cycles)')
        peak_intensity=read_input(1.e-14, '.laser pulse',3,i+1,doc='peak intensity (W/cm^2)')
        pulse_shape=   read_input('cossq','.laser pulse',4,i+1,doc='shapes: cossq,gauss,trapez,file,...')
        phase=         read_input(0.     ,'.laser pulse',5,i+1,doc='carrier-envelope offset phase (optical cycles)')
        delay=         read_input(0.     ,'.laser pulse',6,i+1,doc='time delay of the pulse (fs)')

        # convert to internal units
        Units(au())
        print 'wave length',SI().length(wave_length*1.e-9)
        print 'speed of light',SI().length(speed_of_light)
        print 'duration',OptCyc(wave_length*1.e-9).time(duration)
        laser=SinglePulse(SI().length(wave_length*1.e-9),
                        OptCyc(wave_length*1.e-9).time(duration),
                        SI().intensity(peak_intensity*1.e-4),
                        pulse_shape,
                        phase/(2.*myPi),
                        SI().time(delay*1.e-9)
                        )
        return laser

    @classmethod
    def test(cls):
        read_input_open("laser.inp")
        laser=SinglePulse.read()
        laser.show()

if __name__  == "__main__":
    """standalone: tests"""
    SinglePulse.test()
