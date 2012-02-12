#!/usr/bin/env python
"""conversion between physical unit systems"""

import sys 
 # here is where ALL physical and mathematical constants are defined
from my_constants import *

LocalUnits=None
def if_None_local(to): 
    if to is not None: return to
    else: return LocalUnits

class Units:
    """convert between different unit systems"""

    def __init__(self,local=None):
        global LocalUnits
        if LocalUnits is not None: exit('Units() can be created only once in each code')
        LocalUnits=local
       
    def __str__(self):
        string =self._name
        return string

    def supplement(self):
        self._velocity= self._length/self._time
        self._energy=   self._mass*self._velocity**2
        self._frequency=1/self._time
        self._current  =self._charge*self._velocity
        self._ep0      =self._velocity**2/(speed_of_light**2*self._mu0)
        
        # for SI, the following 2 will be overwritten by Volt=1 and Watt=1
        self._efield=self._charge/self._length**2/(4.*myPi*vacuum_polarizability)
        self._intensity=self._efield**2/2.*speed_of_light*vacuum_polarizability

    def time(self,   x):      return x*self._time     /self.local._time
    def velocity(self,   x):  return x*self._velocity /self.local._velocity
    def length(self,   x):    return x*self._length   /self.local._length
    def mass(self,   x):      return x*self._mass     /self.local._mass
    def charge(self,   x):    return x*self._charge   /self.local._charge
    def intensity(self,   x): return x*self._intensity/self.local._intensity
    def frequency(self,   x): return x*self._frequency/self.local._frequency
    def current(self,   x):   return x*self._current  /self.local._current
    def energy(self,   x):    return x*self._energy   /self.local._energy
    def efield(self,   x):    return x*self._efield   /self.local._efield

    def test(self):
        Units(au()) # set the default output units
        print '\n *** a few known quantities computed by conversions ***\n'
        print '    speed of light (au):',SI().velocity(speed_of_light)
        print '          1au time (SI):',au(SI()).time(1.)
        print '        1au length (SI):',au(SI()).length(1.)
        print '    omega @ 800 nm (au):',nm_inv().energy(1/800.)
        print '800nm optical cyc. (au):',OptCyc(800.e-9).time(1.)
        print '800nm optical cyc. (au):',SI().time(800.e-9/speed_of_light)
        print '       1 au efield (SI):',au(SI()).efield(1.)
        print ' 1 au intensity (W/cm2):',au(SI()).intensity(1.)*1.e-4
        print '       1 au energy (eV):',au(eV()).energy(1.)
        print '       1 au energy (SI):',au(SI()).energy(1.)
        print 'wave length at 1au (SI):',au(SI()).length(2*pi/a_finestructure)
        print ' ----2 x 10e14, 800 nm, ----------------------------------------------'
        print '            field:', ((2.e18/au(SI()).intensity(1.))**0.5)
        print ' quiver amplitude:', ((2.e18/au(SI()).intensity(1.))**0.5/nm_inv().energy(1/800.)**2)
        print '             Pmax:', (2*10*((2.e18/au(SI()).intensity(1.))/nm_inv().energy(1/800.)**2))**0.5
        print '        omega(eV):', au(eV()).energy(nm_inv().energy(1/800.))
        print '  760nm omega(au):', nm_inv().energy(1/760.)
        print '         U_p (eV):', au(eV()).energy((SI().intensity(2.e18)/nm_inv().energy(1/800.)**2)/4.*10.)
        print '        2eV in nm:',1./eV(nm_inv()).energy(2.)
        exit(' --- units tests done ---')

class SI(Units):
    # conversion factors and constants in SI units
    def __init__(self,toUnits=None):
        self.name='SI'
        self.local=if_None_local(toUnits)
        self._length=1.
        self._mass=1.
        self._time=1.
        self._charge=1.
        self._mu0=4.*myPi*1.e-7
        # define the derived units
        self.supplement()

        # SI is special with "Volt" and "Watt"
        self._efield=1.
        self._intensity=1.

class au(Units):
    def __init__(self,toUnits=None):
        self.name='au'
        self.local=if_None_local(toUnits)
        self._length=bohr_radius
        self._mass=electron_mass
        self._time=bohr_radius/(speed_of_light*a_finestructure)
        self._charge=proton_charge
        self._mu0=4.*myPi*a_finestructure**2
        self.supplement() 

class eV(Units):
    def __init__(self,toUnits=None):
        self._name='eV'
        self.local=if_None_local(toUnits)
        self._energy=proton_charge
        self._frequency=self._energy/planck_constant

class nm_inv(Units):
    def __init__(self,toUnits=None):
        self._name='nm^-1'
        self.local=if_None_local(toUnits)
        self._frequency=1.e9*speed_of_light
        self._energy=planck_constant*self._frequency

class OptCyc(Units):
    def __init__(self,wavelength_SI,toUnits=None):
        self._name='OptCyc@'+str(wavelength_SI*1.e-9)+'(nm)'
        self.local=if_None_local(toUnits)
        self._time=wavelength_SI/speed_of_light

if __name__  == "__main__":
    """standalone: tests"""
    Units().test()
