from read_input import *
from math import *
from my_constants import *

class PhysicalSystem:
    """mass, potential, interaction, and time-dependence of the interaction"""

    def __init__(self,name,mass,potential,timedep):
        self.mass=mass
        self.potential=potential
        self.timdep=timedep

    @classmethod
    def read(cls):
        """everything that defines the physical system"""    
        mass=read_input(1.,'mass',doc='mass of the system')
        potential=Potential.read()
        timdep=LaserField.read()
        return PhysicalSystem(name,mass,potential,timedep)

def parameter(name):

    # read parameters when needed

    try: mass
    except NameError: mass=read_input(1.,'.mass',doc='mass(es) in the system')

    try: kappa
    except NameError: kappa=read_input(1.,'.kappa',doc='kappa: harmonic constant')

    try: langl
    except NameError:  langl=read_input(10,'.angular momentum',doc='angular momentum')

    try: peak_field
    except NameError:  peak_field=read_input(0.,'.laser',1,doc='maximal field strength')

    try: frequency
    except NameError: frequency=read_input(1.,'.laser',2,doc='laser frequency')

    def _mass(t): return mass
    def _plus_mass_half(t): return mass/2.
    def _minus_mass_half(t): return -mass/2.
    def _kappa(t): return kappa
    def _kappa_half(t): return kappa/2.
    def _plus_one(t): return 1.
    def _minus_one(t): return -1.
    def _lsq_mass_half(t): return 0.5/mass*langl*(langl+1)
    def _laserpulse(t): return peak_field*cos(2*myPi*frequency*t)

    # assign according to name
    if   name.strip() == 'mass':      return _mass
    elif name.strip() == 'kappa':     return _kappa
    elif name.strip() == 'kappa/2':   return _kappa_half
    elif name.strip() == '-1':        return _minus_one
    elif name.strip() == '+1':        return _plus_one
    elif name.strip() == '1':         return _plus_one
    elif name.strip() == '-1/2m':     return _minus_mass_half
    elif name.strip() == '1/2m':      return _plus_mass_half
    elif name.strip() == '+l(l+1)/2m':return _lsq_mass_half
    elif name.strip() == 'laser(t)':  return _laserpulse
    else: sys.exit('unknown parameter name "'+name.strip()+'"')
