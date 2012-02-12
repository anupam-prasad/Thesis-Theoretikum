# the laser class
class Laser:
    """complete definition of a laser pulse
    shape    ... pulse shape
    intensity... in W/cm2
    FWHM    .... full width half maximum of intensity
    omega   .... carrier frequency(a.u.)
    """
    def __init__(self,shape,intensity,FWHM,omega):
        self.shape=shape
        self.intensity=intensity
        self.FWHM=FWHM
        self.omega=omega

    def laser(shape,intensity,FWHM,omega):
        currentLaser=Laser(shape,intensity,FWHM,omega)
