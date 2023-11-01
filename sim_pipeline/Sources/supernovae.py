from sim_pipeline.Sources.source_base import SourceBase
class SuperNova(SourceBase):
    def __init__(self, cosmo, sky_area):
        """
        :param cosmo: cosmology
        :type cosmo:   ~astropy.cosmology class
        :param sky_area: Sky area over which galaxies are sampled. Must be in units of
            solid angle.
        :type sky_area: `~astropy.units.Quantity`
        """
        super(Supernova, self).__init__(cosmo=cosmo, sky_area=sky_area)
import astropy.units as u

"DTD as a function of x=time where input is either single value or a numPy array"
def delaytimedistribution(x):
    import numpy as np
    import astropy.units as u
    if isinstance(x, (int, float)):
        if x<0.04:
            return 0*u.Gyr
        elif x<0.5:
            return 12.2196*u.Gyr
        else:
            return 1.02*u.Gyr
    elif isinstance(x, np.ndarray):
        result = np.where(x < 0.04, 0, np.where(x < 0.5, 12.2196, 1.02))
        return result * u.Gyr
    else:
        raise ValueError("input must be single value or numpy array")
                    
        
def volumetricSNrate(t):
    "correspondance between time and redshift using z_at_value"
    from astropy.cosmology import Planck13, z_at_value
    import astropy.units as u
    z=z_at_value(Planck13.age, t*u.Gyr)
    
    "SNe volumetric rate integral"
    StarFormationHist=.15*(1+z)**2.7/(1+((1+z)/2.9)**5.6)
    "describes the cosmic star formation rates based on UV and IR band observations"
    " below is integration to obtain SNe volumetric rate at some redshift"
    import scipy.integrate as integrate
    SNVolumetricIntegral= lambda x: StarFormationHist*delaytimedistribution(x)/u.Gyr
    result,_=integrate.quad(SNVolumetricIntegral,0,t)
    return result*u.Unit('Gyr')
 
def supernovatotal(z):
    "main SNe equation as a function of redshift"
    import math
    pi=math.pi
    import scipy.constants
    c=scipy.constants.speed_of_light
    "defining constants pi and speed of light for the main # of SNe equation"
    import astropy.cosmology
    import astropy.cosmology.units
    from astropy.cosmology import FlatLambdaCDM
    cosmo=FlatLambdaCDM(H0=73, Om0=0.3)
    "FlatLambdaCDM model adopted directly from article"
    from astropy.cosmology import WMAP9 as cosmo
    H0=cosmo.H(z)
    "The hubble paramater at some redshift"
    d=cosmo.angular_diameter_distance(z)
    "d is equal to D_a which is the angular diameter distance at some redshift z"
    E=H0/cosmo.H0
    "E is the dimensionless Hubble parameter based on hubble constant of 73 km s^-1 Mpc^-1"
    SuperNovaRate=(((4*pi)*(c/H0)*(1+z)**2*(d**2))/(E))*volumetricSNrate(z)
    "Total Number of Supernova as a function of redshift"
    return SuperNovaRate
