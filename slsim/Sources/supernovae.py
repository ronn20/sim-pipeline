import astropy.units as u
import numpy as np
from scipy import integrate
from astropy.cosmology import Planck13, z_at_value
import math
import scipy.constants
import astropy.cosmology
import astropy.cosmology.units
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import WMAP9 as cosmo

"DTD as a function of x=time where input is either single value or a numPy array"
def delaytimedistribution(x):
    if isinstance(x, (int, float)):
        "checks if input is a single value"
        if x<0.04:
            return 0
        elif x<0.5:
            return 12.2196
        else:
            return 1.02
    elif isinstance(x, np.ndarray):
        "checks if input is an array"
        result = np.where(x < 0.04, 0, np.where(x < 0.5, 12.2196, 1.02))
        return result
    else:
        raise ValueError("input must be single value or numpy array")

"describes volumetric rate of supernova as a function of time in Gyrs"
def volumetricSNrate(t):
    "takes in some numpy array of times representing time since beginning of universe"
    t = np.atleast_1d(t)
    "checks that the input as at least 1 computable value"
    result = np.empty_like(t)
    "creates an empty array with the same shape as our input array"
    val_min = 0
    "sets min value for integration which changes for every iteration of the loop"
    for i, val in enumerate(t):
        "loops run through array of T values to calculate each StarFormationHist and then integral at that T."
        def StarFormationHist(ti):
            "describes the cosmic star formation rates based on UV and IR band observations"
            z = z_at_value(Planck13.age, ti * u.Gyr,zmax=10000.0)
            "correspondance between time and redshift using z_at_value"
            s = 0.15 * (1 + z)**2.7 / (1 + ((1 + z) / 2.9)**5.6)
            "equation of star formation history as a function of redshift"
            return s
        if isinstance(val, (int, float)):
            "checks if input is a single value"
            integral_result, _ = integrate.quad(lambda ti: StarFormationHist(ti) * delaytimedistribution(val), val_min, val)
            result[i] = integral_result
            "integrates and then appends at current i index that result to the array intialized outside the loop"
        elif isinstance(val, np.ndarray):
            "checks if input is an array of value"
            integral_result, _ = integrate.quad(lambda ti: StarFormationHist(ti) * delaytimedistribution(val), val_min, val)
            result[i] = integral_result
            "integrates and then appends at current i index that result to the array intialized outside the loop"
        val_min = val
        "sets new bounds by taking in the old max val as the new min value in order to integrate by parts"
        cumulativeresult=np.cumsum(result)
        "sums up all the array 'parts' at that index i"
    return cumulativeresult                   
 
def supernovatotal(z):
    "main SNe equation as a function of redshift"
    pi=math.pi
    c=scipy.constants.speed_of_light
    "defining constants pi and speed of light for the main # of SNe equation"
    cosmo=FlatLambdaCDM(H0=73, Om0=0.3)
    "FlatLambdaCDM model adopted directly from article"
    H0=cosmo.H(z)
    "The hubble paramater at some redshift"
    d=cosmo.angular_diameter_distance(z)
    "d is equal to D_a which is the angular diameter distance at some redshift z"
    E=H0/cosmo.H0
    "E is the dimensionless Hubble parameter based on hubble constant of 73 km s^-1 Mpc^-1"
    SuperNovaRate=(((4*pi)*(c/H0))*(((1+z)**2*(d**2))/E)*(1/(1+z)))*volumetricSNrate(z)
    "Total Number of Supernova as a function of redshift"
    return SuperNovaRate
