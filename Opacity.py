import numpy as np
import datetime as dt
from astropy import units as u
from astropy import constants as c
from scipy.optimize import curve_fit

import pdb

########################################################################
#
# Opacity: A set of functions related to the atmospheric optical depth
#
# Usage:  from astropy import units as u
#         print(Opacity.pwv(temperature=10*u.Celsius, humidity=10*u.Unit(''), Hh20=1.5*u.km))
#
# Author: @guiguesp - 2021-05-10 in Frozen São Paulo
#                     2021-07-31 Revised formula for Dew point temperature
#                                Added Units
#                                Checked accuracy
#                     2021-10-28 Added refraction formula for radio waves (SST formula)
#                     2022-04-14 Added Extinction Opeacity Calculation
#
########################################################################

def refraction(ele,temp=10,humi=10,pres=570):
    
    ctok = 273.15
    MM_HG_2_HPa=1.333333
    PMIN = 1.00E-6  
    ZMAX = 1.55
    
    z    = np.radians(90 - ele) < ZMAX
    tk   = temp+ctok
    rh   = humi/100.0
    pmb  = pres * MM_HG_2_HPa
    if ( pmb > PMIN ): 
        ps = 10**(((0.03477 * tk - 8.71170) / (0.00412 * tk - 0.12540)) * (1.0 + pmb * 4.5E-06))
        pw = rh * ps / ( 1 - ( 1 - rh ) * ps / pmb ) 
        b  = 4.4474E-06 * tk * ( 1.0 - 0.0074 * pw ) 
        s = np.sin (z)
        c = np.cos (z)
        return  - ((77.6890E-06 * pmb - ( 6.3938E-06 - 0.375463 / tk ) * pw ) / tk ) * \
            ( 1 - b ) * s / np.sqrt ( c * c + 0.001908 + 0.6996 * b - 0.00003117 * s / c ) * (180 / np.pi)
    else:
        return 0
        

def dew(T,r):
    ############################
    #
    # Dew Point Temperature
    #
    #     Inputs:  T, dry-bulb temperature (Celsius)
    #              r, relative humidity in fraction (0...1)
    #
    #     Output:  Dew point temperature in Celsius
    #
    #
    #     Source: https://en.wikipedia.org/wiki/Dew_point (2021-07-31)
    #
    ####################################################################
    
    a = 6.1121 * u.mbar
    b = 18.678 * u.Unit('')
    c = 257.14 * u.Celsius
    d = 234.50 * u.Celsius

    Psm = a * np.exp( (b-T/d) * (T/(c+T)) )
    gm  = np.log(r * np.exp( (b-T/d) * (T/(c+T))) )
    return   c * gm / (b - gm)

def P0(T,h):
    ##############################
    #
    # Water Vapor Partial Pressure
    #
    #    Inputs:  T, dry-bulb temperature in Celsius
    #             h, relative humidity in %
    #
    #    Output:  Water Vapor partial pressure in mbar
    #
    #    Source: ALMA memo 237 - Butler, B (1998)
    #
    ######################################################
    
    D    = dew(T, h/100 )
    expo = 1.81E+00 + (1.727E+01 * D /(D + 237.3 * u.Celsius))
    return  np.exp(expo) * u.mbar
    
def pwv(temperature=20*u.Celsius, humidity=20*u.Unit(''), Hh2o=2.0*u.km):

    #################################
    #
    # Precipitable Water Vapor
    #
    #    Parameter inputs:
    #           temperature = dry-bulb temperature in Celsius
    #           humidity    = relative humidity in %
    #           Hh2o        = Water scale height in km (typically 1.5 - 2 )
    #
    #    Output: Precipitable Water Vapor in mm
    #
    #    Source: ALMA memo 237 - Butler, B (1998)
    #
    ########################################################
    
    amu      = 1.66053906660e-24 * u.g  # Atomic Mass Unit in g
    mw       = 18 * amu                 # Water mole mass
    rhol     = 1 * u.g / u.cm**3        # water density
    constant = mw / ( rhol * c.k_B.cgs)
    return (constant * P0(temperature,humidity).cgs * Hh2o.to('cm') / temperature.to(u.Kelvin, equivalencies=u.temperature())).to('mm')

def Extinction_function(x,*a):

    return a[0] * np.exp(-a[1]/np.sin(np.radians(x)))

def Extinction(x,y):

    est = [100,0.1]
    par, cov = curve_fit(Extinction_function,x,y,p0=est)
    yfit     = Extinction_function(x,*par)
    
    return par,yfit,cov

    
