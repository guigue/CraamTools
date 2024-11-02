import numpy as np
from astropy.time import Time as astrotime
import pdb

##################################################
#
# From Jean Meeus' Astronomical Algorithms
#
#------------------------------------------------
#
# Guiguesp - 2023-03-10T22:55
#
##################################################

def Gregorian(ojd):

    JD = ojd + 0.5

    Z = int(np.floor(JD))
    F = JD % 1

    if (Z < 2299161):
        A = Z
    else:
        alpha = int((Z-1867216.25) / 36524.25)
        A = Z + 1 + alpha - int(alpha/4)
    B = A + 1524
    C = int((B-122.1)/365.25)
    D = int(365.25*C)
    E = int((B-D)/30.6001)

    day = B - D - int(30.6001*E)

    if (E<14):
        month = E - 1
    else:
        month = E -13

    if (month > 2):
        year = C - 4716
    else:
        year = C - 4715

    hh = F*24
    hour = int(np.floor(hh))

    mm = (hh-hour)*60
    minute = int(np.floor(mm))

    second = (mm-minute)*60
    
        
    return [year,month,day,hour,minute,second]


def Julian(y,M,d,h,m,s):

    if (M < 1) or (M>12) or (y < 0) or (d<0) or (d>31) or (h <0) or (m<0) or (s<0):
        return []
    
    df = d + (h+m/60+s/3600)/24
    
    if (M <3):
        y = y -1
        M = M+12

    A = y // 100
    B = 2 - A + A//4

    return int(365.25* (y + 4716)) + int(30.6001 * (M+1)) + df + B - 1524.5
        
