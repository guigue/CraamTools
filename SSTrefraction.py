import numpy as np

#############################
#
# pnt_refract: From sst software.
#
# Input:
#       ele: telescope elevation in degrees
#       tk:  temperature in Celsius degrees
#       pmb: pressure in HPa
#       rh:  humidity in %
#
# Returns:
#       elevation correction in degrees
#
##########################################
#
#  Guiguesp - 2026-03-09
#
########################################

_version_ = '2026-03-09T14:45BRT'

def pnt_refract(ele,tk,pmb,rh):

    tk   += 275.15
    rh   /= 100 
    ztop =  np.radians(90 - ele)
  
    if ( pmb < 1.0E-06 ):
        return 0
   
    if (ztop < 1.55):
        z = ztop
    else:
        ztop = 1.55

    expo = ((0.03477 * tk - 8.71170 ) / ( 0.00412 * tk - 0.12540 ))  * ( 1.0 + pmb * 4.5E-06 )
    ps = 10**expo
    pw = rh * ps / ( 1 - ( 1 - rh ) * ps / pmb )
    b = 4.4474E-06 * tk * ( 1 - 0.0074 * pw )
    s = np.sin ( z )
    c = np.cos ( z )

    dRef =  np.degrees( ( 77.6890E-06 * pmb - ( 6.3938E-06 - 0.375463 / tk ) * pw )/ tk ) * ( 1.0 - b ) * s / np.sqrt ( c**2 + 0.001908 + 0.6996 * b - 3.117E-05 * s / c )
    
    return dRef
