import numpy as np

def angle(ele,temp=20,humi=30,press=590):

    ctok        = 273.15
    MM_HG_2_HPa = 4/3
    PMIN        = 1.00E-6  
    ZMAX        = 1.55

    z    = np.radians(90.0 - ele)
    if (z >= ZMAX):
        z = ZMAX
        
    tk   = temp+ctok
    rh   = humi/100.0
    pmb  = press * MM_HG_2_HPa
    if ( pmb >= PMIN ):

        ps = 10.0E0**(((0.03477 * tk - 8.71170) /  (0.00412 * tk - 0.12540)) * (1.0 + pmb * 4.5E-6))
        pw = rh * ps / ( 1.0 - ( 1.0 - rh ) * ps / pmb ) 
        b  = 4.4474E-6 * tk * ( 1.0 - 0.0074 * pw ) 
        s  = np.sin ( z )
        c  = np.cos ( z )
        
        refr = -np.degrees((77.6890e-6 * pmb - ( 6.3938e-6 - 0.375463 / tk ) * pw ) / tk ) \
            * ( 1.0 - b ) * s / np.sqrt ( c * c + 0.001908 + 0.6996 * b - 0.00003117 * s / c ) 
    else:        
        refr = 0.0

    return refr

    

