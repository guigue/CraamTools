import numpy as np
from astropy import units as u
from astropy import constants as c
import pdb

def flux(KT,Q,en):

    # Compute Bremsstrahlung optically thin emission using output from OSPEX
    # KT : plasma temperature in keV
    # Q  : Emission Measure in 1.0E+49 cm-3
    # en : photon energy limits in keV, ndarray with two elements 

    if en.shape[0] < 2 :
        return
    
    Z   = 1.4  
    D   = 5.7E-12 * Z**2
    EM  = Q * 1.0E+49
    keV = 1000 * c.e.value 
    
    hf  = 10**np.linspace(np.log10(en[0]),np.log10(en[1]),100)
    
    T   = KT * keV / c.k_B.value

    return {'en':hf, 'fx':D * EM * np.exp(-hf/KT) / (hf * keV * u.J.to(u.erg) * np.sqrt(T)) / (4 * np.pi * c.au.to(u.cm).value**2) * keV * u.J.to(u.erg)}
    
