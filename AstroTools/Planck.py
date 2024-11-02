import numpy as np
from CraamTools.AstroTools import constants as ct

def bnu(T,nu):
    C1 = (2 * ct.Phys.Planck * nu**3 / ct.Phys.c**2 )
    C2 = ct.Phys.Planck *nu /(ct.Phys.Boltzmann * T)
    
    return C1 / (np.exp(C2)-1)

def bwl(T,wl):
    C1 = 2*ct.Phys.Planck*ct.Phys.c**2/wl**5
    C2 = ct.Phys.Planck*ct.Phys.c/(ct.Phys.Boltzmann*T*wl)
    return C1 / (np.exp(C2)-1)

def nuM(T):
    x = 2.82144
    C = x*ct.Phys.Boltzmann/ct.Phys.Planck
    return C * T

def RJ(T,nu):
    return 2 * ct.Phys.Boltzmann * T / ct.Phys.c**2 * nu**2

def Wien(T,nu):
    C1 = (2 * ct.Phys.Planck * nu**3 / ct.Phys.c**2 )
    C2 = ct.Phys.Planck * nu /(ct.Phys.Boltzmann * T)
    return C1 * np.exp(-C2)



