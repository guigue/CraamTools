import numpy as np
from astropy import constants as c

def X(T,Q,En):

    """
    T = temperature [K]
    Q = Emission measure [cm-3]
    En = Photon energy [keV]

    Output:
       X-ray flux in Ph/cm2/s/keV

    Reference:
       Tandberg-Hansen & Emslie

    Author:
        Guiguesp @ 2018-05-23 from an old IDL pro
    """

    Z = 1.4                              # Allen (1973)
    D = 5.7E-12 * Z**2
    keV = 1.0E+03 * c.e.value * 1.0E+07  # Kev to ergs

    return D * Q * np.exp(-En * keV / (c.k_B.cgs.value * T) ) / (En*keV * np.sqrt(T) ) / (4 * np.pi * c.au.cgs.value**2) * keV

def Radio(Nmed   = 1.0E+10,
          Tmed   = 1.0E+08,
          Emis   = -1     ,
          Height = 1.0E+09,
          Size   = 1.0E+09,
          quiet  = True)  :

    """
    Radio
    Computes the bremsstrahlung emission of an isothermal source at 1 AU
    in the Rayleigh Jeans approximation.

    INPUT
    none

    OUTPUT
    freq : float array with frequencies (GHz)
    flux : float array with fluxes (s.f.u.)
    tau  :

    PARAMETERS

    nmed    medium  density [cm-3]
    tmed   medium temperature [K]
    Height  height [cm]
    size    size [arc sec]

    HISTORY
    Written by Guigue (guigue@craam.mackenzie.br) using Dulk (1985) formulas
    Buenos Aires, Sept 26, 2006 (1st RIARCHE)

    Ported to python in April 14, 2021 (still not covid-vaccinated)

    """
    erg2sfu = 1.00E+19  # ergs/cm2 -> s.f.u.

    A     = np.pi * Size**2 / 4
    omega = A / c.au.cgs.value**2
    freq  = 1.0E+01**(np.linspace(0,3,1000)) * 1.0E+09
    snu   = 2 * c.k_B.cgs.value * Tmed * freq**2 / c.c.cgs.value**2

    if (Emis > 0 ):
        if (Tmed < 2.0E+05):
            tau = (0.0115 * Emis) / (A * freq**2 * Tmed**1.5) *  (18.2 + 1.5*np.log(Tmed) - np.log(freq))
            flux = snu * omega * erg2sfu * (1 - np.exp(-tau))
        else:
            tau = (0.0115 * Emis ) / (A * freq**2 * Tmed**1.5) * (24.5 + np.log(Tmed) - np.log(freq))
            flux = snu * omega * erg2sfu * (1 - np.exp(-tau) )
    else:
         if (Tmed < 2.0E+05):
            tau  = (0.0115 * Nmed**2 * Height) / (freq**2 * Tmed**1.5) * (18.2 + 1.5 * np.log(Tmed) -np.log(freq))
            flux = snu * omega * erg2sfu * (1 - np.exp(-tau))
         else:
            tau = (0.0115 * Nmed**2 * Height) / (freq**2 * Tmed**1.5) * (24.5 + np.log(Tmed) -np.log(freq))
            flux = snu * omega * erg2sfu * (1 - np.exp(-tau))
    if (not quiet):
        print("\n\nParameters used:\n ")
        print("T plasma = {0:4.2e} K\n".format(Tmed))
        print("Source Size = {0:4.2e} cm, height = {1:4.2e}\n".format(Size,Height))
        if (Emis<0):
            print("Emission Measure = {0:4.2e} cm-3\n".format(Nmed**2* A*Height))
            print("Source Density = {0:4.2e}\n\n".format(Nmed))
        else:
            print("Source Density = {0:4.2e} cm-3\n".format(np.sqrt(Emis/(A*Height))))
            print("Emission Measure = {0:4.2e} cm-3\n\n".format(Emis))

    freq = freq / 1.0E+09

    return freq,flux
