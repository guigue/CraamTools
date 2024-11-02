import numpy as np

_version_ = '20210201T1030'

def eta_s(epsilon,wl,unit='mm'):
    # epsilon: the surface RMS --> in micrometers <--
    # wl:      the wavelength. Default unit is mm,
    # unit:    units of the wavelength. Default is mm
    #          Also accepted cm and um.
    #
    # Author: @guiguesp 2021-02-01 Sampa (at home)

    if (unit=='mm'):
        wl *= 1000

    if (unit=='cm'):
        wl *= 10000

    return np.exp(-(4 * np.pi * epsilon / wl)**2)

def Version():
    return _version_
