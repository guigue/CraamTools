import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd
from astropy import constants as c

def read(fname):

    wn = []
    wl = []
    nu = []
    t  = []

    d = pd.read_csv(fname,sep='\t',header=2)

    wn = np.asarray(d.wave_number)              # wavenumber [cm-1]
    wl = np.asarray(1.0E+04 / wn )              # wavelength [um]
    nu = np.asarray(c.c.value * wn * 1.0E-10)   # Frequency [THz]
    t  = np.asarray(d.Trans)                    # Transmission

    return {'Wavenumber':wn, 'Wavelength':wl, 'Frequency':nu, 'Transmission':t}

def plot(spec):

    fig, ax1 = plt.subplots(1, 1)
    ax1.set_xlabel(r'$\nu$ [THz]')
    ax1.set_ylabel('Transmission')
    ax1.semilogy(spec['Frequency'],spec['Transmission'],color='red')
    plt.show()

    return
