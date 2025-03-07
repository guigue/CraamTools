#########################################################
#
# Goes : some scripts to deal with AIA data in Python
#        returns the filenames of raw data, and a AIA map.
#
#
# @Guiguesp - 2025-02-14BST19:41

__Version__ = '2025-02-14T19:41'

from sunpy import timeseries as ts
from sunpy.net import Fido
from sunpy.net import attrs as a
import sunpy.map

from matplotlib import pyplot as plt
from astropy import units as u

def get(isoStart, isoEnd, wavelength=1600):

    r, = Fido.search(a.Time(isoStart,isoEnd),
                     a.Instrument.aia,
                     a.Wavelength(wavelength*u.angstrom))

    ##
    #
    # How to get the data
    # >>> from sunpy.net import Fido
    # >>> data = Fido.fetch(r[i:j])
    #
    
    return r

def plot(data,savefig=False):

    mapa = sunpy.map.Map(data)
    
    fig = plt.figure()
    fig.set_size_inches(w=10,h=10)
    ax = fig.add_subplot(projection=mapa)
    mapa.plot(axes=ax)

    plt.show(block=False)

    return

    
    
