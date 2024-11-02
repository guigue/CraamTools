# External methods
import sys, string, os, struct, glob
import numpy as np
import datetime as dt
from astropy.io import fits
import pdb

##################################
Version = '20200527T1136BRT'     #
##################################

def ms2dt(base,ms, sst=False):
    bdate = base.split('-')
    year = int(bdate[0])
    month = int(bdate[1])
    day = int(bdate[2])

    if sst:
        ms = ms // 10
        
    hours        =  ms // 3600000
    minutes      = (ms % 3600000) // 60000
    seconds      = ((ms % 3600000) % 60000) / 1.0E+03
    seconds_int  = seconds.astype(int)
    seconds_frac = seconds - seconds_int
    useconds     = (seconds_frac * 1.0E+06).astype('int')
    dttime       = np.array(np.empty(ms.shape[0]),dtype=dt.datetime)

    for i in np.arange(ms.shape[0]):
        dttime[i] = dt.datetime(year,month,day,hours[i],minutes[i],       seconds_int[i],useconds[i])

    return  dttime
