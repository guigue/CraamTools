# External methods
import sys, string, os, struct, glob
import numpy as np
import datetime as dt
from astropy.io import fits
import pdb

##################################
Version = '20200527T1136BRT'     #
##################################

def ms2dt(base,ms):
    bdate = base.split('-')
    year = int(bdate[0])
    month = int(bdate[1])
    day = int(bdate[2])

    hours        =  ms // 3600000
    minutes      = (ms % 3600000) // 60000
    seconds      = ((ms % 3600000) % 60000) / 1.0E+03
    seconds_frac = seconds - seconds.astype('int')
    useconds     = int(seconds_frac * 1e6)

    return  dt.datetime(year,month,day,hours,minutes,seconds.astype('int'),useconds)

    
def ms2iso(base,ms):

    t_dt = ms2dt(base,ms)
    return t_dt.isoformat()
