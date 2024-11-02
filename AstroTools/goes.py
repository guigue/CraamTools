import os, sys, glob
import numpy as np
import datetime as dt
from astropy.io import fits
from astropy.time import Time as astrotime


def read(fname,date_obs):
    g=fits.open(fname)
    
    year  = int(date_obs[0:4])
    month = int(date_obs[5:7])
    day   = int(date_obs[8:10])

    tlist = []
    atime = []
    for t in g[1].data['TIME']:
        hour = int(t/3.6e3)
        mns  = int(t/6.0e1) % 60
        sec  = int(t) % 60
        msec = int((t-int(t))*1.0e6)
        time = dt.datetime(year,month,day,hour,mns,sec,msec)
        tlist.append(time)
        atime.append(astrotime(time,format='datetime',scale='utc'))
    

    goes = {}
    goes.update({'time':g[1].data['TIME']})
    goes.update({'dtime':np.asarray(tlist)})
    goes.update({'atime':np.asarray(atime)})
    goes.update({'flow':g[1].data['FLOW']})
    goes.update({'fhigh':g[1].data['FHIGH']})
    goes.update({'dflow':g[1].data['DFLOW']})
    goes.update({'dfhigh':g[1].data['DFHIGH']})
    goes.update({'emis':g[1].data['EMIS']})
    goes.update({'lrad':g[1].data['LRAD']})
    goes.update({'lxlow':g[1].data['LXLOW']})
    goes.update({'lxhigh':g[1].data['LXHIGH']})
    g.close()
    return goes

