import os, sys, glob
import numpy as np
import datetime as dt
import pandas as pd
from astropy.io import fits
from astropy.time import Time as astrotime

tablename = pd.DataFrame( {'extension' : ['lis'     ,'k7o'          ,'phf'    ,'apl'      ],
                           'station'   : ['San Vito','Sagamore-hill','Palehua','Learmonth']} )
tablename.set_index('extension', inplace=True)

__Version__ = '2024-11-22T1700BST'

def version():
    return __Version__

def read(fname,verbose=True):

    extension = fname.split(sep='.')[-1]
    station = tablename.loc[extension]['station']
    
    d      = []
    time   = []
    mstime = []
    jd     = []
    atime  = []

    nlines=0
    nrlines=0
    with open(fname,'r') as fh:
        for line in fh:
            c=line.strip().split()
            nlines+=1
            if len(c) == 9 :
                nrlines+=1
                year = int(c[0][4:8])
                mon  = int(c[0][8:10])
                day  = int(c[0][10:12])
                hour = int(c[0][12:14])
                mn   = int(c[0][14:16])
                sec  = int(c[0][16:18])
                dtime= dt.datetime(year,mon,day,hour,mn,sec,0)
                ms   = hour*3600000+mn*60000+sec*1000
                t    = astrotime(dtime,format='datetime',scale='utc')

                flux = [float(c[1]),float(c[2]),float(c[3]),float(c[4]),float(c[5]),float(c[6]),float(c[7]),float(c[8])]
                d.append(flux)
                time.append(dtime)
                mstime.append(ms)
                jd.append(t.jd1+t.jd2)
                atime.append(t)

    fh.closed
    if verbose:
        print(' ')
        print('---- Statistics ----')
        print('File :  '+fname+' ,  Station :  '+station)
        print('Number of Records in file = {0:5d}'.format(nlines))
        print('Number of Records read    = {0:5d}'.format(nrlines))
        print('Number of Records missing = {0:5d}'.format(nlines-nrlines))
        print('  ')

    dtime = np.asarray(time)
    ms    = np.asarray(mstime)
    jd    = np.asarray(jd)
    flux  = np.asarray(d)

    freq  = np.array([0.245,0.410,0.610,1.415,2.695,4.995,8.8,15.4])

    data = {}
    data.update({'frequency':freq})
    data.update({'time':dtime})
    data.update({'jd':jd})
    data.update({'mstime':ms})
    data.update({'flux':flux})
    data.update({'station':station})
    data.update({'date_obs':atime[0].fits[0:10]})
    data.update({'t_start':atime[0].fits})
    data.update({'t_end':atime[-1].fits})

    return data

def writeFits(d,fitsname,*comments):

    hdu = fits.PrimaryHDU()
    hdu.header.append(('origin','NOAA/NGDC','ftp://ftp.ngdc.noaa.gov/STP/space-weather/solar-data'))
    hdu.header.append(('telescop','RSTN',''))
    hdu.header.append(('observat',d['station'],''))

    hdu.header.append(('date-obs',d['date_obs'],''))
    hdu.header.append(('t_start',d['t_start'],''))
    hdu.header.append(('t_end',d['t_end'],''))
    hdu.header.append(('data_typ','Flux density','SFU'))

    hdu.header.append(('history','Created by CraamTools.AstroTools.rstn.writeFits, version '+getVersion()))
    hdu.header.append(('history','Created on '+ str(dt.datetime.today())))
    hdu.header.append(('history','Created by '+ os.getlogin() + ' on ' + os.uname()[1]))


    nfreq = str(d['frequency'].shape[0])

    coldefs = fits.ColDefs([fits.Column(name   = 'frequency',
                            format = 'E'   ,
                            unit   = 'GHz'  ,
                            bscale = 1.0   ,
                            bzero  = 0.0   ,
                                        array  = d['frequency'])])
    tbfreq   = fits.BinTableHDU.from_columns(coldefs)

    fits_cols = []
    col       = fits.Column(name   = 'ms',
                            format = 'J'   ,
                            unit   = 'ms'  ,
                            bscale = 1.0   ,
                            bzero  = 0.0   ,
                            array  = d['mstime'])

    fits_cols.append(col)

    col       = fits.Column(name   = 'JulDay',
                            format = 'E'   ,
                            unit   = 'day'  ,
                            bscale = 1.0   ,
                            bzero  = 0.0   ,
                            array  = d['jd'])

    fits_cols.append(col)

    col       = fits.Column(name   = 'flux',
                            format = nfreq+'E'   ,
                            unit   = 'GHz'  ,
                            bscale = 1.0   ,
                            bzero  = 0.0   ,
                            array  = d['flux'])

    fits_cols.append(col)

    coldefs = fits.ColDefs(fits_cols)
    tbdata   = fits.BinTableHDU.from_columns(coldefs)
    # Comments and other headers
    tbdata.header.append(('comment','Time is in milliseconds since 0 UT',''))
    for c in comments:
        tbdata.header.append(('comment',c,''))

    hduList = fits.HDUList([hdu,tbfreq,tbdata])
    hduList.writeto(fitsname,overwrite=True)

    return
