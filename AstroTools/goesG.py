import os, sys, glob
import numpy as np
import pandas as pd
import datetime as dt
from astropy.io import fits
from astropy.time import Time as astrotime
import matplotlib.pyplot as plt
from astropy import constants as c
import pdb

y = 1900
m = 1
d = 1


def correctJD(ojd):

    global y,m,d
    
    JD = ojd + 0.5

    Z = int(np.floor(JD))
    F = JD % 1

    if (Z < 2299161):
        A = Z
    else:
        alpha = int((Z-1867216.25) / 36524.25)
        A = Z + 1 + alpha - int(alpha/4)
    B = A + 1524
    C = int((B-122.1)/365.25)
    D = int(365.25*C)
    E = int((B-D)/30.6001)

    d = B - D - int(30.6001*E)

    if (E<14):
        m = E - 1
    else:
        m = E -13

    if (m > 2):
        y = C - 4716
    else:
        y = C - 4715

    if (y < 1981):
        y+=1900
        
    return Julian(F*3600*24)

def Julian(time):

    f_of_day = time/3600/24
        
    # Gregorian adopted in Oct., 15, 1582
    greg = 15 + 31 * (10 + 12 * 1582) ;
            
    if(m > 2):
        jy = y
        jm = m + 1 
    else:
        jy = y - 1 
        jm = m + 13 

    jd = int(365.25 * jy) + int(30.6001 * jm) + d + 1720995

    if ( d + (31 * (m + 12 * y )) >= greg ) :
        ja = int(0.01 * jy)
        jd = jd + 2 -ja + int(0.25 * ja)

    return float(jd)+f_of_day-0.5 

def getJD(t):

    hh = int(t/3600)
    mm = int((t % 3600)/60)
    ss = int( ( ( (t % 3600)/60) % 1) * 60)
    um = round( ( ( ( ( (t % 3600)/60) % 1) * 60) % 1) * 1.0E+06)
    return astrotime(dt.datetime(y,m,d,hh,mm,ss,um)).jd


def read(lista,version='1981'):

    global y,m,d
    
    go = np.empty((0),dtype=[('jd','float'),('FLow','float'),('FHigh','float')])
    
    for arch in lista:

        g=fits.open(arch)
        print(arch)
        y = int(g[0].header['date-obs'][6:])
        if y < 1900:
            y+=1900
        m = int(g[0].header['date-obs'][3:5])
        d = int(g[0].header['date-obs'][0:2])

        if version=='1981' :
            try:
                t = g[0].data[0,:]
            except:
                print("\n\n\n archivo = {0:s}\n\n\n".format(arch))
            Nrec = g[0].data[0,:].shape[0]
        else:
            t  = g[2].data['time'][0]
            Nrec = g[2].data['time'].shape[1]
        jd = list(map(Julian,t))
        
        __temp__ = np.zeros(Nrec,dtype=[('jd','float'),('FLow','float'),('FHigh','float')])
        __temp__['jd'] = np.asarray(jd)
        if version == '1981':
            __temp__['FLow'] = g[0].data[1,:]
            __temp__['FHigh'] = g[0].data[2,:]
        else:
            __temp__['FLow'] = g[2].data['flux'][0,:,0]
            __temp__['FHigh'] = g[2].data['flux'][0,:,1]

        g.close()
        go = np.concatenate((go,__temp__))

    return go

def aperiod(dirn='1981',satn='14',version='1981',startd='0101',endd='1231'):

    os.chdir(dirn)
    
    wname = 'go' + satn + '*fits'
    listaP = glob.glob(wname)
    listaP.sort()
    
    lista = list()
    if int(dirn) < 1999:
        startfn = 'go' + satn + dirn[2:] + startd + '.fits'
        endfn   = 'go' + satn + dirn[2:] + endd   + '.fits'
    else:
        startfn = 'go' + satn + dirn + startd + '.fits'
        endfn   = 'go' + satn + dirn + endd   + '.fits'
        
    for f in listaP:
        if (f >= startfn) and (f <= endfn):
            lista.append(f)
    lista.sort()
    go = read(lista,version=version)

    os.chdir('..')

    return go
    

def oneyear(dirn='1981',satn='14',version='1981'):

    os.chdir(dirn)
    wname = 'go' + satn + '*fits'
    lista = glob.glob(wname)
    lista.sort()
    go = read(lista,version=version)
    os.chdir('..')

    return go
    
def oneminute_integrated(go):

    # Don't work !!!

    djd = 1 / 1440
    NR = go.shape[0]
    ii = int(0)

    jd = []
    fl = []
    fh = []
    
    while ii < NR:
        jd0 = go['jd'][ii]
        ig,=np.where(go['jd'] >= jd0+djd)
        ie = ig[0]-1

        jd.append((jd0+go['jd'][ie])/2)
        fl.append(go['FLow'][ii:ie].mean())
        fh.append(go['FHigh'][ii:ie].mean())

        ii=ie+1

    igo = np.zeros(len(jd),dtype=[('jd','float'),('FLow','float'),('FHigh','float')])
    igo['jd'] = np.asarray(jd)
    igo['FLow']= np.asarray(fl)
    igo['FHigh']=np.asarray(fh)

    return igo
                  
def oneminute_interp(go):

    Tmin = round((go['jd'][-1]-go['jd'][0])*60*24)
    jd = np.linspace(go['jd'][0],go['jd'][-1],Tmin)
    igo = np.zeros(Tmin,dtype=[('jd','float'),('FLow','float'),('FHigh','float')])
    igo['jd']    = jd
    igo['FLow']  = np.interp(jd,go['jd'],go['FLow'])
    igo['FHigh'] = np.interp(jd,go['jd'],go['FHigh'])
    
    return igo

def all():

    for year in np.arange(1981,2021):
        Ngo=np.load('int_go'+str(year)+'.npy')
        if (year <= 1999):
            print('Year = {0:04d}'.format(year)) 
            for i in np.arange(0,Ngo.shape[0]):
                jd=astrotime(Ngo['jd'][i],format='jd')
                ymdhms=jd.ymdhms
                if (ymdhms[0] < 1980):
                    year=ymdhms[0]+1900
                    
                    iso = str(year) +'-' + str('{0:02d}'.format(ymdhms[1])) +'-' + str('{0:02d}'.format(ymdhms[2])) \
                                    +' ' + str('{0:02d}'.format(ymdhms[3])) +':' + str('{0:02d}'.format(ymdhms[4])) \
                                    +':' + str('{0:06.3f}'.format(ymdhms[5]))                                          
                    Ngo['jd'][i] = astrotime(iso,format='iso').jd
        if (year == 1981):
            go = np.copy(Ngo)
        else:
            go=np.concatenate((go,Ngo))

    return go

    
def all2():

    for year in np.arange(1981,2021):
        Ngo=np.load('int_go'+str(year)+'.npy')
        if (year <= 1999):
            print('Year = {0:04d}'.format(year))
            njd = list(map(correctJD,Ngo['jd']))
            Ngo['jd'] = np.asarray(njd)
        if (year == 1981):
            go = np.copy(Ngo)
        else:
            go=np.concatenate((go,Ngo))

    return go
        
