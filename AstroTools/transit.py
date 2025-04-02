import astropy.time as time
from astropy.time import Time
import datetime
from astropy import units as u
from astropy.coordinates import solar_system_ephemeris, EarthLocation, get_body, AltAz
from matplotlib import pyplot as plt
from matplotlib import dates
import numpy as np

import CASLEO

__version__ = '2025-04-02T1000BRT'

def version():

  return __version__

  
def compute(object='sun',when=Time.now(),plotfig=False,savefig=False):
  
  year=when.datetime.year
  month=when.datetime.month
  day=when.datetime.day
  az=[]
  el=[]
  t=[]
  dt = time.TimeDelta(15*u.minute)
  time0 = Time(str(year)+'-'+str(month)+'-'+str(day)+' '+'09:00:00')
  tt = time0 
  while (tt.datetime.day == when.datetime.day): 
    t.append(tt.datetime)
    with solar_system_ephemeris.set('jpl'):
      sun_radec=get_body(object,tt,CASLEO.Observatory_Coordinates())
    aa=AltAz(location=CASLEO.Observatory_Coordinates(),obstime=tt)
    altaz = sun_radec.transform_to(aa)
    az.append(altaz.az.value)
    el.append(altaz.alt.value)
    tt = tt + dt

  az = np.asarray(az)
  el = np.asarray(el)
  t  = np.asarray(t)

  print('\n\n Maximum Elevation = {0:5.2f} deg at {1:02d}:{2:0d} UT\n\n'.format(el.max(),
                                                                         t[el.argmax()].hour,
                                                                         t[el.argmax()].minute))
                                                                         
  if plotfig:
    up = (el > 0)
    fig,ax=plt.subplots(figsize=(10,7))
    ax.plot(t[up],el[up],'-k')
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%m'))
    ax.set_title('Apparent Elevation for '+ object + ' on ' + time0.strftime("%Y-%m-%d") + ' at CASLEO')
    ax.plot([t[el.argmax()],t[el.argmax()]],[0,el.max()+5],'--r')
    ax.text(t[el.argmax()]+datetime.timedelta(minutes=15),0,'Meridian Transit at: '+ t[el.argmax()].strftime("%H:%M")+' UT')

    if savefig:
      fname = 'Transit-'+object+'-'+time0.strftime("%Y-%m-%d")+'.pdf'
      plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
                  orientation='portrait', papertype=None, format='pdf',
                  transparent=False, bbox_inches=None, pad_inches=0.1,
                  frameon=None, metadata=None)


    
  return t,az,el
    


