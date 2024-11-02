import astropy.time as time
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import solar_system_ephemeris, EarthLocation, get_body, AltAz
from matplotlib import pyplot as plt
from matplotlib import dates
import numpy as np

import CASLEO


def sun_transit(when=Time.now()):
  year=when.datetime.year
  month=when.datetime.month
  day=when.datetime.day
  az=[]
  el=[]
  t=[]
  dt = time.TimeDelta(15*u.minute)
  time0 = Time(str(year)+'-'+str(month)+'-'+str(day)+' '+'09:00:00')
  for i in np.arange(56):
    tt = time0 + i * dt
    t.append(tt.datetime)
    with solar_system_ephemeris.set('jpl'):
      sun_radec=get_body('sun',tt,CASLEO.Observatory_Coordinates())
    aa=AltAz(location=CASLEO.Observatory_Coordinates(),obstime=tt)
    altaz = sun_radec.transform_to(aa)
    az.append(altaz.az.value)
    el.append(altaz.alt.value)

  fig,ax=plt.subplots(figsize=(10,7))
  ax.plot(t,el)
  ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%m'))
    
  return t,az,el
    


