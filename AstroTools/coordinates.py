import numpy as np
from astropy.time import Time

class sun():

    def __init__(self,sst=True):

        if sst:
            self.observatory = self.sst()

        return
        
    def sst(self):

        return {'latitude':-31.79927778,
                'longitude':-69.30208333,
                'height' : 2.491e+03}
                           
    def hel2xy(self,helio0, b0, sunr,arcmin=False):

        arcsec2rad = np.radians(1/3600)
        lon  = np.radians(helio[0])
        lat  = np.radians(helio[1])
        cb0  = np.cos(b0)
        clat = np.cos(lat)
        clon = np.cos(lon)
        sb0  = np.sin(b0)
        slat = np.sin(lat)
        slon = np.sin(lon)

        crho = clat*cb0*clon + slat*sb0
        crho[crho > 1] = 1.0
        crho[crho < -1] = -1.0
        rho   = np.arccos(crho)
        theta = -np.arctan(cb0*clat*slon,slat-sb0*crho)
        ssunr = np.sin(sunr*arcsec2rad)
        rho1  = np.asin(ssunr*sin(rho))
        for i in np.arange(0,5):
            rho1 = np.arcsin(ssunr*np.sin(rho+rho1))

        y = +rho1*np.cos(theta)/arcsec2rad
        x = -rho1*np.sin(theta)/arcsec2rad

        if arcmin:
            x = x/60.0
            y = y/60.0

        return [x,y]

    def  coordinates(self,jd):

 # Adapted from sunpos.pro
 # Guigue - Feb 2014
 
        t = (jd - 2415020) / 36525
        l = (279.696678 + ((36000.768925 * t) % 360)) *3600
        me = 358.475844 + ((35999.049750 * t) % 360)
        ellcor = (6910.1 - 17.2*t)*np.sin(np.radians(me)) + 72.3 * np.sin(np.radians(2*me))
        l += ellcor
        mv = 212.603219 + ((58517.803875 * t) % 360) 
        vencorr = 4.8 * np.cos(np.radians(299.1017 + mv - me))         + \
                  5.5 * np.cos(np.radians(148.3133 +  2.0 * mv  -  2.0 * me )) + \
          2.5 * np.cos(np.radians(315.9433 +  2.0 * mv  -  3.0 * me )) + \
          1.6 * np.cos(np.radians(345.2533 +  3.0 * mv  -  4.0 * me )) + \
          1.0 * np.cos(np.radians(318.15   +  3.0 * mv  -  5.0 * me ))
          
        l += vencorr
        mm = 319.529425  +  (( 19139.858500 * t)  %  360 )
        marscorr = 2.0 * np.cos(np.radians(343.8883 -  2.0 * mm  +  2.0 * me)) + \
                   1.8 * np.cos(np.radians(200.4017 -  2.0 * mm  + me))
        l += marscorr
        mj = 225.328328  +  (( 3034.6920239 * t)  %  360.0 )
        jupcorr = 7.2 * np.cos(np.radians( 179.5317 - mj + me))               + \
                  2.6 * np.cos(np.radians(263.2167 - mj))                     + \
                  2.7 * np.cos(np.radians( 87.1450 - 2.0 * mj  +  2.0 * me )) + \
                  1.6 * np.cos(np.radians(109.4933 - 2.0 * mj  +  me ))
        l += jupcorr
        d = 350.7376814  + (( 445267.11422 * t)  %  360)
        mooncorr  = 6.5 * np.sin(np.radians(d))
        l += mooncorr
        longterm  = + 6.4 * np.sin(np.radians( 231.19  +  20.20 * t ))
        l  += longterm
        l  =  (l + 2592000) % 1296000
        longmed = l/3600
        l  -=  20.5
        omega = 259.183275 - (( 1934.142008 * t ) % 360.0 )
        l  -=  17.2 * np.sin(np.radians(omega))
        oblt  = 23.452294 - 0.0130125*t + (9.2 * np.cos(np.radians(omega)))/3600
        l = l/3600
        ra  = np.arctan( np.sin(np.radians(l)) * np.cos(np.radians(oblt)) , np.cos(np.radians(l)))
        ra[ra < 0.0] += 2*np.pi
        dec = np.arcsin(np.sin(np.radians(l)) * np.sin(np.radians(oblt)))
        
        return [np.degrees(ra),np.degrees(dec)]

    def azel2eq(self, az, el, hsid):
        
        temp  = np.cos(np.radians(az))*np.cos(np.radians(self.observatory['latitude'])) * np.cos(np.radians(el)) + \
                np.sin(np.radians(self.observatory['latitude']))*np.sin(np.radians(el))
        dec   = np.degrees(np.arcsin(temp))
        temp  = (np.sin(np.radians(el)) - np.sin(np.radians(self.observatory['latitude'])) * np.sin(np.radians(dec)))  / \
                (np.cos(np.radians(self.observatory['latitude'])) * np.cos(np.radians(dec)))
        hour_ang = np.degrees(np.arccos(temp))/15
        if (np.sin(np.radians(az)) > 0):
            hour_ang = -hour_ang 
        ra    = hsid-hour_ang 
        ra    = ((ra + 24) % 24)
        return [ra,dec]


    def eq2azel(self,hsidi,rai,deci):

        # From the SST software
        # Guigue - Feb 2014
  
        ra    =  ((rai/15 + 24) % 24)
        dec   = (deci % 360)
        hsid  = ((hsidi+24) % 24) 
  
        hour_ang = hsid - ra 
        temp     = np.sin(np.radians(self.observatory['latitude'])) * np.sin(np.radians(dec)) + \
                   np.cos(np.radians(self.observatory['latitude'])) * np.cos(np.radians(dec)) * np.cos(np.radians(hour_ang*15.0))
        el       = np.arcsin(temp)
        temp     = (np.sin(np.radians(dec))-np.sin(np.radians(self.observatory['latitude'])) * np.sin(np.radians(el))) / \
                   (np.cos(np.radians(self.observatory['latitude'])) * np.cos(np.radians(el)))
        az       = np.degrees(np.arccos(temp))
        
        if (np.sin(np.radians(hour_ang*15.0)) > 0):
            az=-az
            
        return [az,el]

    def hel2eq(self,lat, lon, P, B, R):

        xlp     = np.radians(P)
        xlp1    = np.sin(xlp) 
        xlp2    = np.cos(xlp) 
        xlb0    = np.radians(B) 
        xlb     = -xlb0 
        xlb1    = np.sin(xlb) 
        xlb2    = np.cos(xlb) 
        xlr     = R 
        hb3     = np.radians(self.observatory['latitude'])
        hb4     = np.sin(hb3) 
        hb5     = np.cos(hb3) 
        hl      = -np.radians(self.observatory['longitude']) 
        hl1     = np.sin(hl) 
        hl2     = np.cos(hl) 
        hy0     = hb4*xlb2+hb5*hl2*xlb1 
        off_ra  = (xlp1*hy0+hb5*hl1*xlp2)*xlr # arcmin
        off_dec = (xlp2*hy0-hb5*hl1*xlp1)*xlr # arcmin
        return [off_ra,off_dec]

    def eq2hel(self, off_ra, off_dec, P,  B,  R):
        
# Adapted from radec_hel.pro in /usr/local/sst/idl
# Guigue, Feb 2014

        xlp  = np.radians(P)
        xlp1 = np.sin(xlp) 
        xlp2 = np.cos(xlp) 
        xlb0 = np.radians(B) 
        xlb  =-xlb0 
        xlb1 = np.sin(xlb) 
        xlb2 = np.cos(xlb) 
        if ((R < 15.72) & (R > 16.7)):
            R=16.4 
        xlr  = R 
        xlv1 = off_dec/xlr #  arcmin/arcmin or other
        xlv2 = off_ra/xlr  #  arcmin/arcmin
        xcom = xlv1*xlv1+xlv2*xlv2 

        if (xcom > 1.0):
            return [-999,-999]
        
        xlv3=np.sqrt(1-xcom) 
        xlxx=((xlv2*xlp1+xlv1*xlp2)*xlb2-xlb1*xlv3) 
        xlv4=np.arctan(xlxx/sqrt(1-xlxx*xlxx)) 
        xlv5=np.cos(xlv4) 
        xlv6=0 
        if (xlv5 != 0):
            xlxx  =  ((xlv2*xlp2-xlv1*xlp1)/xlv5) 
            xlv6  =  np.arctan(xlxx/sqrt(1.0-xlxx*xlxx)) 

        lat  =  np.radians(xlv4)   # latitude South is negative
        lon  = -np.radians(xlv6)   #  longitude East is negative

        #For the heliographic longitude the negative sign is East side ****
        return [lon,lat]


    def LST(self,jd):
    # Adapted from ct2lst in /usr/local/sst/idl/
    # Guigue, Feb 2014
   
        c = [280.46061837, 360.98564736629, 0.000387933, 38710000]  # Meeus constants
        jd2000 = 2451545
        t0 = jd - jd2000
        t = t0/36525
        theta = c[0] + (c[1] * t0) + t^2*(c[2] - t/ c[3] )
        lsidt = ( theta + sst_lon)/15
        lsidt[lsidt < 0] = 24 + (lsidt[lsidt < 0] % 24)
        lsidt = lsidt % 24
        return lsidt

    def get_pb0(self,get_pb0, jd):
#
# NAME:
#       GET_pb0
# PURPOSE:
#	Provides the inclination of solar x (p) and sun's center latitude (b0).
#       ; CATEGORY:
#
# CALLING SEQUENCE:
#       OUT = GET_PB0(JD)
# INPUTS:
#       JD -	Julian Day
# KEYWORD INPUTS:
#     
# OUTPUTS:
#	DATA = Vector of solar ephemeris data:
#	  DATA(0) = Distance (AU).
#	  DATA(1) = Semidiameter of disk (sec).
#	  DATA(2) = Longitude at center of disk (deg).
#	  DATA(3) = Latitude at center of disk (deg).
#	  DATA(4) = Position angle of rotation axis (deg).
#	  DATA(5) = decimal carrington rotation number.
# KEYWORD OUTPUTS
# COMMON BLOCKS:
# NOTES:
#       Notes: based on the book Astronomical Formulae
#              for Calculators, by Jean Meeus.
#              Now copied from Yohkoh routine get_sun.pro to here
#	       and modified to give only a short list of params
# MODIFICATION HISTORY:
#	Joaquim E.R Costa 13-Nov-97
# Copyright (C) 1991, Johns Hopkins University/Applied Physics Laboratory
# This software may be used, copied, or redistributed as long as it is not
# sold and this copyright notice is reproduced on each copy made.  This
# routine is provided as is without any express or implied warranties
# whatsoever.  Other limitations apply as described in the file disclaimer.txt.
#-

#Julian Centuries from 1900.0:
        t = (jd - 2415020 )/36525

# Carrington Rotation Number:
        carr = (1/27.2753)*(jd-2398167) + 1

# Geometric Mean Longitude (deg):

        mnl = 279.69668 + 36000.76892*t + 0.0003025*t**2
        mnl = (mnl % 360)

# Mean anomaly (deg):
        mna = 358.47583 + 35999.04975*t - 0.000150*t**2 - 0.0000033*t**3
        mna = mna % 360

# Eccentricity of orbit:
        e = 0.01675104 - 0.0000418*t - 0.000000126*t**2

# Sun's equation of center (deg):
        c = (1.919460 - 0.004789 * t - 0.000014*t**2)*np.sin(np.radians(mna)) + \
            (0.020094 - 0.000100 * t)* np.sin(np.radians(2*mna)) + \
             0.000293 * np.sin(np.radians(3*mna))

# Sun's true geometric longitude (deg)
#   (Refered to the mean equinox of date.  Question: Should the higher
#    accuracy terms from which app_long is derived be added to true_long?)
        true_long = (mnl + c) % 360

# Sun's true anomaly (deg):
        ta = (mna + c) % 360

# Sun's radius vector (AU).  There are a set of higher accuracy
#   terms not included here.  The values calculated here agree with
#   the example in the book:
        dist = 1.0000002*(1 - e**2)/(1 + e*np.cos(np.radians(ta)))

# Semidiameter (arc sec):
        sd = 959.63/dist

# Apparent longitude (deg) from true longitude:
        omega = 259.18 - 1934.142*t		# Degrees
        app_long = true_long - 0.00569 - 0.00479*np.sin(np.radians(omega))

# Latitudes (deg) for completeness.  Never more than 1.2 arc sec from 0,
#   always set to 0 here:
        true_lat = fltarr(n_elements(dist))
        app_lat = fltarr(n_elements(dist))

# True Obliquity of the ecliptic (deg):
        ob1 = 23.452294 - 0.0130125*t - 0.00000164*t**2 + 0.000000503*t**3

# True RA, Dec (is this correct?):
        y = np.cos(np.radians(ob1))*np.sin(np.radians(true_long))
        x = np.cos(np.radians(true_long))
        
        r,true_ra = self.recpol(x, y)
        
        true_ra = true_ra % 360
        true_ra[true_ra < 0] +=  360
        true_ra = true_ra/15
        true_dec = np.degrees(np.arcsin(np.sin(np.radians(ob1)) * np.sin(np.radians(true_long))))

# Apparent  Obliquity of the ecliptic:
        ob2 = ob1 + 0.00256 * np.cos(np.radians(omega))	# Correction.

# Apparent  RA, Dec (agrees with example in book):
        y = np.cos(np.radians(ob2))*np.sin(np.radians(app_long))
        x = np.cos(np.radians(app_long))
        
        r, app_ra = self.recpol(x, y)
        app_ra = app_ra % 360
  
        app_ra[app_ra < 0] += 360
        app_ra = app_ra/15
        app_dec = np.degrees(np.arcsin(np.sin(np.radians(ob2))*np.sin(np.radians(app_long))))

# Heliographic coordinates:
        theta = (jd - 2398220)*360/25.38	# Deg.
        i = 7.25				# Deg.
        k = 74.3646 + 1.395833*t		# Deg.
        lamda = true_long - 0.00569
        lamda2 = lamda - 0.00479*np.sin(np.radians(omega))
        diff = np.degrees(lamda - k)
        x = np.degrees(np.arctan(-np.cos(np.radians(lamda2))*np.tan(np.radians(ob1))))
        y = np.degrees(np.arctan(-np.cos(diff)*np.tan(np.radians(i))))

# Position of north pole (deg):
        pa = x + y

# Latitude at center of disk (deg):
        he_lat = np.degrees(np.arcsin(np.sin(diff)*np.sin(np.radians(i))))

# Longitude at center of disk (deg):
        y = -np.sin(diff)*np.cos(np.radians(i))
        x = -np.cos(diff)
        r,eta = self.rec_pol(x, y)
        he_lon = (eta - theta) % 360
        he_lon[he_lon < 0] += 360
        return [dist,sd,he_lon,he_lat,pa,carr]

    def recpol(self,x,y):

         return np.sqrt(x**2 + y**2), np.degrees(np.arctan2(y, x))
     
    def offeq2xy(self, eq_off, obstime):

#
# NAME     SST2XY
#
# PURPORSE Convert SST horizontal (AZ,EL) coordinates to offsets (X,Y) respect to 
#          the center of Sun disk.  Input are the equatorial offsets of the tracked
#          active region at 0 UT, and the julian day for the outputs.  If a Beam Pos
#          structure is given offsets for the beams are calculated.  Also a solution 
#          of the multibeam system can be given to get offsets of the gyrosynchrotron source.
#
# INPUTS   eq_off:  [right ascention, declination] offsets @ 0 UT of the tracked AR. IN ARC MINUTES
#          jd:      julian day for the outputs
#
# OUTPUT KEYWORDS 
#          ar_hel:  [lon,lat] in Heliographics system of the tracked AR at JD time. DEGREES
#          ar_xy:   [X,Y] offsets of the tracked AR at JD time. ARC SECONDS
#          ar:      Tracked AR
#                   [HLON0UT, HLAT0UT] Heliographics Longitude and Latitud at 0 UT. In DEGREES.
#                   [EQ_RA0UT,EQ_DEC0UT] Equatorial offsets from Sun Center at 0 UT. In ARC MINUTES
#                   [AZ,EL]   horizontal in DEGREES at JD
#                   [RA,DEC]  equatorial offsets from Sun center at JD. In ARC MINUTES
#                   [X,Y]     East-West, North-South offsets from disk center at JD. In ARC SECONDS
#                   [LON,LAT] Heliographics Longitude and Latitud at JD. In DEGREES.
#
# HISTORY
#        First written by Guigue using many different routines, of SST, SSW and JPL Astronomical IDL
#        Library.  See comments.
#        February, 28, 2014
#

        year    = obstime.datetime.year
        month   = obstime.datetime.month
        day     = obstime.datetime.day
        hour    = obstime.datetime.hour
        minute  = obstime.datetime.minute
        second  = obstime.datetime.second
        mstime  = hour * 36000000 + minute * 600000 + second * 10000
        s_eq_coords= self.coordinates(obstime.jd)                      # Sun ephemeris

        jd0hs      = Time(year,month,day,0)                            # compute julday at 0 UT
        pbr        = self.get_pb0(jd0hs)                               # get Sun Radius, P & B angles
        sunr       = pbr[1]
        pang       = pbr[4]
        bang       = pbr[3]
        ar_coord   = self.eq2hel(eq_off[0],eq_off[1],pang,bang,sunr/60.0)  # get heliographics coordinates of AR
        psi        = 14.40 - 1.8*np.sin(np.radians(ar_coord[1]))**2 - 2.4*np.sin(np.radians(ar_coord[1]))**2 # angular speed at AR latitude 
        dlong      = psi * (hour+ minute/60.0+second/3600.0)/24.0                        # rotation angle since 0 UT
        rar_hcoord = ar_coord+[dlong,0.0]                                                # rotated AR Helio coordinates 
        rar_eqoff  = hel2eq(rar_hcoord[1], rar_hcoord[0],pang,bang,sunr/60.0)
        xy_ar      = hel2xy(rar_hcoord,bang*rad,sunr)                                    # XY offsets

        rar_ecoord = [rar_eqoff[0]/60 + s_eq_coords[0] , rar_eqoff[1]/60 + s_eq_coords[1]]
        rar_azel   = self.eq2azel(self.LST(jd),rar_ecoord[0],rar_ecoord[1])
        rar_azel[0]= ((rar_azel[0]+360) % 360)
        sun_azel   = self.eq2azel(self.LST(jd),s_eq_coords[0],s_eq_coords[1])
        sun_azel[0]= ((sun_azel[0]+360) % 360)

        return


