FUNCTION HEL2XY, helio0, b0, sunr

arcsec2rad = 4.84813681d-6

siz = size(helio0)
typ = siz( siz(0)+1 )
if (typ eq 7) then helio = conv_hs2h(helio0) else helio = helio0
we = reform(helio[0,*])
ns = reform(helio[1,*])
nout = n_elements(ns)

lon = we*!dpi/180.0d0
lat = ns*!dpi/180.0d0

we  = 0				;reduce memory requirements
ns  = 0

cb0  = cos(b0)
clat = cos(lat)
clon = cos(lon)
sb0  = sin(b0)
slat = sin(lat)
slon = sin(lon)

crho = clat*cb0*clon + slat*sb0
badplus = where(crho GT 1.0d0,nbadplus)
badminus = where(crho LT -1.0d0,nbadminus)
if nbadplus GT 0 then crho[badplus] = 1.0d0
if nbadminus GT 0 then crho[badminus] = -1.0d0
rho = acos(crho)

theta = -atan(cb0*clat*slon,slat-sb0*crho)

ssunr = sin(sunr*arcsec2rad)
rho1 = asin(ssunr*sin(rho))
for i=0,5 do rho1 = asin(ssunr*sin(rho+rho1))

y = +rho1*cos(theta)/arcsec2rad
x = -rho1*sin(theta)/arcsec2rad

if keyword_set(arcmin) then begin
   x = x/60.0
   y = y/60.0
endif

blss = where(rho GT ((!dpi/2.0d0)-rho1),nblss)
behind = byte(rho) & behind[*] = 0
if nblss GT 0 then behind[blss] = 1b

xy = double(helio)
xy[0,*] = x
xy[1,*] = y

return, xy

END


FUNCTION SUNCOORDS,jd, longmed, oblt

 ; Adapted from sunpos.pro
 ; Guigue - Feb 2014
 
 dtor = !DPI/180.0d       ;(degrees to radian, double precision)
 t = (jd - 2415020.0d)/36525.0d0
 l = (279.696678d0+((36000.768925d0*t) mod 360.0d0))*3600.0d0
 me = 358.475844d0 + ((35999.049750D0*t) mod 360.0d0)
 ellcor  = (6910.1d0 - 17.2D0*t)*sin(me*dtor) + 72.3D0*sin(2.0D0*me*dtor)
 l = l + ellcor
 mv = 212.603219d0 + ((58517.803875d0*t) mod 360.0d0) 
 vencorr = 4.8D0 * cos((299.1017d0 + mv - me)*dtor) + $
          5.5D0 * cos((148.3133d0 +  2.0D0 * mv  -  2.0D0 * me )*dtor) + $
          2.5D0 * cos((315.9433d0 +  2.0D0 * mv  -  3.0D0 * me )*dtor) + $
          1.6D0 * cos((345.2533d0 +  3.0D0 * mv  -  4.0D0 * me )*dtor) + $
          1.0D0 * cos((318.15d0   +  3.0D0 * mv  -  5.0D0 * me )*dtor)
 l = l + vencorr
 mm = 319.529425d0  +  (( 19139.858500d0 * t)  mod  360.0d0 )
 marscorr = 2.0d0 * cos((343.8883d0 -  2.0d0 * mm  +  2.0d0 * me)*dtor ) + $
            1.8D0 * cos((200.4017d0 -  2.0d0 * mm  + me) * dtor)
 l = l + marscorr
 mj = 225.328328d0  +  (( 3034.6920239d0 * t)  mod  360.0d0 )
 jupcorr = 7.2d0 * cos(( 179.5317d0 - mj + me )*dtor) + $
          2.6d0 * cos((263.2167d0  -  MJ ) *dtor) + $
          2.7d0 * cos(( 87.1450d0  -  2.0d0 * mj  +  2.0D0 * me ) *dtor) + $
          1.6d0 * cos((109.4933d0  -  2.0d0 * mj  +  me ) *dtor)
 l = l + jupcorr
 d = 350.7376814d0  + (( 445267.11422d0 * t)  mod  360.0d0 )
 mooncorr  = 6.5d0 * sin(d*dtor)
 l = l + mooncorr
 longterm  = + 6.4d0 * sin(( 231.19d0  +  20.20d0 * t )*dtor)
 l  =    l + longterm
 l  =  ( l + 2592000.0d0)  mod  1296000.0d0 
 longmed = l/3600.0d0
 l  =  l - 20.5d0
 omega = 259.183275d0 - (( 1934.142008d0 * t ) mod 360.0d0 )
 l  =  l - 17.2d0 * sin(omega*dtor)
 oblt  = 23.452294d0 - 0.0130125d0*t + (9.2d0*cos(omega*dtor))/3600.0d0
 l = l/3600.0d0
 ra  = atan( sin(l*dtor) * cos(oblt*dtor) , cos(l*dtor) )
 neg = where(ra LT 0.0d0, Nneg) 
 if Nneg GT 0 then ra[neg] = ra[neg] + 2.0d*!DPI
 dec = asin(sin(l*dtor) * sin(oblt*dtor))
 return, [ra/dtor,dec/dtor]
END


FUNCTION AZEL2EQ, az, el, hsid

  rad   = !DPI/180.0d       
  long  = -69.30208333D
  lat   = -31.79927778D  
  temp  = cos(az*rad)*cos(lat*rad)*cos(el*rad)+sin(lat*rad)*sin(el*rad)
  dec   =  asin(temp)/rad 
  temp  = (sin(el*rad)-sin(lat*rad)*sin(dec*rad))/(cos(lat*rad)*cos(dec*RAD))
  hour_ang = acos(temp)/rad/15.0D
  if (sin(az*rad) ge 0.0) then hour_ang = -hour_ang 
  ra    = hsid-hour_ang 
  ra    = ((ra + 24.0D) mod 24.0D )
  return, [ra,dec]

END



FUNCTION EQ2AZEL,hsidi,rai,deci

  ; From the SST software
  ; Guigue - Feb 2014
  
  rad   = !DPI/180.0d       
  long  = -69.30208333D
  lat   = -31.79927778D  
  ra    =  ((rai/15.0D + 24.0D) mod 24.0D)
  dec   = (deci mod 360.0D)
  hsid  = ((hsidi+24.0D) mod 24.0D) 
  
  hour_ang = hsid - ra 
  temp     = sin(lat*rad)*sin(dec*rad)+cos(lat*rad)*cos(dec*rad)*cos(hour_ang*15.0*rad)
  el       = asin(temp)/rad 
  temp     = (sin(dec*rad)-sin(lat*rad)*sin(el*rad))/(cos(lat*rad)*cos(el*rad))
  az       = acos(temp)/rad 
  for j=0l,n_elements(hour_ang)-1 do if (sin(hour_ang[j]*15.0*rad) ge 0.0) then  az[j]=-az[j] 
  return, [az,el]

END

FUNCTION HEL2EQ,lat, lon, P, B, R

  rad     = !dtor 
  xlp     = P*rad 
  xlp1    = sin(xlp) 
  xlp2    = cos(xlp) 
  xlb0    = B*rad 
  xlb     = -xlb0 
  xlb1    = sin(xlb) 
  xlb2    = cos(xlb) 
  xlr     = R 
  hb3     = lat*rad 
  hb4     = sin(hb3) 
  hb5     = cos(hb3) 
  hl      = -lon*rad 
  hl1     = sin(hl) 
  hl2     = cos(hl) 
  hy0     = hb4*xlb2+hb5*hl2*xlb1 
  off_ra  = (xlp1*hy0+hb5*hl1*xlp2)*xlr ;  arcmin
  off_dec = (xlp2*hy0-hb5*hl1*xlp1)*xlr ; arcmin
  return, [off_ra,off_dec]

END

FUNCTION EQ2HEL, off_ra, off_dec, P,  B,  R
; Adapted from radec_hel.pro in /usr/local/sst/idl
; Guigue, Feb 2014

  rad=!dtor

  xlp  = P*rad 
  xlp1 = sin(xlp) 
  xlp2 = cos(xlp) 
  xlb0 = B*rad 
  xlb  =-xlb0 
  xlb1 = sin(xlb) 
  xlb2 = cos(xlb) 
  if ((R lt 15.72) and (R gt 16.7)) then R=16.4 
  xlr  = R 
  xlv1 = off_dec/xlr ;  arcmin/arcmin or other
  xlv2 = off_ra/xlr  ;  arcmin/arcmin
  xcom = xlv1*xlv1+xlv2*xlv2 

  if (xcom gt 1.0) then return, [-999,-999]

  xlv3=sqrt(1.00000000001-xcom) 
  xlxx=((xlv2*xlp1+xlv1*xlp2)*xlb2-xlb1*xlv3) 
  xlv4=atan(xlxx/sqrt(1.0-xlxx*xlxx)) 
  xlv5=cos(xlv4) 
  xlv6=0.0 
  if (xlv5 ne 0) then begin
    xlxx  =  ((xlv2*xlp2-xlv1*xlp1)/xlv5) 
    xlv6  =  atan(xlxx/sqrt(1.0-xlxx*xlxx)) 
  endif
  lat  =  xlv4/rad   ;   latitude South is negative
  lon  = -xlv6/rad  ;   longitude East is negative

  ;For the heliographic longitude the negative sign is East side ****
   return, [lon,lat]
END

FUNCTION LST,jd
   ; Adapted from ct2lst in /usr/local/sst/idl/
   ; Guigue, Feb 2014
   
   sst_lon = -69.29669444
   c = [280.46061837d0, 360.98564736629d0, 0.000387933d0, 38710000.0 ]  ; Meeus constants
   jd2000 = 2451545.0D0
   t0 = jd - jd2000
   t = t0/36525
   theta = c[0] + (c[1] * t0) + t^2*(c[2] - t/ c[3] )
   lsidt = ( theta + double(sst_lon))/15.0d
   neg = where(lsidt lt 0.0D0, n)
   if n gt 0 then lsidt[neg] = 24.D0 + (lsidt[neg] mod 24)
   lsidt = lsidt mod 24.D0   
   return, lsidt
END

PRO SST2XY, eq_off, jd, sour = s, bpos = b, beams = beams, ar = ar, gss = gss, pntb=pntb, quiet=quiet

;+
;
; NAME     SST2XY
;
; PURPORSE Convert SST horizontal (AZ,EL) coordinates to offsets (X,Y) respect to 
;          the center of Sun disk.  Input are the equatorial offsets of the tracked
;          active region at 0 UT, and the julian day for the outputs.  If a Beam Pos
;          structure is given offsets for the beams are calculated.  Also a solution 
;          of the multibeam system can be given to get offsets of the gyrosynchrotron source.
;
; INPUTS   eq_off:  [right ascention, declination] offsets @ 0 UT of the tracked AR. IN ARC MINUTES
;          jd:      julian day for the outputs
;
; INPUTS KEYWORDS
;          sour:    structure with multibeam solutions, as produced by sstpos.pro. It should have at least
;                   [OFF,EL]  horizontal coordinates of the antenna in DEGREES
;                   TIME      time since 0UT, in 0.1 milliseconds
;          bpos:    structure with beam positions. When given, computes the beam positions in XY coordinates.
;          sour:    structure with the multibeam solution.  When given computes the source positions in XY coordinates.
;          pntb:    Pointing Beam, from 1 to 6. Default is 5
;
;
; OUTPUT KEYWORDS 
;          ar_hel:  [lon,lat] in Heliographics system of the tracked AR at JD time. DEGREES
;          ar_xy:   [X,Y] offsets of the tracked AR at JD time. ARC SECONDS
;          beams:   structure with coordinates of the beams at JD in many different frames
;                   [AZ,EL]  horizontal in DEGREES
;                   [RA,DEC] equatorial offsets from Sun center. In ARC MINUTES
;                   [X,Y]    East-West, North-South offsets from disk center. In ARC SECONDS
;          gss:     structure with the coordinates of the multibeam solution in different frames.
;                   [AZ,EL]   horizontal in DEGREES
;                   [RA,DEC]  equatorial offsets from Sun center. In ARC MINUTES
;                   [X,Y]     East-West, North-South offsets from disk center. In ARC SECONDS
;                   [LON,LAT] Heliographics Longitude and Latitud. In DEGREES.
;          ar:      Tracked AR
;                   [HLON0UT, HLAT0UT] Heliographics Longitude and Latitud at 0 UT. In DEGREES.
;                   [EQ_RA0UT,EQ_DEC0UT] Equatorial offsets from Sun Center at 0 UT. In ARC MINUTES
;                   [AZ,EL]   horizontal in DEGREES at JD
;                   [RA,DEC]  equatorial offsets from Sun center at JD. In ARC MINUTES
;                   [X,Y]     East-West, North-South offsets from disk center at JD. In ARC SECONDS
;                   [LON,LAT] Heliographics Longitude and Latitud at JD. In DEGREES.
;
; HISTORY
;        First written by Guigue using many different routines, of SST, SSW and JPL Astronomical IDL
;        Library.  See comments.
;        February, 28, 2014
;
;+

rad = !dpi/1.80D+02

if keyword_set(pntb) then begin
   if (pntb lt 1) or (pntb gt 6) then pntch=4 else pntch=pntb-1
endif else pntch=4

if not keyword_set(quiet) then print,'Pointing Beam = '+string(pntch+1,format='(i3)')
  
caldat,jd,month,day,year,hour,minute,second                                      ; extract date information from jd
ssttag  = string(year[0]-1900,format='(i3.3)')+string(month[0],format='(i2.2)')+string(day[0],format='(i2.2)')
mstime  = long(hour)*36000000l+long(minute)*600000l+long(second)*10000l

s_eq_coords= SunCoords(jd)                                                       ; Sun ephemeris

jd0hs      = julday(month,day,year,0d0)                                          ; compute julday at 0 UT
pbr        = get_pb0(jd0hs)                                                      ; get Sun Radius, P & B angles
sunr       = pbr[1]
pang       = pbr[4]
bang       = pbr[3]
ar_coord   = eq2hel(eq_off[0],eq_off[1],pang,bang,sunr/60.0)                     ; get heliographics coordinates of AR
psi        = 14.40 - 1.8*sin(ar_coord[1]*rad)^2 - 2.4*sin(ar_coord[1]*rad)^2 ; angular speed at AR latitude 
dlong      = psi * (hour+ minute/60.0+second/3600.0)/24.0                        ; rotation angle since 0 UT
rar_hcoord = ar_coord+[dlong,0.0]                                                ; rotated AR Helio coordinates 
rar_eqoff  = hel2eq(rar_hcoord[1], rar_hcoord[0],pang,bang,sunr/60.0)
xy_ar      = hel2xy(rar_hcoord,bang*rad,sunr)                                    ; XY offsets

rar_ecoord = [rar_eqoff[0]/60.0D + s_eq_coords[0] , rar_eqoff[1]/60.0D + s_eq_coords[1]]
rar_azel   = eq2azel(lst(jd),rar_ecoord[0],rar_ecoord[1])
rar_azel[0]= ((rar_azel[0]+360.0) mod 360.0)
sun_azel   = eq2azel(lst(jd),s_eq_coords[0],s_eq_coords[1])
sun_azel[0]= ((sun_azel[0]+360.0) mod 360.0)

strname    = 'AR'+ssttag
cmd        = 'ar={'+strname+',time:mstime,jd:jd,hlon0ut:ar_coord[0], hlat0ut:ar_coord[1],'     + $
             'eq_ra0ut:eq_off[0],eq_dec0ut:eq_off[1], hlon:rar_hcoord[0],hlat:rar_hcoord[1],' + $
             'ra:rar_eqoff[0],dec:rar_eqoff[1],az:rar_azel[0],el:rar_azel[1],x:xy_ar[0],y:xy_ar[1]}'
err        = execute(cmd)      
              
if keyword_set(b) then begin
  b_ = b
  b_.el = - b_.el

  pntb=4
  
  strname    = 'BEAM'+ssttag
  cmd        = 'beams={'+strname+', az:dblarr(6), el:dblarr(6), ra:dblarr(6), dec:dblarr(6), hlat:dblarr(6),'+ $
                'hlon:dblarr(6), x:dblarr(6), y:dblarr(6)}'
  err        = execute(cmd)
  
  for i=0,5 do begin
    beams.az[i]   = rar_azel[0] - (b_.off[i]-b_.off[pntb])/(60.0 * cos(rar_azel[1]*rad)) 
    beams.el[i]   = rar_azel[1] - (b_.el[i]-b_.el[pntb])/60.0 
    rr            = azel2eq(beams.az[i],beams.el[i],lst(jd))
    beams.dec[i]  = (rr[1]-s_eq_coords[1])*60.0
    beams.ra[i]   = (rr[0]*15.0-s_eq_coords[0]) * cos(rar_ecoord[1]*rad) * 60.0 
    beams.x[i]    =  (-beams.ra[i]*cos(pang*rad) + beams.dec[i]*sin(pang*rad))*60.0
    beams.y[i]    = (beams.ra[i]*sin(pang*rad) + beams.dec[i]*cos(pang*rad))*60.0
  endfor

endif

if (keyword_set(s) and keyword_set(b)) then begin
  s_      = s
  s_.el   = -s_.el
  
  b_      = b
  b_.el   = - b_.el

  pntb    = 4  
  n       = n_elements(s)
  
  strname = 'MB'+ssttag
  cmd     = 'gss={'+strname+',time:0l,az:0.0D+00, el:0.0D+00, ra:0.0D+00, dec:0.0D+00,' + $
             'hlat:0.0D+00, hlon:0.0D+00, x:0.0D+00,y:0.0D+00}'
  err     = execute(cmd)           
  gss     = replicate(gss,n)
  
  for i=0l,n-1 do begin
    jdi         = jd0hs + double(s_[i].time)/8.64D+08
    s_eq_coords = SunCoords(jdi)
    rar_azel   = eq2azel(lst(jdi),rar_ecoord[0],rar_ecoord[1])

    gss[i].time = s_[i].time
    gss[i].az   = rar_azel[0] - (s_[i].off-b_.off[pntb])/(60.0*cos(rar_azel[1]*rad))
    gss[i].el   = rar_azel[1] - (s_[i].el-b_.el[pntb])/60.0
    rr          = azel2eq(gss[i].az,gss[i].el,lst(jdi))
    gss[i].ra   = (rr[0]*1.5D+01-s_eq_coords[0])*cos(s_eq_coords[1]*rad)*60.0
    gss[i].dec  = (rr[1]-s_eq_coords[1])*60.0
    ss          = eq2hel(gss[i].ra,gss[i].dec,pang,bang,sunr/60.0)
    gss[i].hlat = ss[1]
    gss[i].hlon = ss[0]
    xy          = hel2xy([gss[i].hlon,gss[i].hlat],bang*rad,sunr)
    gss[i].x    = xy[0]
    gss[i].y    = xy[1]
  endfor
  
endif

return
end
