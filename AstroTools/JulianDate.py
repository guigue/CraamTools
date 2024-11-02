# External methods
import sys, string, os, struct, glob
import numpy as np
import datetime

def fromIsotime(isotime):
    """
        fromIsotime

        Aim: to return the Julian Date after 1999,
         i.e. it does not computes dates before the
         Gregorian Reformation

        Input: isotime is a string with a time in ISO format
           ISO time format = YYYY-MM-DD hh:mm:ss.ss

        Output: the julian date

        Examples:
        >>> import astronomical_methods as am
        >>> print am.julian_date('2014-10-11 13:05:22.5645678')
        >>>2456942.0454
 
        Change record:
        - Adapted from SST procedure in ephem_lib.c (AM & JERC)
        - First written by Guigue on 19 February 2015, in Sao Paulo
        - Changes on the extraction of hours from the input string - 2017-06-15 (Guigue)

    """
    # Initial Values
    year = 0.0
    month = 0.0
    day = 0.0
    hh = 0.0
    mm = 0.0
    ss = 0.0

    try: 
        ymd,time = isotime.strip().split(" ",1)
    except:
        ymd = isotime.strip()
        time = ''

    year  = float(ymd[0:4])
    month = float(ymd[5:7])
    day   = float(ymd[8:10])

    ltime = len(time)
    if ltime :
        if ltime > 5 : ss = float(time[6:]) 
        if ltime > 2 : mm = float(time[3:5])
        hh = float(time[0:2])  
        f_of_day = (hh+mm/60.0+ss/3600.0)/24.0   
    else:
        f_of_day = 0.0
        
    # Gregorian adopted in Oct., 15, 1582
    greg = 15 + 31 * (10 + 12 * 1582) ;
            
    if(month > 2):
        jy = year 
        jm = month + 1 
    else:
        jy = year - 1 
        jm = month + 13 

    jd = int(365.25 * jy) + int(30.6001 * jm) + day + 1720995

    if ( day + (31 * (month + 12 * year )) >= greg ) :
        ja = int(0.01 * jy)
        jd = jd + 2 -ja + int(0.25 * ja) 
                    

    return float(jd)+f_of_day-0.5 
