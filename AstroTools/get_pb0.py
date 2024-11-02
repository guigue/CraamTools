import numpy as np
import pdb

###################################################
#
# get_pb0
#
# This program computes:
#    he_lon: Longitude at center of disk [deg]
#    he_lat: Latitude at center of disk [deg]
#    Pang:   Position Ang of Rotation Axis [deg]
#    Carr:   Decimal Carrington rotation number
# The only input needed is the julian date (jd)
#
# The code is adapted from original SST (JERCosta) program
# written in IDL (however its origin was a C program)
#
# Author : @guiguesp - on the dark days of SAR-COV-2
#          2020-04-01
#
#################################################

def P(jd):

    t = (jd - 2415020) / 36525
    Carr = (1/27.2753) * (jd-2398167) + 1

    mnl = 279.69668 + 36000.76892 * t + 0.0003025 * t**2
    mnl = mnl % 360

    mna = 358.47583 + 35999.04975 * t - 0.000150 *t**2 - 0.0000033 * t**3
    mna = int(mna) % 360

    e = 0.01675104 - 0.0000418 * t - 0.000000126 * t**2
    c = (1.919460 - 0.004789*t - 0.000014*t**2)*np.sin(np.radians(mna))\
        + (0.020094 - 0.000100*t)*np.sin(np.radians(2*mna)) +  0.000293*np.sin(np.radians(3*mna))
    true_long = (mnl + c) % 360

    ta = (mna + c) % 360

    dist = 1.0000002*(1 - e**2)/(1 + e * np.cos(np.radians(ta)))

    sd = 959.63 / dist

    omega = 259.18 - 1934.142*t
    app_long = true_long - 0.00569 - 0.00479*np.sin(np.radians(omega))

    ob1 = 23.452294 - 0.0130125*t - 0.00000164*t**2 + 0.000000503*t**3

    y = np.cos(np.radians(ob1))*np.sin(np.radians(true_long))
    x = np.cos(np.radians(true_long))
    r = np.sqrt(x**2+y**2)
    true_ra = np.arctan2(y,x) % 360

    if (true_ra < 0):
        true_ra = true_ra + 360
    true_ra = true_ra / 15
    true_dec = np.degrees(np.arcsin(np.sin(np.radians(ob1)) * np.sin(np.radians(true_long))))

    ob2 = ob1 + 0.00256 * np.cos(np.radians(omega))

    y = np.cos(np.radians(ob2)) * np.sin(np.radians(app_long))
    x = np.cos(np.radians(app_long))
    r = np.sqrt(x**2+y**2)
    app_ra = np.arctan2(y,x) % 360
    
    if app_ra < 0 :
        app_ra = app_ra + 360
        
    app_ra = app_ra / 15
    app_dec = np.degrees(np.arcsin(np.sin(np.radians(ob2))*np.sin(np.radians(app_long))))

    theta = (jd - 2398220) * 360/25.38
    i = 7.25
    k = 74.3646 + 1.395833 * t
    lamda = true_long - 0.00569
    lamda2 = lamda - 0.00479 * np.sin(np.radians(omega))
    diff = np.radians((lamda - k))
    x = np.degrees(np.arctan(-np.cos(np.radians(lamda2))*np.tan(np.radians(ob1))))
    y = np.degrees(np.arctan(-np.cos(diff)*np.tan(np.radians(i))))
    Pang = x + y

    he_lat = np.degrees(np.arcsin(np.sin(diff)*np.sin(np.radians(i))))
    y = -np.sin(diff)*np.cos(np.radians(i))
    x = -np.cos(diff)
    r = np.sqrt(x**2+y**2)
    eta = np.degrees(np.arctan2(y,x) )
    he_lon = (eta - theta) % 360
    if (he_lon < 0):
        he_lon = he_lon + 360
    
    return he_lon, he_lat, Pang, Carr
