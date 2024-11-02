# External methods
import sys, string, os, struct, glob
import numpy as np
import datetime
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from CraamTools.filters import RunningMean
import pdb

""" ------------------------------------------------------------------

CIRC_FIT

   A series of procedures to fit a circle from a data set.
   It derives from the old and good SST's 'circle_fit.pro'
   The sequence of use is as follows. Assuming you have an image
   M[x,y]

   limb=circ_fit.Get_Limb(M)
   par,cov = circ_fit.fit(limb)

   par is a 3 element list with the results
   Circle Center: Col = par[0], Row = par[1]
   Circle Radius: par[2]

   cov is covariance matrix.

Author:  @Guiguesp
         2018-02-07 after too many mistakes >-(

------------------------------------------------"""

def circle_func(phi,x0,y0,r):
    """
    CIRCLE

       The fitting function used by curve_fit

    """

    b   = x0 * np.cos(phi) + y0 * np.sin(phi)
    wur = np.sqrt( b*b - x0*x0 - y0*y0 + r*r)
    d   = b + wur

    return d

def fit(x,y,verbose=False):

    ierase   = np.array([])
    Niter    = 0
    Nmin     = 50
    MaxNiter = 5
    Continue = True

    while Continue:

        # Prepare data for polar fitting
        # First guess of circle center
        if (ierase.shape[0]) > 1 :
            x = np.delete(x,ierase)
            y = np.delete(y,ierase)

        Col_0     = x.mean()
        Row_0     = y.mean()

        # Distance of limb to disc center
        dx = x-Col_0
        dy = y-Row_0
        d  = np.sqrt( dx**2 + dy**2)

        # First guess for the circle radius
        R_0    = d.mean()

        phi = np.arctan2(dy,dx)
        so  = np.argsort(phi)
        phi = phi[so]
        x   = x[so]
        y   = y[so]
        d   = d[so]

        # First guess of parameters. Note that, since we subtracted the circle center
        # to the limb, the center first guess is (0,0)
        par0     = (0,0,R_0)

        par, cov = curve_fit(circle_func, phi,d,p0=par0)
        dfit     = circle_func(phi,par[0],par[1],par[2])
        fsigma   = np.std(d-dfit)
        psigma   = np.sqrt(np.diag(cov))
        ierase   = np.ravel(np.where(np.abs(d-dfit) > fsigma))
        Niter    += 1

        if verbose:
            print(' ')
            print(' Iteration # {0:3d}: fit sigma = {1:5f}, R_0 = {2:5f}+-{3:5f}'.format(Niter,fsigma,par[2],psigma[2]))

        if (len(x) == (x.shape[0])) | (Niter > MaxNiter) | (x.shape[0] < 20)  :
            Continue = False

    # Return the disc center in the matrix coordinates adding what we subtracted
    par[0]   += Col_0
    par[1]   += Row_0
    return par,cov

def get_limb(mx,my,m,smooth=0):
    nbins = 200
    h,xh  = np.histogram(np.ravel(m),bins=nbins)
    xh = xh[1:]
    if smooth > 0:
        h = RunningMean.rm1d(h,smooth)
        xh = xh[xh.shape[0]-h.shape[0] : ]
    nn    = nbins//2        
    i     = np.argmax(h[:nn])
    sk    = xh[i]
    i     = np.argmax(h[nn+1:])
    qs    = xh[nn+1+i]
    ll    = 0.5*(sk+qs)
    fig,ax = plt.subplots()
    c   = ax.contour(mx,my,m,[ll])
    plt.close()

    x   = np.ndarray(0)
    y   = np.ndarray(0)

    for seg in np.arange(len(c.allsegs[0])):
        x = np.concatenate((x,c.allsegs[0][seg][:,0]))
        y = np.concatenate((y,c.allsegs[0][seg][:,1]))

    return x,y

def slice(m,a):

    nc     = m.shape[1]
    nl     = m.shape[0]
    
    mx     = np.arange(nc)
    my     = np.arange(nl)
    
    x,y    = get_limb(mx,my,m)
    p,c    = fit(x,y)

    nsteps = int(180/a)
    n      = max(m.shape[0],m.shape[1])
    xc     = p[0]
    yc     = p[1]
    orig   = yc
    ang    = 0

    r = []
    
    for step in np.arange(nsteps):
        x = (xc - (xc - np.arange(n)) * np.cos(np.radians(ang))).astype('int')
        y = (yc - (yc - np.arange(n)) * np.sin(np.radians(ang))).astype('int')
        x[x < 0 ]  = 0
        x[x >= nc] = nc-1
        y[y < 0 ]  = 0
        y[y >= nl] = nl-1
        profile =  np.ascontiguousarray(m[y,x])

        r.append({'x' : x,
                  'y' : y,
                  'profile' : profile, 'ang' : ang})
            
        ang += a

    return r

