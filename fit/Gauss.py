#######################################################################
#
#  Gauss
#
#     A fitting function for Gaussian shapes
#     Adapted from IDL routine GAUSSFIT.PRO
#     There is a small difference with IDL for fittings with number of terms nt > 4.
#     IDL uses a 1st degree polynomial fitting to subtract from the raw data.
#     For some reason, this does not work in python. Gauss.py subtract the mean
#     of the raw data.
#
#  Usage:
#     par,fit,cov = Gauss.fit1d(x,y,[nt=nt])
#         where nt is the number of terms, defaults is 3, maximum is 6
#         returns
#         par : fitting obtained parameters
#               a[0] : amplitude of Gaussian
#               a[1] : center of Gaussian
#               a[2] : standard deviation
#               a[3] : constant
#               a[4] : linear parameter
#               a[5] : quadratic parameter
#  Fitting function:
#     F(x) = a[0] EXP((x-a[1])^2/a[2]^2/2) + a[3] + a[4]*x + a[5]*x^2
#
#  Author: @guiguesp
#  History: started one day 7 (!) years ago.
#  Last modification: 2024-10-16
#     
#     Follows original comments of IDL
#---------------------------------------------------------------
#+
# NAME:
#   GAUSSFIT
#
# PURPOSE:
#   Fit the equation y=f(x) where:
#
#       F(x) = A0*EXP(-z^2/2) + A3 + A4*x + A5*x^2
#           and
#       z=(x-A1)/A2
#
#   A0 = height of exp, A1 = center of exp, A2 = sigma (the width).
#   A3 = constant term, A4 = linear term, A5 = quadratic term.
#   Terms A3, A4, and A5 are optional.
#   The parameters A0, A1, A2, A3 are estimated and then CURVEFIT is
#   called.
#
# CATEGORY:
#   ?? - fitting
#
# CALLING SEQUENCE:
#   Result = GAUSSFIT(X, Y [, A])
#
# INPUTS:
#   X:  The independent variable.  X must be a vector.
#   Y:  The dependent variable.  Y must have the same number of points
#       as X.
#
# KEYWORD INPUTS:
#
#   CHISQ: Set this keyword to a named variable that will contain
#      the value of the chi-square goodness-of-fit.
#
#   ESTIMATES = optional starting estimates for the parameters of the
#       equation.  Should contain NTERMS (6 if NTERMS is not
#       provided) elements.
#
#   MEASURE_ERRORS: Set this keyword to a vector containing standard
#       measurement errors for each point Y[i].  This vector must be the same
#       length as X and Y.
#
#     Note - For Gaussian errors (e.g. instrumental uncertainties),
#        MEASURE_ERRORS should be set to the standard
#        deviations of each point in Y. For Poisson or statistical weighting
#        MEASURE_ERRORS should be set to sqrt(Y).
#
#   NTERMS = Set NTERMS to 3 to compute the fit: F(x) = A0*EXP(-z^2/2).
#      Set it to 4 to fit:  F(x) = A0*EXP(-z^2/2) + A3
#      Set it to 5 to fit:  F(x) = A0*EXP(-z^2/2) + A3 + A4*x
#
#   SIGMA: Set this keyword to a named variable that will contain
#      the 1-sigma error estimates of the returned parameters.
#
#     Note: if MEASURE_ERRORS is omitted, then you are assuming that
#           your model is correct. In this case, SIGMA is multiplied
#           by SQRT(CHISQ/(N-M)), where N is the number of points
#           in X and M is the number of terms in the fitting function.
#           See section 15.2 of Numerical Recipes in C (2nd ed) for details.
#
#   YERROR: The standard error between YFIT and Y.
#
# OUTPUTS:
#   The fitted function is returned.
#
# OPTIONAL OUTPUT PARAMETERS:
#   A:  The coefficients of the fit.  A is a three to six
#       element vector as described under PURPOSE.
#
# COMMON BLOCKS:
#   None.
#
# SIDE EFFECTS:
#   None.
#
# RESTRICTIONS:
#   The peak or minimum of the Gaussian must be the largest
#   or smallest point in the Y vector.
#
# PROCEDURE:
#   The initial estimates are either calculated by the below procedure
#   or passed in by the caller.  Then the function CURVEFIT is called
#   to find the least-square fit of the gaussian to the data.
#
#  Initial estimate calculation:
#   If NTERMS>=4 then a constant term is subtracted first.
#   If NTERMS>=5 then a linear term is subtracted first.
#   If the (MAX-AVG) of Y is larger than (AVG-MIN) then it is assumed
#   that the line is an emission line, otherwise it is assumed there
#   is an absorbtion line.  The estimated center is the MAX or MIN
#   element.  The height is (MAX-AVG) or (AVG-MIN) respectively.
#   The width is found by searching out from the extrema until
#   a point is found less than the 1/e value.
#
# MODIFICATION HISTORY:
#   DMS, RSI, Dec, 1983.
#   DMS, RSI, Jun, 1995, Added NTERMS keyword.  Result is now float if
#               Y is not double.
#   DMS, RSI, Added ESTIMATES keyword.
#   CT, RSI, Feb 2001: Change the way estimates are computed.
#         If NTERMS>3 then a polynomial of degree NTERMS-4 is subtracted
#         before estimating Gaussian coefficients.
#   CT, RSI, Nov 2001: Slight change to above modification:
#         Because a Gaussian and a quadratic can be highly correlated,
#         do not subtract off the quadratic term,
#         only the constant and linear terms.
#         Also added CHISQ, SIGMA and YERROR output keywords.
#   CT, RSI, May 2003: Added MEASURE_ERRORS keyword.
#   CT, RSI, March 2004: If estimate[2] is zero, compute a default value.
#   CT, ITTVIS, Sept 2008: Do all computations in double precision,
#       convert back to single precision if inputs were single precision.
#-
#

# External methods
import numpy as np
from scipy.optimize import curve_fit
import pdb

__version__ = '2024-10-16T12:06BRT'

def version():
    return __version__

def fit1d_func(x,*a):

    nx = x.shape
    n = len(a)

    if a[2] != 0 :
        z = (x-a[1])/a[2]
        ez = np.exp(-z*z/2)
    else:
        ez = np.zeros(nx)

    if n == 3:
        return a[0]*ez
    elif n == 4:
        return a[0]*ez + a[3]
    elif n == 5:
        return a[0]*ez + a[3] + a[4]*x
    elif n == 6:
        return a[0]*ez + a[3] + a[4]*x + a[5]*x*x

def fit1d(x,y,nt=3,est=[]):

    nEst = len(est)
    n    = y.shape[0]
    
    if nEst == 0:
        if (nt >= 4):
            c = y.mean()
            yd = y - c
        else:
            yd = y
            c  = 0.0

        ymax = np.max(yd)
        i0   = np.argmax(yd)
        dy   = yd[i0]
        dx   = dy/np.exp(1)
        i    = 0
        while ((i0+i+1) < n) & ((i0-i) > 0) & (np.abs(yd[i0+i]) > np.abs(dx)) & (np.abs(yd[i0-i]) > np.abs(dx)):
            i+=1

        if (nt == 3):
            est = [ymax,x[i0],np.abs(x[i0]-x[i0+i])]
        elif (nt == 4):
            est = [ymax,x[i0],np.abs(x[i0]-x[i0+i]),c]
        elif (nt == 5):
            est = [ymax,x[i0],np.abs(x[i0]-x[i0+i]),c,c]
        elif (nt == 6):
            est = [ymax,x[i0],np.abs(x[i0]-x[i0+i]),c,c,0]

    par, cov = curve_fit(fit1d_func,x,y,p0=est)
    yfit     = fit1d_func(x,*par)
    
    return par,yfit,cov


def fit2d(z, x=None, y=None, a=None, negative=False, tilt=False):
    
    s=z.shape
    if len(s) != 2 :
        print('\n  ')
        print('z must have two dimensions\n')
        return False

    nx = s[1]
    ny = s[0]
    n  = nx*ny

    if x:
        if nx != len(x):
            print('\n  ')
            print('x array must have size equal to number of columns of z\n')
            return
    else:
        x=np.linspace(0,nx-1,nx)

    if y:
        if ny != len(y):
            print('\n  ')
            print ('\n y array must have size equal to number of rows of z')
            return
    else:
        y=np.linspace(0,ny-1,ny)
        
    weights = np.ones(n)

    if negative:
        i=np.argmin(z)
    else:
        i=np.argmax(z)

    ix = i % nx
    iy = i // nx
    x0 = int(x[ix])
    y0 = int(y[iy])

    if not a:
        parx,yfitx,covx = fit1d(x,z[y0,:],nt=4)
        pary,yfity,covy = fit1d(y,z[:,x0],nt=4)

        a = [ (parx[3]+pary[3])/2, np.sqrt(np.abs(parx[0]*pary[0])), parx[2], pary[2], parx[1], pary[1] ]
        if tilt:
            a.append(0.0)

        a=np.asarray(a)
        
    par, cov = curve_fit(fit2d_func,[x, y],np.ravel(z),p0=a)
    zfit     = (fit2d_func([x,y],*par)).reshape(nx,ny)

    return par,zfit,cov

def fit2d_func(x,*a):
    nx=x[0].shape[0]
    ny=x[1].shape[0]
    xx=x[0]
    yy=x[1]

    tilt = (len(a) == 7)
    if tilt :		      # Rotate?
        xp = np.matmul ( (np.ones(ny)).reshape(ny,1), (xx-a[4]).reshape(1,nx) )  # Expand X values
        yp = np.matmul ( (yy-a[5]).reshape(ny,1), (np.ones(nx)).reshape(1,nx) )	 # Expand Y values
        s = np.sin(a[6])
        c = np.cos(a[6])
        t =  xp * c/a[2] - yp * s/a[2]
        yp = xp * s/a[3] + yp * c/a[3]
        xp = t
    else:
        xp = np.matmul ( (np.ones(ny)).reshape(ny,1)  , (xx-a[4]).reshape(1,nx)/a[2] )   # Expand X values
        yp = np.matmul ( (yy-a[5]).reshape(ny,1)/a[3] , (np.ones(nx)).reshape(1,nx)  )	 # Expand Y values
        s = 0.0
        c = 1.0

    n = nx * ny
    u = np.ravel(np.exp(-0.5 * (xp*xp + yp*yp)))	                                 # Exp() term, Make it 1D

    return a[0] + a[1] * u

    
        
