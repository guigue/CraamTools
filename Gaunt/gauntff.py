import sys, string, os, struct, glob
import numpy as np
from scipy.interpolate import interp2d
import pdb

_Version_ = '20211130T1953BRT'                                   #

#################################################################
#
# gauntff
#
#   Following van Hoof et al (2014), doi:10.1093/mnras/stu1438
#   this class implements the method to interpolate its averaged
#   gaunt factor table for any pair (T,nu).
#
#   It is a port of the IDL  rd_gauntff.pro written by Paulo Simões.
#   Difference between both algorithms are < 1% and due to the
#   differences in the interpolation algorithms in IDL and SCIPY.
#
#   use:
#     >>> import gauntff as g
#     >>> T=2.0E+06  # temperature in K
#     >>> nu=2.0E+10 # frequency in Hz
#     >>> gaunt=g.factor(T,nu)
#     >>> gaunt.compute()
#     >>> print("Gaunt Factor = {0:f}".format(gaunt.value))
#
#   Optional Parameters:
#         filename='filename' : file name of van Hoof ASCII table
#                               default is gauntff.dat, see its location in the system.
#
#   Todo:
#         Implement the other tables as well...
#
#   Author:  @guiguesp on 2021-11-30
#            Adapted from rd_gauntff.pro written by P. Simões
#
##################################################################

class factor(object):

    def __init__(self,T,nu,filename='/monalisa/guigue/Programming/python/CraamTools/Gaunt/gauntff.dat'):
        self.T = T
        self.nu = nu
        self.filename = filename
        return

    @property
    def version(self):
        return self._Version

    def __str__(self):
        return "A Class representing the van Hof's gaunt factor "
    
    def read(self):

        f = open(self.filename,'r')
        for i in np.arange(42):
            n=f.readline()

        gff = np.zeros([81,146],dtype=np.float64)
    
        for i in np.arange(146):
            l=f.readline()
            x=[float(x) for x in l.split()]
            gff[:,i] = x

        f.close()
        self.gff = gff
        return

    def compute(self,verbose=False):

        self.read()
        
        z  = 1.0
        hh = 6.62607551E-27
        kk = 1.38006504E-16
        ry = 2.17987E-11

        dex   = 2.0E-01
        ngd   = 81
        gam20 = -6.0
        g     = np.linspace(gam20,ngd*dex+gam20,ngd)

        nud   = 146
        u0    = -16.0
        u     = np.linspace(u0,nud*dex+u0,nud)

        find_u = np.log10(hh * self.nu / (kk * self.T) )
        find_g = np.log10(z**2 * ry / (kk * self.T) )
        
        if verbose:
            print(" find_u = {0:e}   find_g = {1:e}".format(10.**find_g, 10**find_u))

        self.gffi =  interp2d(u,g,self.gff)
        self.value = self.gffi(find_u,find_g)[0]
        return
    


