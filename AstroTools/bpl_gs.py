#########################################################################################################
#
# bpl_gs : A Bi-Power Law electron distribution and its synchrotron emission.
#
# Aim: To compute the synchrotron emission from an electron distribution described by
#      a double (bi) power-law at a 1 AU distance from the observer.
#      It is based on the original Ramaty's FORTRAN code for a single power-law
#      This code was translated to IDL and added the second power law. You may
#      the original code in $HOME/solar/RadioFlux/bplgs.pro
#
# Usage:
#      >>> import bpl_gs
#      >>> b=bpl_gs.BrokenPowerLawGS(
#                 Ntotal=1.0E+30,                      # Total electron number
#                 delta=np.array([3.0,4.0]),           # delta array, 2 elements needed
#                 Edges=np.array([0.02,0.5,30.]),      # Energy edges: low cutoff, brea-energy, high cutoff [MeV]
#                 bmag=1.0E+03,                        # Magnetic Field intensity
#                 viewangle=45.0,                      # Viewing angle [degrees]
#                 size=14E+08,                         # Source size [cm]
#                 height=1.0E+08,                      # Source height [cm]
#                 kf=121,                              # Number of harmonics in calculation
#                 etr=3.0,                             # Transition energy between full gyrosynchrotron calculation and synchrotron aproximation [MeV]
#                 alpha=1.5,                           # Ramaty's alpha factor for the razin effect (see his 1969 paper)
#                 distance=1.0)                        # distance in AU
#
# Output:
#     'b' is an object with two elements 'version' and 'flux'. The latter has the main
#     outputs: fghz (frequency array in GHz), jy (total flux in Jy) q,v (Stokes Q,V)
#
# Author:  Guigue - On the day the World became too much crazy...
#
#######################################################################################################

import numpy as np
from scipy import special as sp
import pdb                                       # add pdb.set_trace() in the line where you want to stop
from datetime import datetime
import astropy.constants as cts
from CraamTools.AstroTools import constants as ct

class Flux:

    def __init__(self,source):

        e1d  = np.zeros(source.kf) ; e2d  = np.zeros(source.kf) ; a1d  = np.zeros(source.kf) ; a2d  = np.zeros(source.kf)
        rccd = np.zeros(source.kf) ; ffbd = np.zeros(source.kf) ; phi1 = np.zeros(source.kf) ; phi2 = np.zeros(source.kf)
        phit = np.zeros(source.kf) ; rcd  = np.zeros(source.kf) ; q    = np.zeros(source.kf) ; v    = np.zeros(source.kf)
        scd  = np.zeros(source.kf) ; ff   = np.zeros(source.kf)

# -----------------Constants -------------------------------------------------------------

        m0     = cts.m_e.cgs.value                 # electron mass [g]
        c      = cts.c.cgs.value                   # [cm/s] speed of light
        e      = cts.e.esu.value                   # electron charge
        em_fac = e**3 / (m0 * c**2)                #  Emissivity factor
        ab_fac = 4 * np.pi**2 * e                  # Absorption factor
        E0     = (cts.m_e*cts.c**2).to('MeV')      # electron rest energy [MeV]
        nuB    = 0.5 / np.pi * e / (m0 * c)        # x bmag = Fundamental frequency
        au     = cts.au.cgs.value                  # Astronomical Unit [cm]
        au2    = au * au
        erg2Jy = 1.0E+23                           # erg/cm2/s -> Jy

# ----------------------------------------------------------------------------------------

        d     = source.distance * au                 # distance in cm
        d2    = d*d                                # square of the distance
        cs    = np.cos(np.radians(source.viewangle))
        ss    = np.sin(np.radians(source.viewangle))
        if (cs < 0.1) or (cs > 0.95):
            print('  ')
            print('Viewangle is beyond limits')
            print('  ')
            return

        area  = (source.size / 2)**2 * np.pi         # source area
        SolAng= area / d2                          # source solid angle
        vol   = (source.size / 2 )**2 * np.pi * source.height
        dens  = source.Ntotal / vol
        ffp   = 1.5 / source.alpha
        gtr   = source.etr / E0.value + 1.0

#----------------------------------------------------------------------------------------

        if (np.shape(source.delta)[0] != 2) and (np.shape(Edges)[0] != 3):
            print('   ')
            print('delta should have 2 elements and Energies 3.')
            print('Please correct input and run again.')
            print('   ')
            return


        t1 = datetime.now()

        gEdges = source.Edges / E0.value + 1.0

# ------------compute the normalized electron distribution factors (a.k.a. "anors")----------
        ce1    = np.array([source.Edges[0]**(1.0-source.delta[0]),source.Edges[1]**(1.0-source.delta[1])])
        ce2    = np.array([source.Edges[1]**(1.0-source.delta[0]),source.Edges[2]**(1.0-source.delta[1])])
        ce     = (ce1-ce2)/(source.delta-1.0)
        f      = [1.0,source.Edges[1]**(source.delta[1]-source.delta[0])]
        aa     = f * ce
        anor   = 1.0 / np.sum(aa) * np.array([1.0 , source.Edges[1]**(source.delta[1]-source.delta[0])])

#------------compute the j1 and j2 Ramatys integer energy limits----------------------------
#            In the fortran code, j1 and j2 are inputs
        j1 = int(( 2.0 + np.log10(source.Edges[0]) ) * 10.0)
        j2 = int(( 2.0 + np.log10(source.Edges[2]) ) * 10.0)

# ------------start loop on ffb-------------------------------
        k = 0
        while (k < source.kf):
            if (k <= 100):
                ffb = 10**(0.01*k)
            else:
                ffb = 10**(0.1*(k-90))
            ffbd[k]=ffb
            an1,an2,ath1,ath2 = self.refr(ffb,ffp,cs)
            e1=0.0 ; e2=0.0 ; a1=0.0 ; a2=0.0

# ----------start integration over energy--------------------
            j = j1
            while (j <= j2):

                el    = 10**(0.1*(j-1)-2)
                eu    = 10**(0.1*j-2)
                em    = 10**(0.1*(j-0.5)-2)
                gamma = em/E0.value + 1.0
                g1    = 0.0
                g2    = 0.0

                if (gamma < gtr):
                    if (ffb > ffp):
                        an = np.sqrt(an1)
                        g1 = self.gsy(gamma,ffb,cs,an,ath1)

                    if (ffb > (np.sqrt(ffp**2+0.25)+0.5)):
                        an = np.sqrt(an2)
                        g2= self.gsy(gamma,ffb,cs,an,ath2)
                else:
                    if (ffb > (np.sqrt(ffp**2+0.25)+0.5)):
                        g1,g2 = self.ssy(gamma,ffb,cs,source.alpha)

#------------compute the normalized number of electrons in the energy band --------------------
                if ( (el < source.Edges[1]) and (eu > source.Edges[1]) ):
                    Nel = anor[0] * (source.Edges[1]**(1.0-source.delta[0]) - el**(1.0-source.delta[0])) / (1.0-source.delta[0]) + \
                          anor[1] * (eu**(1.0-source.delta[1]) - source.Edges[1]**(1.0-source.delta[1])) / (1.0-source.delta[1])
                else:
                    if (el < source.Edges[1]):
                        interval = 0
                    else:
                        interval = 1
                    ex = source.delta[interval]
                    Nel = anor[interval] * (eu**(1-ex) - el**(1-ex)) / (1-ex)

#-------------------------------------------------------------------------------------------------
                e1 = e1 + g1 * Nel * source.bmag * em_fac
                e2 = e2 + g2 * Nel * source.bmag * em_fac
                a1 = a1 + g1 * Nel / ffb / ffb / source.bmag * ab_fac * \
                    (ex * gamma * (gamma + 1.0) + 2.0 * gamma**2 - 1.0) / gamma / (gamma**2 - 1.0)
                a2 = a2 + g2 * Nel / ffb / ffb / source.bmag * ab_fac * \
                    (ex * gamma * (gamma + 1.0) + 2.0 * gamma**2 - 1.0) / gamma / (gamma**2 - 1.0)

                j=j+1


#--------end of integration over energy----------------------
            a1 = a1 * dens
            a2 = a2 * dens
            e1 = e1 * dens
            e2 = e2 * dens

            e1d[k] = e1
            e2d[k] = e2
            a1d[k] = a1
            a2d[k] = a2

            if (e2 > 0.0):
                rccd[k] = (e2-e1) / (e2+e1)
            else:
                rccd[k] = 0.0

            arg1 = source.height * a1
            if (arg1 < 1.0E-03):
                phi1[k] = e1 * vol / d2
            else:
                phi1[k] = SolAng * e1 / a1 * (1.0 - np.exp(-arg1))

            arg2 = source.height * a2
            if (arg2 < 1.0E-03):
                phi2[k] = e2 * vol / d2
            else:
                phi2[k] = SolAng * e2 / a2 * (1.0 - np.exp(-arg2))

#-----calculate Stokes parameters------------------
            phit[k] = phi1[k] + phi2[k]
            q[k]    = phi1[k] * (1.0 - ath1**2) / (1.0 + ath1**2) + \
                phi2[k] * (1.0 - ath2**2) / (1.0 + ath2**2)
            v[k]    = 2.0 * (phi1[k] * ath1 / (1.0 + ath1**2) + phi2[k] * ath2 / (1.0 + ath2**2))
#-----calculate polarization-------------------
            if (phit[k] > 0):
                rcd[k]  = (phi2[k] - phi1[k] ) / phit[k]
                scd[k]  = v[k] / abs(v[k]) * np.sqrt(q[k]**2 + v[k]**2) / phit[k]
            else:
                rcd[k] = 0.0
                scd[k] = 0.0

            ff[k]    = ffbd[k] * nuB * source.bmag
            k = k + 1

# ------------end of loop over ffb----------------------------

        self.fghz = ff / 1.0E+09
        self.f_o = phi1 * erg2Jy
        self.f_x = phi2 * erg2Jy
        self.jy  = phit * erg2Jy
        self.a_o = a1d
        self.a_x = a2d
        self.e_o = e1d
        self.e_x = e2d
        self.rc  = rcd
        self.q   = q
        self.v   = v

        t2 = datetime.now()


        print("  ")
        print("Total time elapsed = ",t2-t1)
        print("-----------------------------------------")


        return

    def bespr(self,n,x):
        b1 = sp.jv(n+1,x)
        b  = sp.jv(n,x)
        bpr = -b1 + n / x * b
        return bpr

    def refr(self,ffb,ffp,cs):

        ss    = np.sqrt(1 - cs*cs)
        diff4 = (ffp*ffp - ffb*ffb)

        anum  = 2 * ffp*ffp * (ffp*ffp - ffb*ffb)
        dnum1 = np.sqrt(ffb**4 * ss**4 + 4 * ffb*ffb * diff4*diff4 * cs*cs) - \
            2 * ffb*ffb * diff4 - ffb*ffb * ss*ss
        dnum2 =-np.sqrt(ffb**4 * ss**4 + 4 * ffb*ffb * diff4*diff4 * cs*cs) - \
            2 * ffb*ffb * diff4 - ffb*ffb * ss*ss

        if (dnum1 == 0):
            an1 = 1.0
        else:
            an1   = 1 + anum / dnum1
        if (dnum2 == 0):
            an2 = 1.0
        else:
            an2   = 1 + anum / dnum2
        aknum = 2 * ffb * diff4 * cs
        dknum1= np.sqrt(ffb**4 * ss**4 + 4 * ffb*ffb * diff4*diff4 * cs*cs) - ffb*ffb * ss*ss
        dknum2=-np.sqrt(ffb**4 * ss**4 + 4 * ffb*ffb * diff4*diff4 * cs*cs) - ffb*ffb * ss*ss
        ath1  =-aknum / dknum1
        ath2  =-aknum / dknum2
        return an1, an2, ath1, ath2

    def ssy(self,gamma,ffb,cs,alpha):

        xd = np.array([0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, \
                       0.300, 0.400, 0.50, 0.600, 0.70, 0.800, 0.9, 1.00, 1.2, 1.40, \
                       1.600, 1.800, 2.00, 2.500, 3.00, 3.500, 4.0, 4.50, 5.0, 6.00, \
                       7.0, 8.0, 9.0, 10.0])

        fd = np.array([0.213 , 0.358 , 0.445 , 0.583 , 0.702 , 0.722 , 0.818   , 0.874, 0.904,\
                       0.917 , 0.919 , 0.901 , 0.872 , 0.832 , 0.788 , 0.742   , 0.694, 0.655,\
                       0.566 , 0.486 , 0.414 , 0.354 , 0.301 , 0.200 , 0.130   , 0.0845, 0.0541,\
                       0.0339, 0.0214, 0.0085, 0.0033, 0.0013, 0.0005, 0.00019 ] )

        ss  = np.sqrt(1 - cs*cs)

        ffc = (2/3) * ffb / (ss * gamma * gamma) * (1 + (9/4) * (gamma*gamma - 1) / (alpha**2 * ffb * ffb) )**1.5

        if (ffc <= xd[0]):
            f = 2.15 * ffc**(1/3) * (1 - 0.844 * ffc**(2/3) )
        if (ffc > xd[33]):
            f = 1.253 * np.exp(-ffc) * np.sqrt(ffc) * (1 + 55 / (72 * ffc) )
        if (ffc > xd[0] and ffc <= xd[33]):
            #            The following sentence is equivalent to IDL i = min(where(ffc le xd))
            i = np.where(ffc <= xd)[0][0]
            f = ( fd[i] * (ffc - xd[i-1]) + fd[i-1] * (xd[i] - ffc) )  / ( xd[i] - xd[i-1] )

        gtot = 0.138 * ss * f
        g1   = gtot / 2

        return g1,g1

    def gsy(self,gamma,ffb,cs,an,ath):

        beta = np.sqrt(gamma*gamma-1) / gamma
        ss   = np.sqrt(1 - cs*cs)
        ffc  = ffb * 2 / ( 3 * ss * gamma * gamma)
        if (ffc >= 20.0):
            g12 = 0.0
            return g12
        else:
            is1 = int(ffb * gamma * (1 - an * beta * cs) + 1)
            is2 = int(ffb * gamma * (1 + an * beta * cs))
            sum12 = 0.0

        i = is1
        while (i <= is2):
            cphis  = (1 - i / (ffb * gamma) ) / (beta * cs * an)
            sphis  = np.sqrt(1.0 - cphis*cphis)
            xs     = i * an * beta * ss * sphis / (1.0 - an * beta * cs * cphis)
            xstr   = xs / i

            if (ffb > 250 ):
                xstr = 0.9
                if (gamma > 5):
                    xstr = 0.92
                if (gamma > 10):
                    xstr = 0.97
                if (gamma > 15):
                    xstr = 0.98

            if (ffb > 50):
                xstr = 0.8
                if (gamma > 5):
                    xstr = 0.9
                if (gamma > 10):
                    xstr = 0.95
                if (gamma > 15):
                    xstr = 0.96

            if (xs >= xstr*i):
                b = sp.jv(i,xs)
                bpr=self.bespr(i,xs)

                f12    = (-beta * sphis * bpr + ath * (cs / ss / an - beta * cphis / ss) * b)**2
                s12old = sum12
                sum12  = sum12 + f12
                if ( (s12old > 0) and ( (sum12-s12old)/s12old < 1.0e-04) ):
                    break

            i = i+1

        g12 = sum12 / (2 * beta * cs) * ffb / (1 + ath*ath)
        return g12

class BrokenPowerLawGS:

    def __init__(self,
                 Ntotal=1.0E+30,
                 delta=np.array([3.0,4.0]),
                 Edges=np.array([0.02,0.5,30.]),
                 bmag=1.0E+03,
                 viewangle=45.0,
                 size=14E+08,
                 height=1.0E+08,
                 kf=121,
                 etr=3.0,
                 alpha=1.5,
                 distance=1.0):

        self.version   = '2020-05-20T1825BST'
        self.Ntotal    = Ntotal
        self.delta     = delta
        self.Edges     = Edges
        self.bmag      = bmag
        self.viewangle = viewangle
        self.size      = size
        self.height    = height
        self.kf        = kf
        self.etr       = etr
        self.alpha     = alpha
        self.distance  = distance
        self.nuB       = 0.5 * cts.e.esu * self.bmag /(np.pi * cts.m_e.cgs * cts.c.cgs)
        self.nuP       = 1.5 * self.nuB / self.alpha
        self.ntherm    = self.nuP**2 * cts.m_e.cgs * np.pi / cts.e.esu**2


        self.flux      = Flux(self)

        return
