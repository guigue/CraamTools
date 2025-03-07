import os
import numpy as np
from astropy import constants as c
from astropy import units as u
import pandas as pd
import pdb

from CraamTools import vertex70 as v
from CraamTools.fit import Gauss

_Version = '2025-03-07T1432BRT'

def hum_to_ppmv(humidity=20,   # %
                temperature=10,# C
                pressure=760): # mbar


    pressure /= 1.01325E+03           # atm
    T0  = 273.15                   # 0 C in Kelvin
    temperature +=T0
    humidity /=100    
    f = T0/temperature             # Temperature ratio
    molecular_weight = 18.01528    # g/mole
    standard_volume  = 22.41383E+03 # cm3 atm / mole    
    r_sat            = f * np.exp(18.9766 - (14.9595+2.43882*f)*f) # Water vapor saturation density

    return r_sat * humidity / molecular_weight * standard_volume / pressure / f * u.Unit('')

def h2o_column(fname='pwv03-mod.tp6'):

    # This function integrates H2O along the atmosphere
    # and returns the precipitable water vapor in mm

    # The input file is the <>.tp6 output from MODTRAN6
    # One must take care about the top lines to throw away.

    # The conversion constant C between atm-cm to pr mm is taken
    # from Wei et al 2019.

    C = 10 * 2.69E+19 / 3.34E+22
    
    f=open(fname,'r')

    h2o = []
    z   = []

    l = ' '
    header_section = 'ATMOSPHERIC PROFILES (AFTER COLUMN SCALING)'
    while (l.rstrip().lstrip() != header_section):
        l=f.readline()

    f.readline()
    f.readline()
    f.readline()
    
    l=f.readline()
    while (l.rstrip().lstrip() != ''):
        t= l.split()
        h2o.append(float(t[3]))
        z.append(float(t[1]))
        l = f.readline()

    h2o = np.asarray(h2o)
    z   = np.asarray(z)
    
    f.close()
    
    return np.trapz(h2o,z) * C * u.Unit('mm'), {'pwv':h2o*C*u.Unit('mm')/u.Unit('km'),'z':z*u.Unit('km')}

        
class HATS(object):

    def __init__(self):
        
        self.version = _Version
        bpf = v.read(os.getenv("HOME")+'/HATS/Documents/CRAAM/Materiales/Archivos/20-TYDEX-BPF15.0-24_01045.0-v70.csv')
        bpf['Transmission']*=bpf['Transmission']

        gaussP,yfit,cov=Gauss.fit1d(bpf['Frequency'],bpf['Transmission'],nt=4)
        hpbw = 2*gaussP[2]*np.sqrt(2*np.log(2))
        self.bpf = {'Gaussian_Parameters':gaussP, 'Gaussian_fit':yfit, 'Gaussian_covariance':cov,
                                 'Filter':bpf}

        return

class AR30T(object):

    def __init__(self):
        
        self.version = _Version
        bpf = pd.read_csv(os.getenv("HOME")+'/IR/documentos/ulis_3elements.csv')
        self.bpf = {'Wavelength':bpf['wl'],'Wavenumber': 1E+04/bpf['wl'],
                    'Frequency': c.c.value*1E-06/bpf['wl'],'Transmission':bpf['T']}
        
            # A Supergaussian with exponent 10 for the IR camera
            #hwidth = 8.5
            #sgp    = int(10)
            #ss     = hwidth**sgp/np.log(2)
            #nu_central = 31.5
            #nu_bottom  = 20
            #nu_upper   = 50
            #x = (self.Data['Frequency'] >= nu_bottom) & (self.Data['Frequency'] <= nu_upper)
            #bpf = np.exp(-(self.Data['Frequency'][x]-nu_central)**sgp/ss)

        return

class MIRI(object):

    def __init__(self):
        
        self.version = _Version
        bpf = pd.read_csv(os.getenv("HOME")+'/IR/QWIP_Response.csv')
        self.bpf = {'Wavelength':bpf['wl'],'Wavenumber': 1E+04/bpf['wl'],
                    'Frequency': c.c.value*1E-06/bpf['wl'],
                    'Transmission':{'blue':bpf['chan_blue'],'red':bpf['chan_red']}}
        
        return
    
class model(object):

    def __init__(self,root):
        
        self.version = _Version

        fname = root+'.tp7'

        try:
            fh = open(fname)
            self.read(fh)
            self.File_Name = fname
        except:
            print("\n\n File "+fname+"  not found\n\n")

        self.HATS = HATS()
        self.AR30T = AR30T()
        self.MIRI = MIRI()
        
        return

    def compute_optical_depth(self,instrument='hats',type='real'):

        self.Transmission = {'Instrument':instrument,'Type':type}

        if (instrument.lower() == 'hats'):
            if (type.lower() == 'model'):
                nu_bottom = 10
            elif (type.lower() == 'real'):
                nu_bottom = self.HATS.bpf['Filter']['Frequency'][self.HATS.bpf['Filter']['Frequency'].shape[0]-1]
            nu_upper  = 20
            x = (self.Data['Frequency'] >= nu_bottom) & (self.Data['Frequency'] <= nu_upper)

            if (type.lower() == 'real'):
                bpf = np.interp(self.Data['Frequency'][x],np.flip(self.HATS.bpf['Filter']['Frequency']),np.flip(self.HATS.bpf['Filter']['Transmission']))
            elif (type.lower() == 'model'):
                bpf = np.exp(-(self.Data['Frequency'][x]-self.HATS.bpf['Gaussian_Parameters'][1])**2/(2*self.HATS.bpf['Gaussian_Parameters'][2]**2))

        elif (instrument.lower() == 'ar30t'):
            nu_bottom = self.AR30T.bpf['Frequency'].min()
            nu_upper = self.AR30T.bpf['Frequency'].max()
            x = (self.Data['Frequency'] >= nu_bottom) & (self.Data['Frequency'] <= nu_upper)
            bpf = np.interp(self.Data['Frequency'][x],np.flip(self.AR30T.bpf['Frequency']),np.flip(self.AR30T.bpf['Transmission']))

        elif (instrument.lower() == 'miri'):
            nu_bottom = self.MIRI.bpf['Frequency'].min()
            nu_upper = self.MIRI.bpf['Frequency'].max()
            x = (self.Data['Frequency'] >= nu_bottom) & (self.Data['Frequency'] <= nu_upper)
            bpf = {}
            bpf.update({'red':np.interp(self.Data['Frequency'][x],np.flip(self.MIRI.bpf['Frequency']),np.flip(self.MIRI.bpf['Transmission']['red']))})
            bpf.update({'blue':np.interp(self.Data['Frequency'][x],np.flip(self.MIRI.bpf['Frequency']),np.flip(self.MIRI.bpf['Transmission']['blue']))})
            
        else:
            return
        


        if (instrument.lower() == 'miri'):
            self.Transmission.update({'bpf':bpf,'Frequency':self.Data['Frequency'][x],'Value':{'red':0.0, 'blue':0.0}})
            self.Transmission['Value']['red']  = np.trapz(self.Data['Transmission'][x]*bpf['red'],self.Data['Frequency'][x])  / np.trapz(bpf['red'],self.Data['Frequency'][x] )
            self.Transmission['Value']['blue'] = np.trapz(self.Data['Transmission'][x]*bpf['blue'],self.Data['Frequency'][x]) / np.trapz(bpf['blue'],self.Data['Frequency'][x])

        else:
            self.Transmission.update({'bpf':bpf,'Frequency':self.Data['Frequency'][x],'Value':0.0})
            self.Transmission['Value'] = np.trapz(self.Data['Transmission'][x]*bpf,self.Data['Frequency'][x]) / np.trapz(bpf,self.Data['Frequency'][x])

        return
        
    def read(self,fh):

        for i in np.arange(11):
            fh.readline()

        wn = []
        wl = []
        nu = []
        tr = []
        od = []
        br = []

        for line in fh:

            token = line.strip().split()
            if len(token) > 1 :
                wn.append(float(token[0]))                          # wavenumber [cm-1]
                wl.append(1.0E+04 / float(token[0]))                # wavelength [um]
                nu.append(c.c.value * float(token[0]) * 1.0E-10)    # Frequency [THz]
                tr.append(float(token[1]))                          # Transmission
                if (float(token[1]) > 0):
                    od.append(-np.log(float(token[1])))
                else:
                    od.append(100)

                br.append(float(token[14]))                         # Solar Brightness top of the Atmosphere

        fh.close()

        self.Data = np.zeros(len(wn),dtype=[('Wavenumber',np.float64), ('Wavelength',np.float64)   ,
                                            ('Frequency',np.float64) , ('Transmission',np.float64) ,
                                            ('Opt_Depth',np.float64) , ('Brightness',np.float64)   ] )

        self.Data['Wavenumber']   = np.asarray(wn)
        self.Data['Wavelength']   = np.asarray(wl)
        self.Data['Frequency']    = np.asarray(nu)
        self.Data['Transmission'] = np.asarray(tr)
        self.Data['Opt_Depth']    = np.asarray(od)
        self.Data['Brightness']   = np.asarray(br)

        return
    
            

    
