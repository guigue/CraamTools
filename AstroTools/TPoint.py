import numpy as np
from matplotlib import pyplot as plt

lista = ["HATS-OpT-TP-2025-12-08_in.dat",
         "HATS-OpT-TP-2025-12-09_in.dat",
         "HATS-OpT-TP-2025-12-10_in.dat",
         "HATS-OpT-TP-2025-12-12_in.dat",
         "HATS-OpT-TP-2025-12-13_in.dat",
         "HATS-OpT-TP-2025-12-15_in.dat",
         "HATS-OpT-TP-2026-01-09_in.dat",
         "HATS-OpT-TP-2026-01-10_in.dat",
         "HATS-OpT-TP-2026-01-11_in.dat",
         "HATS-OpT-TP-2026-01-12_in.dat",
         "HATS-OpT-TP-2026-01-13_in.dat",
         "HATS-OpT-TP-2026-01-14_in.dat",
         "HATS-OpT-TP-2026-01-15_in.dat",
         "HATS-OpT-TP-2026-01-16_in.dat",
         "HATS-OpT-TP-2026-01-17_in.dat"]

def ascii2sec(d,m,s):
    
    arcdeg = float(d)
    arcmin = float(m)
    arcsec = float(s)
    sec    = (np.abs(arcdeg)*3600 + arcmin*60 + arcsec)
    
    if arcdeg < 0:
        return -sec
    else:
        return sec
    
def read(lnames=lista, telescope='HATS'):

    ephra  = []
    ephdec = []
    telra  = []
    teldec = []
    lsid   = []
    
    for fname in lnames:
        with open(fname,'r',encoding='ascii') as f:
            if telescope.upper() == 'HATS':
                for n in np.arange(5):
                    z=f.readline()
                for l in f:
                    s=l.split()
                    ephra.append((ascii2sec(s[0],s[1],s[2]))*15)
                    dec = ascii2sec(s[3],s[4],s[5])
                    ephdec.append(dec)
                    telra.append((ascii2sec(s[6],s[7],s[8]))*15)
                    dec = ascii2sec(s[9],s[10],s[11])
                    teldec.append(dec)
                    lsid.append(float(s[12])+float(s[13])/60)

    return {'EphRA' : np.asarray(ephra)  ,
            'EphDec': np.asarray(ephdec) ,
            'TelRA' : np.asarray(telra)  ,
            'TelDec': np.asarray(teldec) ,
            'lSID'  : np.asarray(lsid)
            }


def plot(c,savefig=False):

    hCirc = 180*3600

    dRA  = (c['EphRA'] - c['TelRA'])
    dDec = (c['EphDec'] - c['TelDec']) 
    
    fig  = plt.figure()
    fig.set_size_inches(w=10,h=6)
    Bpos = [0.1,0.1,0.88,0.88]
    ax   = fig.add_subplot(1,1,1, position=Bpos)
    ax.tick_params(axis='both', labelsize=14)     
    ax.plot(dRA,dDec,'o',color='red')
    ax.set_xlabel('RA [arc sec]',fontsize=14)
    ax.set_ylabel('Dec [arc sec]',fontsize=14)

    if savefig:
        fname = 'OpT_Errors.pdf'
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w',
                    orientation='portrait', format='pdf',
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    metadata=None)

    return
    


            
