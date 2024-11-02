import numpy as np
import pandas as pd

def read_out(fname='amc.out'):

    fq = []
    tau = []
    ts = []

    with open(fname,'r') as file:
        lines = file.readlines()
        for line in lines:
            c=line.strip().split()
            fq.append(float(c[0]))
            tau.append(float(c[1]))
            ts.append(float(c[2]))

    return {'frequency':np.asarray(fq), 'tau':np.asarray(tau), 'Ts': np.asarray(ts)}

def integrate(tau,frange,type='tau'):
    
    x = ( (tau['frequency'] >= frange[0]) & (tau['frequency'] <= frange[1]) )

    if (type == 'tau'):
        return np.trapz(tau['tau'][x],tau['frequency'][x])/(frange[1]-frange[0])
    else:
        return np.trapz(tau['Ts'][x],tau['frequency'][x])/(frange[1]-frange[0])
