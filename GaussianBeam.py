import numpy as np

def sigma(hpbw):
    return hpbw/np.sqrt(np.log(256))

def hpbw(sigma):
    return sigma * np.sqrt(np.log(256))


