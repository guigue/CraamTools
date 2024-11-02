import numpy as np
import pdb

def contiguo(index):

    nu = len(index)
    if (nu < 2):
        return np.zeros([2,1],int)
    
    n = nu-1
    ad = 0
    
    dif = index[1:] - index[:-1]

#    pdb.set_trace()
    
    if ((dif[0] != 1) & (dif[1] == 1)):
        index = index[1:]
        dif  = index[1:] - index[:-1]
        ad   = 1
        nu   = nu - 1
        n    = n - 1
    
    if ((dif[n-1] != 1) & (dif[n-2] == 1)):
        index = index[0:n-1]
        dif   = index[1:] -index[:-1]
        nu    = nu-1
        n     = n-1
    
    tt, = np.where(dif == 1)
#    pdb.set_trace()    
    if (len(tt) == 0): 
        return np.zero([2,1],int)
        
    if (len(tt) == nu-1):
        pp = np.zeros([2,1],int)
        pp[0,0] = 0
        pp[1,0] = nu-1
        return pp+ad
        
    k, = np.where(dif != 1)
    dif[k] = 0
    dif2 = dif[1:]-dif[:-1]
        
    pi, = np.where(dif2 == -1)
    pf, = np.where(dif2 == 1)
        
    if (pi[0] < pf[0]):
        tmp = np.zeros(len(pi)+1,int)
        for i in range(len(pi)):
            tmp[i+1] = pi[i]
        pi = tmp
        
    if (pi[len(pi)-1] < pf[len(pf)-1]):
        tmp = np.zeros(len(pf)+1,int)
        for i in range(len(pf)):
            tmp[i] = pf[i]
        tmp[len(pf)] = nu-1
        pf = tmp        
            
    pi[1:]  = np.array(pi[1:]) + 2 
    pp      = np.zeros([2,len(pi)],int)
    pp[0,:] = pi
    pp[1,:] = pf
        
    return (pp+ad).T
