import numpy as np
from collections import deque,Counter
from bisect import insort, bisect_left
from itertools import islice
import pdb


d = deque()

def SliceMean(item):
    global d
    old = d.popleft()            # pop oldest from left
    d.append(item)               # push newest in from right
    try :
        m = np.mean(d)
        return m
    except:
        print(d)
        return 0

def rm1d(y,M,truncate=False,wrap=False,zero=False,mirror=False):
    """
    rm1d: Finds the mean for the points in a sliding window (fixed size)
          as it is moved from left to right by one point at a time.
    Inputs:
          y: ndarray containing items for which a mean (in a sliding window) is
             to be calculated (N items)
          M: number of items in sliding window
    Options: The options follow the IDL 8.7 smooth options. Only one is valid.
             Boolean variables:
              truncate
              wrap
              zero
              mirror
            If no option is given, it fills with the y values at beginning and
            end.
    Otput:
          means: ndarray the same size of y

    Author: Adapted from a public code.
            Guiguesp @ Sao Paulo - 2018-03-01
            Introduced the map() function to speed the process


    """
    if not (M % 2) :
        M+=1        # must be odd

    N = y.shape[0]
    edge=True

    M2 = M//2

    if wrap:
        y=np.concatenate((y[-M2+1:],y,y[0:M2]))
    elif zero:
        y=np.concatenate((np.zeros(M2),y,np.zeros(M2)))
    elif mirror:
        y=np.concatenate((y[M2:0:-1],y,y[N-2:N-2-M2:-1]))
    elif truncate:
        y=np.concatenate((np.zeros(M2)+y[0],y,np.zeros(M2)+y[-1]))
    else:
        edge=False

    N = y.shape[0]

    # Load deque (d) with first window of seq
    global d
    d = deque(y[0:M])
    means = [np.mean(d)]             # contains mean of first window

    meansF = list((map(SliceMean,y[M:])))
    means  = means + meansF
    means  = np.asarray(means)

    if (not edge) :
        means = np.concatenate((y[0:M2+1],means,y[-M2+1:]))

    return means                     # return an ndarray
