from scipy.special._ufuncs import gammainc, gamma
import numpy as np
from scipy.optimize import fminbound
__author__ = 'Evgeniya Predybaylo'


def wavelet(Y, dt, pad=0, dj=-1, s0=-1, J1=-1, mother=-1, param=-1):
    
    n1 = len(Y)

    if s0 == -1:
        s0 = 2 * dt
    if dj == -1:
        dj = 1. / 4.
    if J1 == -1:
        J1 = np.fix((np.log(n1 * dt / s0) / np.log(2)) / dj)
    if mother == -1:
        mother = 'MORLET'

    #....construct time series to analyze, pad if necessary
    x = Y - np.mean(Y)
    if pad == 1:
        base2 = np.fix(np.log(n1) / np.log(2) + 0.4999)  # power of 2 nearest to N
        x = np.concatenate((x, np.zeros((2 ** (base2 + 1) - n1).astype(np.int64))))

    n = len(x)

    #....construct wavenumber array used in transform [Eqn(5)]
    kplus = np.arange(1, np.fix(n / 2 + 1))
    kplus = (kplus * 2 * np.pi / (n * dt))
    kminus = (-(kplus[0:-1])[::-1])
    k = np.concatenate(([0.], kplus, kminus))

    #....compute FFT of the (padded) time series
    f = np.fft.fft(x)  # [Eqn(3)]

    #....construct SCALE array & empty PERIOD & WAVE arrays
    j = np.arange(0, J1+1)
    scale = s0 * 2. ** (j * dj)
    wave = np.zeros(shape=(int(J1 + 1), n), dtype=complex)  # define the wavelet array

    # loop through all scales and compute transform
    for a1 in range(0, int(J1+1)):
        daughter, fourier_factor, coi, dofmin = wave_bases(mother, k, scale[a1], param)
        wave[a1, :] = np.fft.ifft(f * daughter)  # wavelet transform[Eqn(4)]

    period = fourier_factor * scale  # [Table(1)]
    coi = coi * dt * np.concatenate((np.insert(np.arange(int((n1 + 1) / 2 - 1)), [0], [1E-5]),
        np.insert(np.flipud(np.arange(0, n1 / 2 - 1)), [-1], [1E-5])))  # COI [Sec.3g]
    wave = wave[:, :n1]  # get rid of padding before returning

    return wave, period, scale, coi
