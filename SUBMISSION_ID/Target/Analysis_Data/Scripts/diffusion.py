from numpy import *
import sys
from numpy.fft import fft, ifft, fftshift

def autocorr(xi,yi):
    """FFT based autocorrelation, x is a nx1 array. The code can be easily extended for an m by n matrix"""

    length = xi.shape[0]
    x =	zeros((1,length))
    x[0,:] = xi[:]

    y = zeros((1,length))
    y[0,:] = yi[:]

    fftx = fft(x, n=(length*2-1), axis=1) #padding with zeros affects the final outcome
    ffty = fft(y, n=(length*2-1), axis=1)

    corr_xy = ifft(fftx * conjugate(ffty), axis=1)
    corr_xy = real(fftshift(corr_xy, axes=1)) #assumes no complex part, remove real for complex arrays
    
    corr_yx = ifft(ffty * conjugate(fftx), axis=1)
    corr_yx = real(fftshift(corr_yx, axes=1))

    corr = 0.5 * (corr_xy[:,length:] + corr_yx[:,length:]) / range(1,length)[::-1]
    return corr[0,:]
    
def diffusion(vi, vj, dt, usedFraction = .5):
    """
    calculate the diffusion coefecient for a pair of time series vi, vj
    param vi: an n * nt array of time series
    param vj: an n * nt array of time series
    dt: timestep
    usedFraction: fraction of data points to use
    sampled at interval dt
    """
    
    #number of points to use for integration
    npts = int(usedFraction * vi.shape[1])

    # the velocity (auto)correlation function for a time average.
    corr = zeros(npts-1)

    for i in arange(vi.shape[0]):
        corr += Autocorr(vi[i,:npts], vj[i,:npts])
        
    #average correlation func        
    correlation = corr / float(vi.shape[0])

    return trapz(correlation, dx=dt)
