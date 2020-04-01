import scipy
from numpy import *
from scipy.optimize import leastsq
from scipy.optimize import fmin_slsqp
import sys
from numpy.fft import fft, ifft, fftshift

DIFFDEBUG = False 

if DIFFDEBUG:
    import pylab
    import pdb
    

def d(vi, vj, dt):
    """
    calculate the diffusion coefecient for a pair of time series vi, vj
    @param vi: an n * nt array of time series
    @param vj: an n * nt array of time series
    sampled at interval dt
    """
    
    # fraction of time series to use for fitting correlation func
    usedFraction = .25
    
    #number of points to use for integration
    npts = int(floor(vi.shape[1]*usedFraction))
    
    #pylab.cla()
    
    # the velocity (auto)correlation function for a time average.
    correlation = None
    for i in arange(vi.shape[0]):
        corr = array(map(lambda t: trapz(roll(vi[i], t) * vj[i], dx=dt), arange(npts))) * 1./(npts * dt)
        #pylab.plot(corr)
        if correlation == None:
            correlation = corr
        else:
            correlation = correlation + corr
    
    #average correlation func        
    correlation = correlation / float(vi.shape[0])        
    
    #pylab.plot(correlation, 'b-', linewidth=8)
    
    x = scipy.linspace(0, dt*correlation.size, correlation.size)
    
    print('x,corr: \n', x.shape, correlation.shape)
    
    # linear fit A and B, used as guess to least squares fit
    B,A=polyfit(x,correlation,1)
    
    print('polyfit: ', A, B)
    
    ## Error function
    # e = lambda v, t, y: (v[0] * exp(-v[2] * t) + v[1] * t * exp(-v[3] * t) - y)
    
    def e(v,t,y):
        return (v[0] * exp(-v[2] * t) + v[1] * t * exp(-v[3] * t) - y)
    
    
    ## Initial parameter value
    v0 = [A,B,1,1]
          
    ## Fitting
    v, success = leastsq(e, v0, args=(x,correlation), maxfev=10000)
    
    #pylab.plot(map(lambda t: v[0] * exp(-v[2] * t) + v[1] * t * exp(-v[3] * t), x), 'ro')
    
    print(v,success)
    
    return coorelation_integral(v)

def d2(vi, vj, dt):
    """
    calculate the diffusion coefecient for a pair of time series vi, vj
    @param vi: an n * nt array of time series
    @param vj: an n * nt array of time series
    sampled at interval dt
    """
    
    # fraction of time series to use for fitting correlation func
    usedFraction = .25
    
    #number of points to use for integration
    npts = int(usedFraction * vi.shape[1])
    #npts = int(floor(vi.shape[1]/2))
    
    #pylab.cla()
    
    
    
    # the velocity (auto)correlation function for a time average.
    correlation = None
    for i in arange(vi.shape[0]):
        corr = array(map(lambda t: trapz(vi[i,t:t+npts]*vj[i,0:npts], dx=dt), arange(vi.shape[1]-npts))) * 1./(npts * dt)
        #corr = array(map(lambda t: trapz(roll(vi[i], t) * vj[i], dx=dt), arange(npts))) * 1./(npts * dt)
        pylab.plot(corr)
        if correlation == None:
            correlation = corr
        else:
            correlation = correlation + corr
    
    #average correlation func        
    correlation = correlation / float(vi.shape[0])        
    
    

    pylab.plot(correlation, 'b-', linewidth=8)
    
    x = scipy.linspace(0, dt*correlation.size, correlation.size)
    
    print('x,corr: \n', x.shape, correlation.shape)
    
    # linear fit A and B, used as guess to least squares fit
    B,A=polyfit(x,correlation,1)
    
    print('polyfit: ', A, B)
    
    ## Error function
    # e = lambda v, t, y: (v[0] * exp(-v[2] * t) + v[1] * t * exp(-v[3] * t) - y)
    
    def e(v,t,y):
        #print('v: ', v)
        #print('t: ', t)
        #print('y: ', y)
        return (v[0] * exp(-v[2] * t) + v[1] * t * exp(-v[3] * t) - y)
    

       
       
    
    
    ## Initial parameter value
    v0 = [A,B,1,1]
          
    ## Fitting
    v, success = leastsq(e, v0, args=(x,correlation), maxfev=10000)
    
    pylab.plot(map(lambda t: v[0] * exp(-v[2] * t) + v[1] * t * exp(-v[3] * t), x), 'ro')
    
    print(v,success)
    
    return coorelation_integral(v)

def d3(vi, vj, dt):
    """
    calculate the diffusion coefecient for a pair of time series vi, vj
    @param vi: an n * nt array of time series
    @param vj: an n * nt array of time series
    sampled at interval dt
    """
    
    # fraction of time series to use for fitting correlation func
    usedFraction = .25
    
    #number of points to use for integration
    npts = int(usedFraction * vi.shape[1])
    
    pylab.cla()
    
    # the velocity (auto)correlation function for a time average.
    correlation = None
    for i in arange(vi.shape[0]):
        corr = array(map(lambda t: trapz(vi[i,t:] * vj[i,:vi.shape[1]-t], dx=dt) / ((vi.shape[1]-t)*dt), arange(npts))) 
        pylab.plot(corr)
        if correlation == None:
            correlation = corr
        else:
            correlation = correlation + corr
    
    #average correlation func        
    correlation = correlation / float(vi.shape[0])        
    
    pylab.plot(correlation, 'b-', linewidth=8)
    
    x = scipy.linspace(0, dt*correlation.size, correlation.size)
    
    print('x,corr: \n', x.shape, correlation.shape)
    
    # linear fit A and B, used as guess to least squares fit
    B,A=polyfit(x,correlation,1)
    
    print('polyfit: ', A, B)
    
    ## Error function
    # e = lambda v, t, y: (v[0] * exp(-v[2] * t) + v[1] * t * exp(-v[3] * t) - y)
    
    def e(v,t,y):
        return (v[0] * exp(-v[2] * t) + v[1] * t * exp(-v[3] * t) - y)
    
    
    ## Initial parameter value
    v0 = [A,B,1,1]
          
    ## Fitting
    v, success = leastsq(e, v0, args=(x,correlation), maxfev=10000)
    
    pylab.plot(map(lambda t: v[0] * exp(-v[2] * t) + v[1] * t * exp(-v[3] * t), x), 'ro')
    
    print(v,success)
    
    return coorelation_integral(v)

def d4(vi, vj, dt):
    """
    calculate the diffusion coefecient for a pair of time series vi, vj
    @param vi: an n * nt array of time series
    @param vj: an n * nt array of time series
    sampled at interval dt
    """
    
    # fraction of time series to use for fitting correlation func
    usedFraction = .5
    
    #number of points to use for integration
    npts = int(usedFraction * vi.shape[1])
    
    if DIFFDEBUG:
        pylab.cla()
    
    # the velocity (auto)correlation function for a time average.
    correlation = None
    for i in arange(vi.shape[0]):
        corrf = array([trapz(vi[i,t:] * vj[i,:vi.shape[1]-t], dx=dt) for t in arange(npts)])
        corrb = array([trapz(vj[i,t:] * vi[i,:vi.shape[1]-t], dx=dt) for t in arange(npts)])
        den = linspace(2.0 * dt * vi.shape[1], 2.0 * dt * (vi.shape[1] - npts + 1), npts)
        corr = (corrf + corrb) / den
        
        if DIFFDEBUG:
            pylab.plot(scipy.linspace(0, dt*corr.size, corr.size), corr)
            
        if correlation == None:
            correlation = corr
        else:
            correlation = correlation + corr
    
    #average correlation func        
    correlation = correlation / float(vi.shape[0])     
    x = scipy.linspace(0, dt*correlation.size, correlation.size)   
    
    scale = exp(-50.0*linspace(0,1.0,x.size)**2)
    rescaled_correlation = scale * correlation
    
    if DIFFDEBUG:
        pylab.plot(x, correlation, 'b-', linewidth=8)
        pylab.plot(x, rescaled_correlation)
        pylab.plot(x, correlation[0] * scale)
    
    
    
    ## Initial parameter value
    v0 = [0, correlation[0], 1]
    
    #mfloat = sys.float_info.max
    #epsilon = sys.float_info.epsilon
    
    #if DIFFDEBUG:
    #    iprint = 3
    #else:
    #    iprint = 0
        
    def leastsq_err(p,x,y):
        # idea here is that the constrained minimization function fmin_slsqp
        # does not work correctly, it does not obey bounds, or at least I can't 
        # get it to obey bounds. Furthermore, it does not optimize correctly as
        # the leastsq seems to. So, we fake the constraint of only postive exponents
        # by having the err function return at least 2 times the data function if 
        # the exponent is negative. The idea is, that an all zero fit function will
        # always be better than a negative exponent. 
        if p[2] <= sys.float_info.epsilon:
            #print('warning, negative exponent {}'.format(p[2]))
            return y * (p[2] - 2.0)
        #print('v: {}'.format(p))
        return y - correlation_func(p,x)
    
    v, success = leastsq(leastsq_err, v0, args=(x,rescaled_correlation), maxfev=10000)
    
          
    #v = fmin_slsqp(obj, v0, 
    #                     bounds=[(-mfloat,mfloat), (-mfloat, mfloat), (epsilon, mfloat)], 
    #                     args=(x,correlation),
    #                     iprint=iprint,
    #                     acc=1e-2)
    
    if DIFFDEBUG:
        pylab.plot(x, map(lambda t: correlation_func(v,t), x), 'ro')
        print('v: '.format(v))
        print('verify: {}'.format(trapz(map(lambda t: correlation_func(v,t), scipy.linspace(0, dt*10000, 10000)), dx=dt)))
        
    return coorelation_integral(v)

def Autocorr(xi,yi):
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
    
def d5(vi, vj, dt):
    """
    calculate the diffusion coefecient for a pair of time series vi, vj
    @param vi: an n * nt array of time series
    @param vj: an n * nt array of time series
    sampled at interval dt
    """
    
    # fraction of time series to use for fitting correlation func
    usedFraction = .5
    
    #number of points to use for integration
    npts = 5 #int(usedFraction * vi.shape[1])
    
    if DIFFDEBUG:
        pylab.cla()
    
    # the velocity (auto)correlation function for a time average.
    corr = zeros(npts-1)

    for i in arange(vi.shape[0]):
        corr += Autocorr(vi[i,:npts], vj[i,:npts])
        
    #average correlation func        
    correlation = corr / float(vi.shape[0])     
    #index = where(correlation <= .0)[0]
    #if index.shape[0] < 1:
	#index = correlation.shape[0]
    #else:
	#index = index[0]

    return trapz(correlation, dx=dt)

def correlation_func(v,t):
    return (v[0] * t + v[1]) * exp(-v[2] * t) 

def correlation3(v,t):
    return v[0] * exp(-v[3] * t) + v[1] * t * exp(-v[4] * t) + v[2] * (t**2) * exp(-v[5] * (t**2))  

def coorelation_integral(v):
    return (v[0] + v[2]*v[1]) / (v[2]*v[2])

def correlation_deriv(v,t):
    dfdA        = exp(-v[2] * t) 
    dfdB        = t * exp(-v[3] * t)
    dfdalpha    = -v[0] * t * exp(-v[2] * t) 
    dfdbeta     = -v[1] * (t**2) * exp(-v[3] * t)
    return array([dfdA, dfdB, dfdalpha, dfdbeta], 'f')


def dmatrix(v, dt, box):
    """
    @param src: a nframe * nsubsystem * nop * 3 array
    @return: the diffusion tensor
    """
    dtensor = zeros([v.shape[2], v.shape[3], 3, v.shape[2], v.shape[3], 3], 'f')
    
    #for ri in arange(dtensor.shape[0]):
    #    for rj in arange(dtensor.shape[1]):
    #        for rk in arange(dtensor.shape[2]):
    #            for ci in arange(ri, dtensor.shape[3]):
    #                for cj in arange(rj, dtensor.shape[4]):
    #                    for ck in arange(rk, dtensor.shape[5]):
    #                        dtensor[ri,rj,rk,ci,cj,ck] = d4(v[:,:,ri,rj,rk],v[:,:,ci,cj,ck],dt)
    #                        print('{},{},{},{},{},{}={}'.format(ri,rj,rk,ci,cj,ck, dtensor[ri,rj,rk,ci,cj,ck]))
    
    for ri in arange(dtensor.shape[0]):
        for rj in arange(dtensor.shape[1]):
            for rk in arange(dtensor.shape[2]):
                dtensor[ri,rj,rk,ri,rj,rk] = d5(v[:,:,ri,rj,rk],v[:,:,ri,rj,rk],dt * 0.001 * box[rk] * box[rk])

    size=dtensor.shape[0]*dtensor.shape[1]*dtensor.shape[2]
    return reshape(dtensor,(size,size))

def test(v, dt, i,j,k):
    return d4(v[:,:,i,j,k],v[:,:,i,j,k],dt)
    
    

    
