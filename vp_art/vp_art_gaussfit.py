#
#  DESCRIPTION:
#    This is a crude script for fitting a gaussian function
#    to input data arrays.
#
#  EXAMPLE:
#    >>>  xout, yout = mypy.gaussfit( xin, yin )
#
#  INPUTS:
#    (1) x --- data array for x values
#    (2) y --- data array for y values
#
#  OUTPUTS:
#    (1) fitx --- fit array x values
#    (2) fity --- fit array y values
#
#  NOTES:
#    (1) Like I said, "crude". Try to make sure the inputs
#        are not too large and centered near the peak
#    (2) The output fitx has 10x the sampling as the input x

import numpy as np
from scipy import optimize as opt

class Parameter:
    def __init__(self, value): self.value = value
    def set(self, value): self.value = value
    def __call__(self): return self.value

def fitter(function, parameters, y, x = None):
    def f(params):
        i=0
        for p in parameters:
            p.set(params[i])
            i+=1
        return y - function(x)
    if x is None: x = arange(y.shape[0])
    p = [param() for param in parameters]
    opt.leastsq(f, p)

def gaussfit( x, y ):

# GUESSING THE GAUSSIAN PARAMETERS
    mu0 = Parameter( (x[-1]+x[0])/2.0 )
    sig0 = Parameter( (x[-1]-x[0])/2.0 )
    h0 = Parameter( max(y)-min(y) )

    def gau(x): return h0() * np.exp(-((x-mu0())/sig0())**2)

# PERFORMING FIT AND PRODUCING OUTPUTS
    fitter( gau, [mu0, sig0, h0], y, x )
    fitx = np.linspace( x[0], x[-1], 10*len(x) )
    fity = gau( fitx )

    return fitx, fity, [mu0.value, sig0.value, h0.value]
    
    
    
