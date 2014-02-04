### WRITES DATA FILES FOR THE SPECTRA IN A SINGLE IMAGE

import numpy as np
import pyfits as pf

import vp_art_record_spectrum as spec


def filenum( n, length=4  ):
    """
    #  DESCRIPTION:
    #    Takes an integer and returns a string
    #    of it of the specified length.
    #
    #  EXAMPLES:
    #    >>> mypy.filenum(34)
    #    >>>    '0034'
    #    >>>
    #    >>> mypy.filenum(8,6)
    #    >>>    '000008'
    #    >>>
    """
    if type(n)==float or type(n)==np.float64: n = int(n)
    if type(n)!=int: raise TypeError('Input needs to be an integer')

    s0 = str(n)
    if len(s0)>length: raise ValueError('Input is already larger than specified length')

    outer = ''
    for i in range( length - len(s0) ): outer+='0'

    return outer+s0




def write_spectra(data_im, trace_im, fibers):

    dat = pf.getdata(data_im)
    trc = pf.getdata(trace_im)

    outer_files = []

    for fib in fibers:

        wid = len(fib.xy[0][1])
        ys = np.arange( fib.tracey-wid/2 , fib.tracey+wid/2+1 )
        ys = np.array( [ int(i) for i in ys ] )

        if 0<fib.id<1000: filename = 'spec_'+data_im[:-5]+'_fib'+spec.filenum(fib.id,3)+'.dat'
        else: raise ValueError('Fiber ID out of range, must be <1000')

        outer = open( filename, 'w' )
        outer_files.append( filename )

        for i in range(len(fib.xy)):

            xy = fib.xy[i]
            lam = fib.lamaxis[i]

            weights = trc[ys,xy[0]]
            fluxes  = dat[ys,xy[0]]

            wgt_ave = np.average( fluxes, weights=weights )

            outer.write( str(lam)+'  '+str(wgt_ave)+'\n' )
        outer.close()
    return outer_files
        
