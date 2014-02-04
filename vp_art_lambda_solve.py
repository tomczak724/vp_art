###  THIS FILE COMBINES A SUITE OF SCRIPTS FOR
###  CALCULATING THE WAVELENGTH SOLUTION FOR A
###  NIGHT.

import glob as gl
import numpy as np
import pylab as pl
import pyfits as pf
import decimal as deci
import vp_art_gaussfit as gf

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
    if type(n)==float or type(n)==np.float64 or type(n)==np.int64: n = int(n)
    if type(n)!=int: raise TypeError('Input needs to be an integer')

    s0 = str(n)
    if len(s0)>length: raise ValueError('Input is already larger than specified length')

    outer = ''
    for i in range( length - len(s0) ): outer+='0'

    return outer+s0


def num2sci(n, sigfig=2):
    if n==0: return '0.000e0'
    if abs(n)>=1: pow = int(np.log10(abs(n)))
    if abs(n)<1: pow = int(np.log10(abs(n)))-1
    p1 = str(round(n/10.0**pow,sigfig))
    while len(p1)<sigfig+2: p1+='0'
    p2 = str('e'+str(pow))
    return p1+p2



###  reads in the txt-files containing the wavelength
###  solutions info if it was already done.
def read_lamsoln_from_txt(file_prefix,fibers):
    for fib in fibers:
        fileinfo = np.loadtxt( gl.glob( file_prefix+filenum(fib.id,3)+'.*' )[0] )
        fib.lines_lam_x = fileinfo
    return fibers


###  This script finds the peaks associated
###  with the line-list for each fiber. 
def find_lines(comp, fibers, linesdat, params):

    xs = linesdat[:,1]
    min_dx = min( xs[1:] - xs[:-1] )

    dat = pf.getdata(comp)

    ### SORTING FIBERS BY ASCENDING ID
    fibs = []
    for i in range( 0 , params.FIBERS_TOTAL+1 ):
        for fib in fibers:
            if fib.id==i:
                fibs.append(fib)
                break
            
    fiber0 = fibs[0]

    #
    # FOR EXAMPLE:  xy = [ 0, [3,4,5,6,7] ]
    #

    for lam_x in linesdat:

        x = np.arange( int(lam_x[1])-int(min_dx)/2 , int(lam_x[1])+int(min_dx)/2 )
        f = dat[ fiber0.tracey ][ x ]

        fitx,fitf,msh = gf.gaussfit(x,f)

        line_x = fitx[ fitf.tolist().index(max(fitf)) ]
        fiber0.lines_lam_x.append( [ lam_x[0] , line_x ] )

    fib_cnt = -1
    for fib in fibs[1:]:
        fib_cnt += 1
        
        for lam_x in fibs[fib_cnt].lines_lam_x:

            x = np.arange( int(lam_x[1])-int(min_dx)/2 , int(lam_x[1])+int(min_dx)/2 )
            f = dat[ fib.tracey ][ x ]

            fitx,fitf,msh = gf.gaussfit(x,f)

            line_x = fitx[ fitf.tolist().index(max(fitf)) ]
            fib.lines_lam_x.append( [ lam_x[0] , line_x ] )

##  PLOTTING INITIAL SPECRAL FIT
        if fib_cnt==0:
            xaxis = [ i[0] for i in fib.xy ]
            flux = dat[ fib.tracey ][ xaxis ]
            fig = pl.figure(figsize=(15.025, 6.0))
            sp = fig.add_subplot(111)
            p = sp.plot( xaxis, flux, color='k', lw=2 )
            xmin,xmax = min(xaxis),max(xaxis)
            ymin,ymax = min(flux),max(flux)
            dx,dy = xmax-xmin, ymax-ymin
            ax = sp.axis([xmin-0.05*dx,xmax+0.05*dx,
                          ymin-0.05*dy,ymax+0.05*dy])

            for line in linesdat[:,1]: vline = sp.axvline(line, color='r', lw=1)

            stall = raw_input('\n\tHow does the initial guess look? (press enter to continue)\n')
            pl.close()


    for fib in fibers:
        outer1 = open('lambda_solve/solutions/lambda_soln_id'+filenum(fib.id,3)+'.txt','w')
        outer2 = open('lambda_solve/comp_spectra/comp_spectrum_id'+filenum(fib.id,3)+'.txt','w')
        for lam_x in fib.lines_lam_x: outer1.write( str(lam_x[0])+'  '+str(lam_x[1])+'\n' )
        for xn,fn in zip(xaxis,flux): outer2.write(str(xn)+' '+str(fn)+'\n')
        outer1.close()
        outer2.close()






def calc_wavlength_soln(comp, fiber, lines_lam_x, xaxis, polyorder):
    '''
    #  comp  ---------->  String of the traced comp image
    #
    #  fiber  --------->  Fiber class defined in vp_art_trace_fibers
    #
    #  lines_lam_x  --->  File containing the pixel locations of found peaks
    #                     with the corresponding wavelength that it should
    #                     represent (from lines.dat)
    #
    #  xaxis  --------->  Full range of pixels along the image to evaluate
    #                     the wavelength axis
    #
    #  polyorder  ----->  Order of the polynomial fit for the wavelength
    #                     solution specified in vp_art.param
    '''
    compdat = pf.getdata(comp)
    flux = compdat[ fiber.tracey ][ xaxis.astype(int) ]

    lams_true = lines_lam_x[:,0]
    xs_peak = lines_lam_x[:,1]

    poly = np.polyfit(xs_peak, lams_true, polyorder)

    lamaxis = np.zeros( len(xaxis) )
    for i in range(len(poly)): lamaxis += poly[-1-i]*xaxis**i
    xs_fit = np.interp( lams_true, lamaxis, xaxis )
    dx = xs_peak - xs_fit
    mn = num2sci( np.average(dx) )
    st = num2sci( np.std(dx) )
    print '\t'+filenum(fiber.id,3)+'  -----  '+mn+'  '+st

    fig = pl.figure(figsize=(15.025, 8.0))
    sp1 = fig.add_subplot(211)
    sp2 = fig.add_subplot(212)

    p = sp1.plot( xaxis, flux, color='k', lw=2)
    for x in xs_fit: vline = sp1.axvline( x, color='r', ls=':', lw=2 )
    xmin,xmax = min(xaxis),max(xaxis)
    ymin,ymax = min(flux),max(flux)
    xlen,ylen = xmax-xmin, ymax-ymin
    ax = sp1.axis([xmin-0.05*xlen,xmax+0.05*xlen,
                   ymin-0.05*ylen,ymax+0.05*ylen])


    p = sp2.scatter( xs_peak, dx, s=50, facecolor='r', linewidth=2 )
    zero = sp2.axhline(0, color='k', lw=2)
    meanline = sp2.axhline( np.float(mn), color='r', lw=1.5, ls='-', \
                           label='AVE = '+mn )
    stdline1 = sp2.axhline( np.float(mn)+np.float(st), \
                           color='r', lw=1.5, ls='--', label='STD = '+st )
    stdline2 = sp2.axhline( np.float(mn)-np.float(st), \
                           color='r', lw=1.5, ls='--' )

    ymin,ymax = np.float(mn)-3.5*np.float(st),np.float(mn)+3.5*np.float(st)
    ax = sp2.axis( [xmin-0.05*xlen,xmax+0.05*xlen , ymin,ymax] )
    sp2.set_xlabel('X_peaks   (pixel)')
    sp2.set_ylabel('X_peak - X_fit   (pixels)')
    sp2.legend(loc='upper left')

    fig.subplots_adjust(hspace=0)
    fig.savefig('lambda_solve/figures/fig_lambda_soln_id'+filenum(fiber.id,3)+'.png')
    pl.close()


    return lamaxis

