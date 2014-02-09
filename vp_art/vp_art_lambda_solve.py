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






def interactive_initial_guess(comp_im, trace_im, comp_names):

    comp = pf.getdata(comp_im)
    trc = pf.getdata(trace_im)

    # reading in reference spectrum
    ref_lam = np.loadtxt('/Users/tomczak/pydir/vp_art/data/spectrum_'+comp_names[0]+'.dat')[:,0]
    ref_flux = np.zeros(len(ref_lam))
    for comp0 in comp_names:
        compspecdat = np.loadtxt('/Users/tomczak/pydir/vp_art/data/spectrum_'+comp0+'.dat')
        ref_flux += compspecdat[:,1]
    ref_spec = np.array(zip(ref_lam, ref_flux))


    # reading in reference spectrum's line list
    lines = []
    for comp0 in comp_names:
        complinesdat = np.loadtxt('/Users/tomczak/pydir/vp_art/data/lines_'+comp0+'.dat')[:,0]
        lines += complinesdat.tolist()
    lines.sort()



    # extracting spectrum for the first fiber
    fibnums = pl.unique(trc)
    fibnums.sort()
    fibnums = fibnums[pl.find(fibnums > 0)]

    first_fiber_inds = pl.where(trc == fibnums[0])
    first_fiber_spec = pl.zeros(max(first_fiber_inds[1]) + 1)
    for y,x in zip(first_fiber_inds[0], first_fiber_inds[1]):
        first_fiber_spec[x] += comp[y][x]

    xaxis = range(max(first_fiber_inds[1])+1)

    class compspec:

        def __init__(self, xaxis, flux, ref_spec, ref_lines):
            self.xaxis, self.flux = xaxis, flux
            self.ref_spec = ref_spec
            self.ref_lines = ref_lines
            self.soln_data = []
            self.counter = 0

            self.fig = pl.figure(figsize=(15, 7))
            self.sp1 = self.fig.add_subplot(211)
            self.sp2 = self.fig.add_subplot(212)
            self.zooms = []

            self.sp1.plot(self.xaxis, self.flux/max(self.flux))
            self.sp1.set_ylim(-0.03, 1.2)
            self.sp1.set_title('ArcLamp Spectrum: x = identify line,  i = skip line,  z/a = zoom in/out,  r = reset')
            self.sp1.set_xlabel('Pixel')

            self.sp2.plot(self.ref_spec[:,0], self.ref_spec[:,1]/max(self.ref_spec[:,1]))
            self.sp2.set_ylim(-0.03, 1.2)
            x0, x1, y0, y1 = self.sp2.axis()
            self.ytext = y0 + (y1-y0)*0.85
            self.show_lines = [[self.sp2.axvline(self.ref_lines[0], color='r', lw=1, ls='--'),
                                self.sp2.annotate(str(int(self.ref_lines[0]+0.5)),
                                                  (self.ref_lines[0], self.ytext),
                                                  rotation='vertical', size=11)]]
            self.sp2.set_title('Reference Spectrum')
            self.sp2.set_xlabel('Wavelength (Angstroms)')

            self.fig.subplots_adjust(hspace=0.6, right=0.95)
            self.fig.canvas.mpl_connect('key_press_event', self.key)
            pl.show()

        def key(self, event):
            '''
            #  x = identify line
            #  i = skip to next line
            #  z = zoom in
            #  a = zoom out
            #  r = reset
            '''
            k = event.key
            x = event.xdata
            y = event.ydata
            sp = event.inaxes

            if k == 'x' and sp.get_subplotspec() == self.sp1.get_subplotspec():
                self.sp1.axvline(x, color='r', lw=1)
                self.show_lines[self.counter][0].set_ls('-')
                self.counter += 1
                if self.counter == len(self.ref_lines):
                    self.sp1.text(0.35, 0.35, "\n  Complete: close  \n     this window \n",
                                  transform=self.sp1.transAxes, size=22, bbox=dict(fc='w'))
                    self.sp2.text(0.35, 0.35, "\n  Complete: close  \n     this window \n",
                                  transform=self.sp2.transAxes, size=22, bbox=dict(fc='w'))
                    pl.savetxt('initial_lambda_soln.dat', self.soln_data)
                else:
                    self.soln_data.append([self.ref_lines[self.counter], x])
                    self.show_lines.append([self.sp2.axvline(self.ref_lines[self.counter],
                                                             color='r', lw=1, ls='--'),
                                            self.sp2.annotate(str(int(self.ref_lines[self.counter]+0.5)),
                                                              (self.ref_lines[self.counter], self.ytext),
                                                              rotation='vertical', size=11)])

            if k == 'i':
                self.show_lines[self.counter][0].set_lw(0)
                self.show_lines[self.counter][1].set_visible(False)
                self.counter += 1
                if self.counter == len(self.ref_lines):
                    self.sp1.text(0.35, 0.35, "\n  Complete: close  \n     this window \n",
                                  transform=self.sp1.transAxes, size=22, bbox=dict(fc='w'))
                    self.sp2.text(0.35, 0.35, "\n  Complete: close  \n     this window \n",
                                  transform=self.sp2.transAxes, size=22, bbox=dict(fc='w'))
                    pl.savetxt('initial_lambda_soln.dat', self.soln_data)
                else:
                    self.show_lines.append([self.sp2.axvline(self.ref_lines[self.counter],
                                                             color='r', lw=1, ls='--'),
                                            self.sp2.annotate(str(int(self.ref_lines[self.counter]+0.5)),
                                                              (self.ref_lines[self.counter], self.ytext),
                                                              rotation='vertical', size=11)])

            if k == 'z':
                self.zooms.append([self.sp1.axis()[:2], self.sp2.axis()[:2]])
                axis = sp.axis()
                dx = (axis[1] - axis[0])/10.
                sp.set_xlim(x-dx, x+dx)

            if k == 'a':
                if len(self.zooms):
                    self.sp1.set_xlim(self.zooms[-1][0])
                    self.sp2.set_xlim(self.zooms[-1][1])
                    self.zooms.pop(-1)

            if k == 'r':
                self.fig.clf()
                self.soln_data = []
                self.counter = 0

                self.sp1 = self.fig.add_subplot(211)
                self.sp2 = self.fig.add_subplot(212)
                self.zooms = []

                self.sp1.plot(self.xaxis, self.flux/max(self.flux))
                self.sp1.set_ylim(-0.03, 1.2)
                self.sp1.set_title('ArcLamp Spectrum: x = identify line,  i = skip line,  z/a = zoom in/out,  r = reset')
                self.sp1.set_xlabel('Pixel')

                self.sp2.plot(self.ref_spec[:,0], self.ref_spec[:,1]/max(self.ref_spec[:,1]))
                self.sp2.set_ylim(-0.03, 1.2)
                x0, x1, y0, y1 = self.sp2.axis()
                self.ytext = y0 + (y1-y0)*0.85
                self.show_lines = [[self.sp2.axvline(self.ref_lines[0], color='r', lw=1, ls='--'),
                                    self.sp2.annotate(str(int(self.ref_lines[0]+0.5)),
                                                      (self.ref_lines[0], self.ytext),
                                                      rotation='vertical', size=11)]]
                self.sp2.set_title('Reference Spectrum')
                self.sp2.set_xlabel('Wavelength (Angstroms)')

            pl.draw()

    interactive_go = compspec(xaxis, first_fiber_spec, ref_spec, lines)
