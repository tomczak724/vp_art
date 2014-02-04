###  THIS SCRIPT USES THE FLAT-STACK TO TRACE OUT THE
###  FIBERS. EACH FIBER WILL BE 5 PIXELS WIDE AND
###  "STRAIGHTENED" ALONG THE X-AXIS.
###
###  A GAUSSIAN IS FIT TO THE FIBER AT EACH X COORDINATE
###  AND 

import os
import glob as gl
import pylab as pl
import pyfits as pf
import shutil as shu
from scipy import optimize as opt
import vp_art_gaussfit as gf

import vp_art_record_spectrum as spec
import os
import glob as gl
import pylab as pl
import pyfits as pf
import shutil as shu
from scipy import optimize as opt
import vp_art_gaussfit as gf

import vp_art_record_spectrum as spec


def parabola(x, a, b, c): return a*x**2 + b*x + c


###  id = the absolute ID for the fiber (ie if 121 is dead it will be excluded)
###  xy  = xaxis of the fibers and the corresponding pixels along the y axis
class fiber:
    def __init__(self, id, x0, y0):
        self.id = id
        self.xy = [ [x0,y0] ]

    ###  gaussfit parameters rescaling array
        self.h = []
        self.mu = []
        self.sig = []
    ###  format for "rescale"...
    ###  [ [ 0.99, 1.01, 0.98, 1.00, 0.99, 1.02 ], [...], ... ]
    ###  where arrays are normalized scalings for
    ###  each width along the x-axis
        self.rescale = []

        self.tracey = 0
        self.lines_lam_x = []
        self.lamaxis = 0


###  Takes a single linear array of the peaks and returns an
###  array of arrays of the peaks.
###    EXAMPLE:  [1,2,3,7,8,9] --> [ [1,2,3], [7,8,9] ]
###  If dx between adjacent values is >1, starts a new peak
def sep_peaks( inx, iny ):
    outer_x = []
    outer_y = []
    tmp_x = [ inx[0] ]
    for i in inx[1:]:
        if i-tmp_x[-1]<1.1:
            tmp_x.append(i)
        else:
            outer_x.append( pl.array(tmp_x) )
            outer_y.append( iny[pl.array(tmp_x)] )
            tmp_x = [i]
    #outer_x.append(tmp_x)
    outer_x.append( pl.array(tmp_x) )
    outer_y.append( iny[pl.array(tmp_x)] )
    return outer_x , outer_y




###  reads in the txt-files containing the
###  tracing info if it was already done.
def read_fibers_from_txt(files, fiber_ys):
    fibers  = []
    for txt in files:
        dat = pl.loadtxt(txt)
        width = (len(dat[0])-2)/2
        id = int(dat[0][0])
        x0 = dat[0][1]
        y0 = dat[0][ 2 : 2+width ]
        rescale0 = dat[0][ 2+width : ]
        fib = fiber( id , x0 , y0  )
        fib.rescale.append( rescale0 )
        for line in dat[1:]:
            fib.xy.append( [ line[1], line[ 2 : 2+width ] ] )
            fib.rescale.append( line[ 2+width : ] )
        for id_y in pl.loadtxt(fiber_ys):
            if id_y[0]==id: fib.tracey=int(id_y[1])
        fibers.append(fib)
    return fibers



###  input "flatim" is a string
###  This script traces out the fibers and returns
###  and array of fiber-structures containing the
###  pixel-corrdinate information.
def trace_fibers(flatim, params):

    trace_im1 = pf.getdata(flatim)*0

    imdat, imhead = pf.getdata(flatim), pf.getheader(flatim)

    ###  separating fibers in first column and assigning fibers ids
    print '\n\tSEARCHING FOR FIBERS BETWEEN x=0 & y=[' +str(params.LO_BUFFER) +':'+str(len(imdat[:,0])-params.HI_BUFFER)+']'
    fiber_peaks_pix = pl.find( imdat[:,0] > pl.median(imdat[:,0]) )
    fiber_peaks_pix,fiber_peaks_flx = sep_peaks( fiber_peaks_pix, imdat[:,0] )
    if pl.median(fiber_peaks_pix[0])<=params.LO_BUFFER:
        fiber_peaks_pix = fiber_peaks_pix[1:]
        fiber_peaks_flx = fiber_peaks_flx[1:]
    if (len(imdat[:,0])-pl.median(fiber_peaks_pix[-1]))<=params.HI_BUFFER:
        fiber_peaks_pix = fiber_peaks_pix[:-1]
        fiber_peaks_flx = fiber_peaks_flx[:-1]
    print '\t  --> FOUND ', len(fiber_peaks_pix), ' FIBER PEAKS'

    ###  creating array for fibers
    fibers0 = []
    id_cnt = 1
    for f in range(len(fiber_peaks_pix)):
        while params.FIBERS_EXCLUDE.tolist().count(id_cnt)==1: id_cnt+=1

        fibx,fiby = fiber_peaks_pix[f],fiber_peaks_flx[f]
        peakx = fibx[ fiby.tolist().index(max(fiby)) ]
        yrange = pl.arange(  peakx-params.FIBER_WIDTH/2  ,  peakx+params.FIBER_WIDTH/2+1  )

        fibers0.append( fiber(id_cnt, 0,     yrange     ))
        id_cnt+=1

##  TRACING FIBERS ALONG X-AXIS INDIVIDUALLY
    for fib in fibers0:
        for x in range(1,len(imdat)):
##  FIRST, TAKE THE FLUXES IN THE PIXELS AT x
##  THAT ENCOMPASSED THE PEAK AT x-1
            fluxes = imdat[ fib.xy[-1][1] , x ]
##  NEXT, FIND THE VERTICAL SHIFT TO CENTER ON
##  THE PEAK FLUX AT x. MAXIMUM SHIFT IS DETERMINED
##  FROM THE FIBER_WIDTH PARAMETER.
            deltay = range( -len(fluxes)/2+1 , len(fluxes)/2+1 )[ fluxes.tolist().index(max(fluxes)) ]

##  RECORD THE NEW Y-PIXELS THAT ARE CENTERD ON
##  THE FIBER AT x.
            fib.xy.append( [ x, fib.xy[-1][1]+deltay ] )

##  FLAG PIXELS FOR FIBER IN FIRST-PASS TRACE IMAGE
            trace_im1[fib.xy[-1][1],x] = fib.id



    trc0 = 'trace_pass1.fits'
    print '\n\tWRITING INITIAL TRACING TO ', trc0
    try: pf.writeto(trc0, trace_im1, header=imhead)
    except:
        os.remove(trc0)
        pf.writeto(trc0, trace_im1, header=imhead)
        
    return fibers0



###  This routine runs through the given width of a given fiber at a given x-coord...
###      (1) Fits a gaussian
###      (2) Shifts the fit to the peak of the flat data
###      (3) Reshapes the data to fit the gaussian
###      (4) Rescales to conserve flux
def calc_trace_align(flatim, biasim, fibers1, fiber_folder, params):

    imdat, imhead = pf.getdata(flatim), pf.getheader(flatim)
    pass2 = pf.getdata(flatim)*0
    weight_im = pf.getdata(flatim)*0
    pass3 = pf.getdata(flatim)*0
    fibers_y = open('fibertrace_ys.txt','w')

    print '\n\tBEGINNING TRACE-ALIGNMENT'
    regfile = open('trace.reg','w')
    for fib in fibers1:

        fib.tracey = (len(fib.xy[0][1])+1)*fib.id+len(fib.xy[0][1])/2
        fibers_y.write( str(fib.id)+'  '+str(fib.tracey)+'\n' )

        if 0<fib.id<1000: tracename = 'fibertrace_id'+spec.filenum(fib.id,3)
        else: raise ValueError('Fiber ID out of range, must be <1000')

        tracename += str(fib.id)+'.txt'
        tracetxt = open(tracename,'w')
        tracetxt.write('#  id  x_ind  [y_inds]  [rescale]\n')


# update 2014
        regx, regy = [],[]
        for xy in fib.xy:
            regx.append(xy[0])
            regy.append( (xy[1][0]+xy[1][-1])/2. )
        regx = pl.array(regx)*1.0
        regy = pl.array(regy)*1.0
        fit = opt.curve_fit( parabola, regx, regy, p0=[(2.0e-7), -(2.5e-3), 1000.0], sigma=pl.zeros(len(regx))+1 )

        for x in range(0, len(regx)-5, 5):
            regfile.write( 'line( '+str(x+1)+' , '+str(parabola(x,*fit[0])+1)+' , '+str(x+5+1)+' , '+str(parabola(x+5,*fit[0])+1)+' ) # line=0 0\n' )
# update 2014


        for xy in fib.xy:

            tracetxt.write('  '+str(fib.id)+'  '+str(xy[0]) )
            for y in xy[1]: tracetxt.write('  '+str(y))

            col = imdat[:,xy[0]]

            dy = (params.FIBER_WIDTH+1)/2
            inx = xy[1]
            iny = col[inx]

            fitx, fity, msh = gf.gaussfit(inx,iny)

            fib.mu.append(msh[0])
            fib.sig.append(msh[1])
            fib.h.append(msh[2])

            peaks_dif = inx[dy-1] - msh[0]
            fityvals = msh[2]*pl.exp( - ( ((inx-msh[0]-peaks_dif)/msh[1])**2 ) )
            rescale = fityvals/iny
            rescalenorm = rescale * pl.trapz(iny,inx) / pl.trapz(iny*rescale,inx)
            fib.rescale.append( rescalenorm )

            for y in range( len(iny) ):
                pass2[(len(iny)+1)*fib.id+y][xy[0]] = iny[y]*rescalenorm[y]
                weight_im[(len(iny)+1)*fib.id+y][xy[0]] = iny[y] / max(iny)
                pass3[(len(iny)+1)*fib.id+y][xy[0]] = fib.id
                tracetxt.write( '  '+str(rescalenorm[y]) )
            tracetxt.write('\n')

        tracetxt.close()
        shu.move(tracename,fiber_folder)
        if fib.id%(len(fibers1)/10)==0: print '\t  --> TRACE-ALIGNED ', fib.id, ' FIBERS OF ', fibers1[-1].id
    regfile.close()
    fibers_y.close()
    shu.move('fibertrace_ys.txt',fiber_folder)

    try: pf.writeto('trace_pass2.fits', pass2, header=imhead)
    except:
        os.remove('trace_pass2.fits')
        pf.writeto('trace_pass2.fits', pass2, header=imhead)

    try: pf.writeto('trace_weights.fits', weight_im, header=imhead)
    except:
        os.remove('trace_weights.fits')
        pf.writeto('trace_weights.fits', weight_im, header=imhead)

    try: pf.writeto('trace_pass3.fits', pass3, header=imhead)
    except:
        os.remove('trace_pass3.fits')
        pf.writeto('trace_pass3.fits', pass3, header=imhead)

    return fibers1



def do_trace_align(im,fibers2,params):


    imdat, imhead = pf.getdata(im), pf.getheader(im)

    outname = im[:-5]+'_t.fits'
    #outer = pl.array( [ pl.arange(len(dat[0]))*0 for i in range(params.FIBER_WIDTH*len(fibers2)) ] )
    outer = pf.getdata(im)*0

    # DONT DO IT IF THE FILE EXISTS
    globber = gl.glob(outname)
    if len(globber)==1: return outname

    for fib in fibers2:

        fib_cnt = 1

        for f in range(len(fib.xy)):

            xy = fib.xy[f]
            rescale = fib.rescale[f]
            #rescaled_vals = dat[ xy[1] , xy[0] ] * rescale

            #            print 'xy = ',xy
            #            print 'rescale = ',rescale
            #            print 'rescaled = ',rescaled_vals

            for y in range( len(rescale) ):
                rescaled_val = imdat[ xy[1][y] ][ xy[0] ] ###* rescale[y]
                outer[(len(rescale)+1)*fib.id+y][xy[0]] = rescaled_val

        fib_cnt += 1

    pf.writeto( outname, outer, header=imhead )
    return outname



